using System.Globalization;
using Gurobi;
/*
Mathematical programming model for masonry stability check.
Primary variables: fx[i], fy[i] (global forces at vertex i).
Constraints: Block equilibrium, Vertex No-Tension, Vertex Friction, Face Strength/Eccentricity.
Goal: Maximize load factor lambda.

Input Files:
1. Equilibrium Matrix File (*.txt): Contains MatrixA and base load vector B.
   Format:
   matrix_A_size
   [n = rows = num_blocks * 3]
   [p = cols = num_vertices * 2]
   matrix_A_values
   [n*p values...]
   vector_B_values
   [n values...]
2. Face Geometry File (*.txt):
   Format (per line): # faceId, v1, v2, length, thickness ex: 0, 10, 11, 2.0, 0.15
       # lines are disregarded

*/
namespace InterlockingMasonryLocalForces
{
    // -----------------------------------------------------------------
    //  Data Classes (ProblemData, GeometryModel, Face, ContactPoint)
    // -----------------------------------------------------------------

    /// Global problem parameters and data arrays.
    public class ProblemData
    {
        // Matrix-based equilibrium: A * f_all = B * lambda
        public double[,] MatrixA { get; set; }  // Size (n x p)
        public double[] VectorB { get; set; }   // Size (n)

        public int NumRows => (MatrixA == null) ? 0 : MatrixA.GetLength(0);   // n
        public int NumCols => (MatrixA == null) ? 0 : MatrixA.GetLength(1);   // p



        // friction, compressive strength, etc.
        public double Mu { get; set; } = 0.4;    // friction
        public double SigmaC { get; set; } = 20000;  // compressive strength
        // For partial contact: thickness, etc., read from Face definitions.
        public double Cohesion { get; set; } = 0.0;
    }

    /// Represents a single contact point (vertex).
    public class ContactPoint
    {
        public int Id { get; set; }
    }

    /// Represents a 2D contact face, assumed 2 vertexes, but adaptable
    public class Face
    {
        public int Id { get; set; }
        public List<int> VertexIds { get; set; } = new List<int>();
        public double Length { get; set; }    // contact length (height) in 2D
        public double Thickness { get; set; }
    }

    /// A container for all geometry in the problem (faces, vertices, etc.).
    public class GeometryModel
    {
        public Dictionary<int, ContactPoint> Vertices { get; } = new Dictionary<int, ContactPoint>();
        public Dictionary<int, Face> Faces { get; } = new Dictionary<int, Face>();
    }

    /// Simple struct to identify one face-vertex pair
    public struct FaceVertexPair
    {
        public int FaceId;
        public int VertexId;
        public FaceVertexPair(int f, int v)
        {
            FaceId = f;
            VertexId = v;
        }
    }

    // -----------------------------------------------------------------
    // 2) TheLocalOptimizer class
    // -----------------------------------------------------------------
 
    public class LocalOptimizer
    {
            // We'll store the list of face-vertex pairs, in order
            private List<FaceVertexPair> faceVertexPairs;

            // Gurobi variable array (size = 2 * faceVertexPairs.Count):
            // the even indices are f_n, the odd are f_t, for that pair.
            private GRBVar[] fAll;
            private GRBVar lambda;


            /// The main solve routine
            public void SolveProblem(GeometryModel geometry, ProblemData data)
            {
                // 1) Collect face-vertex pairs in a consistent order
                faceVertexPairs = new List<FaceVertexPair>();
                foreach (var fkvp in geometry.Faces)
                {
                    int faceId = fkvp.Key;
                    var face = fkvp.Value;
                    foreach (int vId in face.VertexIds)
                    {
                        faceVertexPairs.Add(new FaceVertexPair(faceId, vId));
                    }
                }

                // Sort them if desired by faceId, then vertexId
                faceVertexPairs.Sort((a, b) => {
                    int cmp = a.FaceId.CompareTo(b.FaceId);
                    if (cmp == 0) cmp = a.VertexId.CompareTo(b.VertexId);
                    return cmp;
                });

                // 2) Make Gurobi environment and model
                using (GRBEnv env = new GRBEnv(true))
                {
                    env.Start();
                    using (GRBModel model = new GRBModel(env))
                    {
                        model.ModelName = "PureLocalVariables";

                        // 3) Create local variables fAll, plus lambda
                        CreateVariables(model, data);

                        // 4) Add equilibrium constraints:
                        AddEquilibriumConstraints(model, data);

                        // 5) Add friction & no-tension constraints
                        AddContactConstraints(model, data);

                        // 6)   Add face eccentricity constraints 
                        AddFaceEccConstraints(model, geometry, data);


                    // 7) Objective: maximize lambda
                    GRBLinExpr obj = 0.0;
                        obj.AddTerm(1.0, lambda);
                        model.SetObjective(obj, GRB.MAXIMIZE);
                        model.Write("debugModel.lp");
                    // 8) Solve
                        model.Optimize();

                        // 9) Print solution
                        PrintSolution(model);
                    }
                }
            }

            /// Create 2 local variables (f_n, f_t) for each face-vertex pair,
            /// plus the load factor lambda.
            /// The matrix A is expected to have #cols = 2 * faceVertexPairs.Count.
            private void CreateVariables(GRBModel model, ProblemData data)
            {
                int m = faceVertexPairs.Count; // number of pairs
                fAll = new GRBVar[2 * m];

                for (int j = 0; j < m; j++)
                {
                    // f_n >= 0
                    fAll[2 * j] = model.AddVar(
                        0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS,
                        $"fN_{j}"
                    );

                    // f_t unbounded
                    fAll[2 * j + 1] = model.AddVar(
                        -GRB.INFINITY, GRB.INFINITY, 0.0, GRB.CONTINUOUS,
                        $"fT_{j}"
                    );
                }

                // Also create lambda
                lambda = model.AddVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, "lambda");

                model.Update();
            }

            /// <summary>
            /// A * fAll = B * lambda
            /// Where A has #rows = data.NumRows, #cols = data.NumCols = 2 * (faceVertexPairs.Count).
            /// </summary>
            private void AddEquilibriumConstraints(GRBModel model, ProblemData data)
            {
                int n = data.NumRows;
                int p = data.NumCols;  // should match fAll.Length

                if (n == 0 || p == 0)
                    return;

                if (p != fAll.Length)
                {
                    Console.WriteLine($"Error: A-matrix has {p} columns but we have {fAll.Length} local variables. Mismatch!");
                    return;
                }

                for (int i = 0; i < n; i++)
                {
                    GRBLinExpr lhs = 0.0;
                    for (int j = 0; j < p; j++)
                    {
                        double valA = data.MatrixA[i, j];
                        if (Math.Abs(valA) > 1e-15)
                            lhs.AddTerm(valA, fAll[j]);
                    }
                    // Right side: B[i] * lambda
                    GRBLinExpr rhs = 0.0;
                    rhs.AddTerm(data.VectorB[i], lambda);

                    model.AddConstr(lhs == rhs, $"Equil_{i}");
                }
            }

            /// <summary>
            /// For each pair j, friction constraints: -mu * fN_j <= fT_j <= mu * fN_j
            /// plus no tension fN_j >= 0 (already in the variable bound).
            /// </summary>
            private void AddContactConstraints(GRBModel model, ProblemData data)
            {
                double mu = data.Mu;
                int m = faceVertexPairs.Count;

                for (int j = 0; j < m; j++)
                {
                    GRBVar fN = fAll[2 * j];
                    GRBVar fT = fAll[2 * j + 1];

                    // - mu * fN <= fT <= mu * fN
                    model.AddConstr(fT <= mu * fN, $"FricPos_{j}");
                    model.AddConstr(fT >= -mu * fN, $"FricNeg_{j}");
                }
            }

            private void AddFaceEccConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
            {
            double sigmaC = data.SigmaC;

               foreach (var fKvp in geometry.Faces)
               {
                Face face = fKvp.Value;
                // We assume each face in 2D has exactly 2 vertices.
                // If there's only 1 or 0, skip. If more than 2, handle carefully or skip.
                if (face.VertexIds.Count != 2) continue;

                int v1 = face.VertexIds[0];
                int v2 = face.VertexIds[1];

                // 1) Find the local indices in faceVertexPairs
                int j1 = FindFaceVertexIndex(face.Id, v1);
                int j2 = FindFaceVertexIndex(face.Id, v2);
                if (j1 < 0 || j2 < 0)
                    {
                    // If we can't find them, skip or throw an error
                    Console.WriteLine($"Error: missing face-vertex pair for face={face.Id}, v1={v1}, v2={v2}");
                    continue;
                    }

                // 2) Retrieve the local normal variables
                //    f_n(v1) = fAll[2*j1],  f_n(v2) = fAll[2*j2]
                GRBVar fn1 = fAll[2 * j1];
                GRBVar fn2 = fAll[2 * j2];

                // 3) Define a new variable f_n,k = fn1 + fn2
                GRBVar fnk = model.AddVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS,
                    $"fnk_face{face.Id}");
                  {
                    GRBLinExpr sumExpr = 0.0;
                    sumExpr.AddTerm(1.0, fn1);
                    sumExpr.AddTerm(1.0, fn2);
                    model.AddConstr(fnk == sumExpr, $"Def_fnk_{face.Id}");
                  }

                double Lk = face.Length;       // half-joint length in 2D
                double t_k = face.Thickness;

                // 4) Eccentricity variable e_k with range in [-Lk/2, Lk/2]
                GRBVar eK = model.AddVar(-Lk / 2.0, Lk / 2.0, 0.0, GRB.CONTINUOUS,
                    $"eK_face{face.Id}");

                // 5) Moment equilibrium (bilinear eqn):
                //    fnk * eK = (fn2 - fn1)*(Lk/2) => fnk*eK - (Lk/2)*fn2 + (Lk/2)*fn1 == 0
                   {
                    GRBQuadExpr mq = new GRBQuadExpr();
                    // + fnk * eK
                    mq.AddTerm(1.0, fnk, eK);
                    // - (Lk/2)*fn2
                    mq.AddTerm(-Lk / 2.0, fn2);
                    // + (Lk/2)*fn1
                    mq.AddTerm(Lk / 2.0, fn1);

                    model.AddQConstr(mq == 0.0, $"MomentEq_face_{face.Id}");
                    }

                    // 6) Strength limit constraints, if sigmaC > 0 and thickness > 0
                    //    eK + fnk/(2*sigmaC*t_k) <= Lk/2
                    //    eK - fnk/(2*sigmaC*t_k) >= -Lk/2
                    if (sigmaC > 1e-9 && t_k > 1e-9)
                    {
                    double denom = 2.0 * sigmaC * t_k;

                       {
                        // eK + (1/denom)*fnk <= Lk/2
                        GRBLinExpr lhsUp = 0.0;
                        lhsUp.AddTerm(1.0, eK);
                        lhsUp.AddTerm(1.0 / denom, fnk);
                        model.AddConstr(lhsUp <= Lk / 2.0, $"StrengthUp_{face.Id}");
                       }

                       {
                        // eK - (1/denom)*fnk >= -Lk/2
                        GRBLinExpr lhsLo = 0.0;
                        lhsLo.AddTerm(1.0, eK);
                        lhsLo.AddTerm(-1.0 / denom, fnk);
                        model.AddConstr(lhsLo >= -Lk / 2.0, $"StrengthLo_{face.Id}");
                       }
                    }
               }
            }

             /// Find the index j in faceVertexPairs for (faceId, vertexId).
             private int FindFaceVertexIndex(int faceId, int vertexId)
        {
            for (int j = 0; j < faceVertexPairs.Count; j++)
            {
                var fv = faceVertexPairs[j];
                if (fv.FaceId == faceId && fv.VertexId == vertexId)
                    return j;
            }
            return -1;
        }
        /// <summary>
        /// Example of printing the solution
        /// </summary>
        private void PrintSolution(GRBModel model)
            {
                int status = model.Status;
                if (model.SolCount > 0)
                {
                    Console.WriteLine($"Feasible solution found. Status={status}");
                    double lamVal = lambda.X;
                    Console.WriteLine($"lambda = {lamVal:F4}");

                    for (int j = 0; j < fAll.Length; j += 2)
                    {
                        double fn = fAll[j].X;
                        double ft = fAll[j + 1].X;
                        // faceVertexPairs[j/2] is the pair info
                        var pair = faceVertexPairs[j / 2];
                        Console.WriteLine($"Pair (face={pair.FaceId}, v={pair.VertexId}): fN={fn:F3}, fT={ft:F3}");
                    }
                }
                else
                {
                    Console.WriteLine($"No feasible solution. Status={status}");
                }
            if (model.Status == 4)
                //LOADED = 1,  OPTIMAL = 2, INFEASIBLE = 3, INF_OR_UNBD = 4, UNBOUNDED = 5

            {
                // Sometimes turning off presolve helps to distinguish infeasible from unbounded
                model.Parameters.Presolve = 0;
                model.Optimize();

                if (model.Status == 3)
                {
                    Console.WriteLine("Model is infeasible. Computing IIS...");
                    model.ComputeIIS();
                    model.Write("infeasible.ilp");
                }
                else if (model.Status == 5)
                {
                    Console.WriteLine("Model is unbounded!");
                    // Possibly do something else
                }
            }
        }
        
    } //end public class LocalOptimizer


    // -----------------------------------------------------------------
    // 3) Program.Main: file I/O for matrix and faces
    // -----------------------------------------------------------------

    internal class Program
    {
        static void Main(string[] args)
        {
            try
            {
                // 1) Create ProblemData and load eternal files 
                ProblemData pData = new ProblemData();
                string matrixFilePath = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\matrix_A_parallel.txt";
                LoadMatrixAndVector(pData, matrixFilePath);

                // 2) Create a geometry
                GeometryModel geometry = new GeometryModel();

                // faces
                string faceFilePath = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\face_parallel.txt";
                LoadFacesFromFile(geometry, faceFilePath);
                // 4) Possibly ensure we have all vertices that appear in faces
                AutoPopulateVertices(geometry);

                // 5) Solve with the local approach
                LocalOptimizer optimizer = new LocalOptimizer();
                optimizer.SolveProblem(geometry, pData);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error: " + ex.Message + "\n" + ex.StackTrace);
            }
            Console.WriteLine("Done geometry.");
            Console.ReadKey();
        }  // end of Main

        // load   Matrix A and Vector B from the specified format
        static void LoadMatrixAndVector(ProblemData pData, string filePath)
        {
            Console.WriteLine($"Loading matrix A, vector B from: {filePath}");
            if (!File.Exists(filePath)) throw new FileNotFoundException(filePath);

            bool readingMatrixSize = false;
            bool readingMatrixValues = false;
            bool readingVectorValues = false;
            int sizeLinesRead = 0;

            int n_local = 0, p_local = 0;
            List<double> matrixAValues = new List<double>();
            List<double> vectorBValues = new List<double>();

            var lines = File.ReadAllLines(filePath);
            for (int i = 0; i < lines.Length; i++)
            {
                string line = lines[i].Trim();
                if (string.IsNullOrEmpty(line) || line.StartsWith("#")) continue;
                if (line == "matrix_A_size")
                {
                    readingMatrixSize = true;
                    readingMatrixValues = false;
                    readingVectorValues = false;
                    sizeLinesRead = 0;
                    continue;
                }
                if (line == "matrix_A_values")
                {
                    readingMatrixSize = false;
                    readingMatrixValues = true;
                    readingVectorValues = false;
                    continue;
                }
                if (line == "vector_B_values")
                {
                    readingMatrixSize = false;
                    readingMatrixValues = false;
                    readingVectorValues = true;
                    continue;
                }
                // parse lines based on flags
                if (readingMatrixSize)
                {
                    if (sizeLinesRead == 0)
                        n_local = int.Parse(line);
                    else if (sizeLinesRead == 1)
                        p_local = int.Parse(line);
                    sizeLinesRead++;
                    if (sizeLinesRead == 2)
                    {
                        if (n_local <= 0 || p_local <= 0)
                            throw new InvalidDataException("Invalid matrix dimensions read.");

                        pData.MatrixA = new double[n_local, p_local];
                        Console.WriteLine($"Matrix dimension: {n_local} x {p_local}");
                        readingMatrixSize = false;
                    }
                    continue;
                }
                if (readingMatrixValues)
                {
                    // each line can have multiple values separated by spaces
                    var parts = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                    foreach (var part in parts)
                    {
                        if (!double.TryParse(part, NumberStyles.Any, CultureInfo.InvariantCulture, out double val))
                            throw new InvalidDataException($"Could not parse matrix value '{part}' on line {i + 1}.");
                        matrixAValues.Add(val);
                    }
                }
                else if (readingVectorValues)
                {
                    var parts = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                    foreach (var part in parts)
                    {
                        if (!double.TryParse(part, NumberStyles.Any, CultureInfo.InvariantCulture, out double val))
                            throw new InvalidDataException($"Could not parse vector B value '{part}' on line {i + 1}.");
                        vectorBValues.Add(val);
                    }
                }
            }

            // Fill pData.MatrixA
            if (pData.MatrixA == null)
                throw new InvalidDataException("Matrix A not allocated. Did we miss 'matrix_A_size' lines?");

            if (matrixAValues.Count != (pData.MatrixA.GetLength(0) * pData.MatrixA.GetLength(1)))
                throw new InvalidDataException("Mismatch in # of matrix A values vs. dimension.");

            int idx = 0;
            for (int r = 0; r < pData.MatrixA.GetLength(0); r++)
            {
                for (int c = 0; c < pData.MatrixA.GetLength(1); c++)
                {
                    pData.MatrixA[r, c] = matrixAValues[idx++];
                }
            }

            // Fill pData.VectorB
            if (vectorBValues.Count == 0)
                throw new InvalidDataException("No vector B values found in file.");

            if (vectorBValues.Count != pData.MatrixA.GetLength(0))
            {
                Console.WriteLine($"Warning: #B-values={vectorBValues.Count}, expected {pData.MatrixA.GetLength(0)}. Using them anyway...");
            }

            pData.VectorB = vectorBValues.ToArray();

            Console.WriteLine("Matrix A and Vector B loaded successfully.");

            // CHECK first few rows and columns of A
            Console.WriteLine("Some rows of Matrix A:");
            int numRowsToPrint = Math.Min(5, pData.NumRows);   // e.g., print up to 5 rows
            int numColsToPrint = Math.Min(15, pData.NumCols);   // e.g., print up to 5 columns
            for (int r = 0; r < numRowsToPrint; r++)
            {
                string rowStr = "";
                for (int c = 0; c < numColsToPrint; c++)
                {
                    rowStr += pData.MatrixA[r, c].ToString("F3") + " ";
                }
                Console.WriteLine($"Row {r}: {rowStr}");
            }
        }


        // Load Faces from file:   Format (per line): "face_ID_k", vertexId1, vertexId2, h_k, t_k
        static void LoadFacesFromFile(GeometryModel geometry, string faceFilePath)
        {
            Console.WriteLine($"Loading faces from: {faceFilePath}");
            if (!File.Exists(faceFilePath))
                throw new FileNotFoundException(faceFilePath);

            int lineNo = 0;
            var lines = File.ReadAllLines(faceFilePath);
            foreach (var rawLine in lines)
            {
                lineNo++;
                string line = rawLine.Trim();
                if (string.IsNullOrEmpty(line) || line.StartsWith("#"))
                    continue; // ignore blank or commented lines

                // split by comma
                string[] parts = line.Split(',');
                if (parts.Length < 5)
                {
                    Console.WriteLine($"Warning: line {lineNo}, expected 5 values (faceID,v1,v2,length,thickness). Found {parts.Length}. Skipping.");
                    continue;
                }

                if (!int.TryParse(parts[0], out int faceId))
                {
                    Console.WriteLine($"Cannot parse faceId at line {lineNo}. Skipping.");
                    continue;
                }
                if (!int.TryParse(parts[1], out int v1))
                {
                    Console.WriteLine($"Cannot parse vertex1 ID at line {lineNo}. Skipping.");
                    continue;
                }
                if (!int.TryParse(parts[2], out int v2))
                {
                    Console.WriteLine($"Cannot parse vertex2 ID at line {lineNo}. Skipping.");
                    continue;
                }

                if (!double.TryParse(parts[3], NumberStyles.Any, CultureInfo.InvariantCulture, out double length))
                {
                    Console.WriteLine($"Cannot parse face length at line {lineNo}. Skipping.");
                    continue;
                }
                if (!double.TryParse(parts[4], NumberStyles.Any, CultureInfo.InvariantCulture, out double thick))
                {
                    Console.WriteLine($"Cannot parse face thickness at line {lineNo}. Skipping.");
                    continue;
                }

                // Create the Face object
                Face newFace = new Face
                {
                    Id = faceId,
                    VertexIds = new List<int> { v1, v2 },
                    Length = length,
                    Thickness = thick
                };

                // Check for duplicate face ID
                if (geometry.Faces.ContainsKey(faceId))
                {
                    Console.WriteLine($"Warning: face ID={faceId} already exists, skipping duplicate line {lineNo}.");
                    continue;
                }

                geometry.Faces.Add(faceId, newFace);
            }

            Console.WriteLine($"Done loading faces. Found total {geometry.Faces.Count} faces.");

            // Quick listing:
            foreach (var kvp in geometry.Faces)
            {
                Face f = kvp.Value;
                Console.WriteLine($"Face ID={f.Id} with vertices ({f.VertexIds[0]}, {f.VertexIds[1]}) " +
                                  $"Length={f.Length}, Thickness={f.Thickness}");
            }
        }
        // Insert vertices if not already present
        static void AutoPopulateVertices(GeometryModel geometry)
        {
            foreach (var face in geometry.Faces.Values)
            {
                foreach (int vId in face.VertexIds)
                {
                    if (!geometry.Vertices.ContainsKey(vId))
                    {
                        geometry.Vertices[vId] = new ContactPoint { Id = vId };
                    }
                }
            }
            Console.WriteLine($"AutoPopulateVertices: now have {geometry.Vertices.Count} vertices in the model.");
        }

    } // end of Program
} //end namespace
