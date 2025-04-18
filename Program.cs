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
   Format (per line):   faceID, length, thickness, cohesionFlag, vertex1, vertex2, ....
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

        public double[] G { get; set; }
        public int NumRows => (MatrixA == null) ? 0 : MatrixA.GetLength(0);   // n
        public int NumCols => (MatrixA == null) ? 0 : MatrixA.GetLength(1);   // p



        // friction, compressive strength, etc.
        public double Mu { get; set; } = 0.4;    // friction
        public double SigmaC { get; set; } = 8200;  // compressive strength
        // For partial contact: thickness, etc., read from Face definitions.
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
        public double Depth { get; set; }
        public double Thickness { get; set; }

        // If you want a boolean for cohesion
        public double CohesionValue { get; set; }

        public List<int> VertexIds { get; set; } = new List<int>();
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
            //  storing Gurobi variables (eccVars) mapped by faceId
            private Dictionary<int, GRBVar> faceEccVars = new Dictionary<int, GRBVar>();
            // Similarly,  face-level normal sums
            private Dictionary<int, GRBVar> faceFnkVars = new Dictionary<int, GRBVar>();

        private int FindFaceVertexPairIndex(int faceId, int vertexId)
        {
            // faceVertexPairs is  List<FaceVertexPair>, presumably defined as a field in this class
            for (int j = 0; j < faceVertexPairs.Count; j++)
            {
                var pair = faceVertexPairs[j];
                if (pair.FaceId == faceId && pair.VertexId == vertexId)
                    return j;
            }
            return -1; // not found
        }

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
                        AddContactConstraints(model, geometry, data);

                        // 6)   Add face eccentricity constraints 
                        AddFaceEccConstraints(model, geometry, data);

                        // 7) Objective: maximize lambda
                        GRBLinExpr obj = 0.0;
                        obj.AddTerm(1.0, lambda);
                        model.SetObjective(obj, GRB.MAXIMIZE);
                        model.Write("debugModel.lp");
                       
                        // 8) Solve
                        model.Optimize();
                        SaveResultsToFile(model, @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\results.txt");
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


            /// A * fAll - G = B * lambda
            /// Where A has #rows = data.NumRows, #cols = data.NumCols = 2 * (faceVertexPairs.Count).
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
                // Subtract the gravity portion G[i]
                if (data.G != null && i < data.G.Length)
                    lhs.AddConstant(-data.G[i]);
                // Right side: B[i] * lambda
                GRBLinExpr rhs = 0.0;
                if (data.VectorB != null && i < data.VectorB.Length)
                    rhs.AddTerm(data.VectorB[i], lambda);

                model.AddConstr(lhs == rhs, $"Equil_{i}");
                }
            }


        /// For each pair j, friction constraints: -mu * fN_j <= fT_j <= mu * fN_j
        /// plus no tension fN_j >= 0 (already in the variable bound).
        private void AddContactConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            double mu = data.Mu;

            foreach (var kvp in geometry.Faces)
            {
                Face face = kvp.Value;
                if (face.VertexIds.Count < 2)
                    continue; // skip if face has <2 vertices (unexpected)


                // Use face-specific cohesion value 
                double cohesion = face.CohesionValue;
                double contactArea = face.Thickness * face.Depth;

                foreach (int vId in face.VertexIds)
                {
                    int pairIndex = FindFaceVertexPairIndex(face.Id, vId);
                    if (pairIndex < 0)
                    {
                        Console.WriteLine($"Warning: Missing face-vertex pair index for face={face.Id}, vertex={vId}");
                        continue;
                    }

                    GRBVar fN = fAll[2 * pairIndex];
                    GRBVar fT = fAll[2 * pairIndex + 1];

                    // Split cohesion contribution across the two vertices
                    double cohesionShare = 0.5 * cohesion * contactArea;

                    // Frictional constraints with cohesion
                    model.AddConstr(fT <= mu * fN + cohesionShare, $"FricPos_{face.Id}_{vId}");
                    model.AddConstr(fT >= -(mu * fN + cohesionShare), $"FricNeg_{face.Id}_{vId}");

                    // No tension constraint
                    model.AddConstr(fN >= 0.0, $"NoTension_{face.Id}_{vId}");
                }
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
                int idx1 = FindFaceVertexPairIndex(face.Id, v1);
                int idx2 = FindFaceVertexPairIndex(face.Id, v2);
                if (idx1 < 0 || idx2 < 0)
                {
                    // cannot proceed if the face–vertex pair wasn't found
                    Console.WriteLine($"Missing face–vertex pair for Face={face.Id}. Skipping ecc constraints.");
                    continue;
                }

                GRBVar fn1 = fAll[2 * idx1];   // local normal var for (face, v1)
                GRBVar fn2 = fAll[2 * idx2];   // local normal var for (face, v2)

                // 3) Define a new variable f_n,k = fn1 + fn2
                GRBVar fnk = model.AddVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, $"fnk_face{face.Id}");
                {
                    GRBLinExpr sumExpr = new GRBLinExpr();
                    sumExpr.AddTerm(1.0, fn1);
                    sumExpr.AddTerm(1.0, fn2);
                    model.AddConstr(fnk == sumExpr, $"Def_fnk_face{face.Id}");
                }

                double Lk = face.Depth;
                GRBVar eK = model.AddVar(-Lk / 2.0, Lk / 2.0, 0.0, GRB.CONTINUOUS, $"eK_face{face.Id}");

                faceEccVars[face.Id] = eK;

                // 3) Moment equilibrium:  fnk * eK = (fn2 - fn1)*(Lk/2)
                GRBQuadExpr mq = new GRBQuadExpr();
                mq.AddTerm(1.0, fnk, eK);        // fnk*eK
                mq.AddTerm(-Lk / 2.0, fn2);      // -(Lk/2)*fn2
                mq.AddTerm(Lk / 2.0, fn1);       // +(Lk/2)*fn1
                model.AddQConstr(mq == 0.0, $"MomentEq_face{face.Id}");

                //    fnk * eK = (fn2 - fn1)*(Lk/2) => fnk*eK - (Lk/2)*fn2 + (Lk/2)*fn1 == 0
                //  Strength limit constraints, if sigmaC > 0 and thickness > 0
                // 4) Strength limit constraints
                double t_k = face.Thickness;
                if (sigmaC > 1e-9 && t_k > 1e-9)
                {
                    double denom = 2.0 * sigmaC * t_k;

                    // eK + fnk/denom <= Lk/2
                    {
                        GRBLinExpr lhsUp = 0.0;
                        lhsUp.AddTerm(1.0, eK);
                        lhsUp.AddTerm(1.0 / denom, fnk);
                        model.AddConstr(lhsUp <= Lk / 2.0, $"StrengthUp_face{face.Id}");
                    }

                    // eK - fnk/denom >= -Lk/2
                    {
                        GRBLinExpr lhsLo = 0.0;
                        lhsLo.AddTerm(1.0, eK);
                        lhsLo.AddTerm(-1.0 / denom, fnk);
                        model.AddConstr(lhsLo >= -Lk / 2.0, $"StrengthLo_face{face.Id}");
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

        ///  printing the solution
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

        private void SaveResultsToFile(GRBModel model, string resultsFilePath)
        {
            // If no solution, do nothing
            if (model.SolCount == 0)
            {
                Console.WriteLine("No feasible solution, nothing to save.");
                return;
            }

            // Ensure the output directory exists
            string directoryPath = Path.GetDirectoryName(resultsFilePath);
            if (!Directory.Exists(directoryPath))
            {
                Directory.CreateDirectory(directoryPath);
            }

            using (StreamWriter writer = new StreamWriter(resultsFilePath))
            {
                // 1) Save lambda
                try
                {
                    double lamVal = lambda.Get(GRB.DoubleAttr.X);
                    writer.WriteLine($"lambda = {lamVal:F4}");
                    Console.WriteLine($"lambda = {lamVal:F4}");
                }
                catch (GRBException e)
                {
                    writer.WriteLine($"Could not retrieve lambda: {e.Message}");
                }

                // 2) Save local (f_n, f_t) for each face–vertex pair
                for (int j = 0; j < faceVertexPairs.Count; j++)
                {
                    int indexFn = 2 * j;
                    int indexFt = 2 * j + 1;

                    double fnVal, ftVal;
                    try
                    {
                        fnVal = fAll[indexFn].Get(GRB.DoubleAttr.X);
                        ftVal = fAll[indexFt].Get(GRB.DoubleAttr.X);
                    }
                    catch (GRBException e)
                    {
                        writer.WriteLine($"Could not retrieve forces for pair j={j}: {e.Message}");
                        continue;
                    }

                    var pair = faceVertexPairs[j]; // (faceId, vertexId)
                    writer.WriteLine($"Face {pair.FaceId}, Vertex {pair.VertexId}: " +
                                     $"fn={fnVal:F3}, ft={ftVal:F3}");
                }
                // 4) If we have face eccentricities, print them too
                foreach (var kvp in faceEccVars)
                {
                    int faceId = kvp.Key;
                    GRBVar eccVar = kvp.Value;
                    double eccVal = eccVar.X;  // or eccVar.Get(GRB.DoubleAttr.X)
                    writer.WriteLine($"Face {faceId}: eccentricity = {eccVal:F3}");
                }

                // 3) Compute and save total forces per face
                writer.WriteLine();
                writer.WriteLine("=== Face-level Total Forces ===");

                var faceTotalForces = new Dictionary<int, (double fnSum, double ftSum)>();

                for (int j = 0; j < faceVertexPairs.Count; j++)
                {
                    int indexFn = 2 * j;
                    int indexFt = 2 * j + 1;

                    double fnVal = fAll[indexFn].Get(GRB.DoubleAttr.X);
                    double ftVal = fAll[indexFt].Get(GRB.DoubleAttr.X);

                    var pair = faceVertexPairs[j]; // (faceId, vertexId)

                    if (!faceTotalForces.ContainsKey(pair.FaceId))
                        faceTotalForces[pair.FaceId] = (0.0, 0.0);

                    var current = faceTotalForces[pair.FaceId];
                    faceTotalForces[pair.FaceId] = (current.fnSum + fnVal, current.ftSum + ftVal);
                }

                // Write totals
                foreach (var kvp in faceTotalForces)
                {
                    writer.WriteLine($"Face {kvp.Key}: total_fn = {kvp.Value.fnSum:F3}, total_ft = {kvp.Value.ftSum:F3}");
                }

            }

            Console.WriteLine($"Results saved to {resultsFilePath}");
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
                ProblemData data = new ProblemData();
                string matrixFilePath = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\matrix_A_parallel.txt";
                LoadMatrixAndVector(data, matrixFilePath);

                // 2) Create a geometry
                GeometryModel geometry = new GeometryModel();

                // faces
                string faceFilePath = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\face_parallel.txt";
                LoadFacesFromFile(geometry, faceFilePath);
                if (!string.IsNullOrEmpty(faceFilePath))
                {
                    LoadFacesFromFile(geometry, faceFilePath);
                    // Then AddContactConstraints, EccConstraints, etc.
                }
                else
                {
                    Console.WriteLine("No face file. Skipping face constraints entirely.");
                }


                // 4) Possibly ensure we have all vertices that appear in faces
                AutoPopulateVertices(geometry);

                // 5) Solve with the local approach
                LocalOptimizer optimizer = new LocalOptimizer();
                optimizer.SolveProblem(geometry, data);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error: " + ex.Message + "\n" + ex.StackTrace);
            }
            Console.WriteLine("Done geometry.");
            Console.ReadKey();
        }  // end of Main

        // load   Matrix A and Vector B from the specified format
        static void LoadMatrixAndVector(ProblemData data, string filePath)
        {
            Console.WriteLine($"Loading matrix A, vector B from: {filePath}");
            if (!File.Exists(filePath)) throw new FileNotFoundException(filePath);

            bool readingMatrixSize = false;
            bool readingMatrixValues = false;
            bool readingVectorValues = false;
            bool readingGValues = false;
            int sizeLinesRead = 0;

            int n_local = 0, p_local = 0;
            List<double> matrixAValues = new List<double>();
            List<double> vectorBValues = new List<double>();
            List<double> gravityValues = new List<double>();

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
                if (line == "vector_G_values")
                {
                    readingMatrixSize = false;
                    readingMatrixValues = false;
                    readingVectorValues = false;
                    readingGValues = true;
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

                        data.MatrixA = new double[n_local, p_local];
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
                else if (readingGValues)
                {
                    var parts = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
                    foreach (var part in parts)
                    {
                        if (!double.TryParse(part, NumberStyles.Any, CultureInfo.InvariantCulture, out double val))
                            throw new InvalidDataException($"Could not parse gravity G value '{part}' on line {i + 1}.");
                        gravityValues.Add(val);
                    }
                }
            }
                // Fill data.MatrixA
                if (data.MatrixA == null)
                    throw new InvalidDataException("Matrix A not allocated. Did we miss 'matrix_A_size' lines?");

                if (matrixAValues.Count != (data.MatrixA.GetLength(0) * data.MatrixA.GetLength(1)))
                    throw new InvalidDataException("Mismatch in # of matrix A values vs. dimension.");

                int idx = 0;
                for (int r = 0; r < data.MatrixA.GetLength(0); r++)
                {
                    for (int c = 0; c < data.MatrixA.GetLength(1); c++)
                    {
                        data.MatrixA[r, c] = matrixAValues[idx++];
                    }
                }

                // Fill data.VectorB
                if (vectorBValues.Count == 0)
                    throw new InvalidDataException("No vector B values found in file.");

                if (vectorBValues.Count != data.MatrixA.GetLength(0))
                {
                    Console.WriteLine($"Warning: #B-values={vectorBValues.Count}, expected {data.MatrixA.GetLength(0)}. Using them anyway...");
                }
                data.VectorB = vectorBValues.ToArray();

                // Fill data.G
                if (gravityValues.Count != data.MatrixA.GetLength(0))
                {
                    Console.WriteLine($"Warning: #G-values={gravityValues.Count}, expected {data.MatrixA.GetLength(0)}. Using them anyway...");
                }
                data.G = gravityValues.ToArray();

                Console.WriteLine("Matrix A, Vector B, and Gravity G loaded successfully.");

                // CHECK first few rows and columns of A
                Console.WriteLine("Some rows of Matrix A:");
                int numRowsToPrint = Math.Min(5, data.NumRows);   // e.g., print up to 5 rows
                int numColsToPrint = Math.Min(15, data.NumCols);   // e.g., print up to 5 columns
                for (int r = 0; r < numRowsToPrint; r++)
                {
                    string rowStr = "";
                    for (int c = 0; c < numColsToPrint; c++)
                    {
                        rowStr += data.MatrixA[r, c].ToString("F3") + " ";
                    }
                    Console.WriteLine($"Row {r}: {rowStr}");
                }
           
        }

        // Load Faces from file:   Format (per line): faceID, length, thickness, cohesionFlag, vertex1, vertex2, ...
        static void LoadFacesFromFile(GeometryModel geometry, string faceFilePath)
        {
            // If the user didn't provide a faceFilePath, or it's empty, skip or return
            if (string.IsNullOrEmpty(faceFilePath))
            {
                Console.WriteLine("No face file provided. Skipping face constraints altogether.");
                return;
            }
            Console.WriteLine($"Loading faces from: {faceFilePath}");
            if (!File.Exists(faceFilePath))
                throw new FileNotFoundException(faceFilePath);

            var lines = File.ReadAllLines(faceFilePath);
            int lineNo = 0;
            foreach (var rawLine in lines)
            {
                lineNo++;
                string line = rawLine.Trim();
                // Skip empty or commented lines
                if (string.IsNullOrEmpty(line) || line.StartsWith("#"))
                    continue;

                // Split by comma
                string[] parts = line.Split(',');
                // We expect at least 5 columns for: faceID, length, thickness, cohesionValue, plus >=1 vertex
                if (parts.Length < 5)
                {
                    Console.WriteLine($"Warning: line {lineNo}, expected at least 5 columns. Found {parts.Length}. Skipping this line.");
                    continue;
                }

                // 1) Parse the first 4 fields
                if (!int.TryParse(parts[0], out int faceId))
                {
                    Console.WriteLine($"Could not parse faceId at line {lineNo}. Skipping.");
                    continue;
                }
                if (!double.TryParse(parts[1], NumberStyles.Any, CultureInfo.InvariantCulture, out double faceLength))
                {
                    Console.WriteLine($"Could not parse face length at line {lineNo}. Skipping.");
                    continue;
                }
                if (!double.TryParse(parts[2], NumberStyles.Any, CultureInfo.InvariantCulture, out double faceThickness))
                {
                    Console.WriteLine($"Could not parse face thickness at line {lineNo}. Skipping.");
                    continue;
                }

                if (!double.TryParse(parts[3], NumberStyles.Any, CultureInfo.InvariantCulture, out double cohesionValue))
                {
                    Console.WriteLine($"Could not parse cohesion value at line {lineNo}. Skipping.");
                    continue;
                }

                // 2) Extract the vertex IDs from the remaining fields
                var vertexList = new List<int>();
                for (int i = 4; i < parts.Length; i++)
                {
                    if (int.TryParse(parts[i], out int vId))
                    {
                        vertexList.Add(vId);
                    }
                    else
                    {
                        Console.WriteLine($"Could not parse vertex ID '{parts[i]}' at line {lineNo}. Skipping that ID.");
                    }
                }
                if (vertexList.Count < 2)
                {
                    Console.WriteLine($"Warning: face {faceId} has fewer than 2 vertices listed. This might cause issues. We proceed anyway.");
                }

                // Create and register the Face object
                if (geometry.Faces.ContainsKey(faceId))
                {
                    Console.WriteLine($"Warning: face ID={faceId} already exists, skipping duplicate line {lineNo}.");
                    continue;
                }

                // 3) Create the Face object
                Face newFace = new Face
                {
                    Id = faceId,
                    Depth = faceLength,
                    Thickness = faceThickness,
                    CohesionValue = cohesionValue,
                    VertexIds = vertexList,
                };

                geometry.Faces.Add(faceId, newFace);
            }

            Console.WriteLine($"Done loading faces. Found total {geometry.Faces.Count} faces.");

            // Quick listing:
            foreach (var kvp in geometry.Faces)
            {
                var f = kvp.Value;
                string vStr = string.Join(",", f.VertexIds);
                Console.WriteLine($"Face ID={f.Id}, L={f.Depth:F3}, t={f.Thickness:F3},Cohesion={f.CohesionValue:F3}, Vertices=[{vStr}]");
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
