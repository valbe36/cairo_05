using System.Globalization;
using Gurobi;
using System.IO; 
using System.Collections.Generic;
using static InterlockingMasonryLocalForces.ProblemData;

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
   Format (per line):   faceID, length, thickness, cohesion, friction, vertex1, vertex2.normal x, normal y, tangent x, tangent y
The VertexIds should be ordered consistently, such that when moving along the tangent 
from the first vertex to the second vertex, the normal points from block j to block j+1.
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
        public double[,] MatrixA { get; set; }
        public double[] B { get; set; }
        public double[] G { get; set; }
        public int NumBlocks { get; set; }
        public double Mu { get; set; } = 0.4;
        public double SigmaC { get; set; } = 8200;

        public int NumRows => MatrixA?.GetLength(0) ?? 0;
        public int NumCols => MatrixA?.GetLength(1) ?? 0;


        // Add this helper to track column indices
        private static Dictionary<(int FaceId, int VertexId), int> _columnMap;

        public static double[,] BuildEquilibriumMatrix(GeometryModel geometry, int numBlocks)
        {
            // Filter out support blocks (ID ≤ 0)
            var nonSupportBlocks = geometry.Blocks.Values
                .Where(b => b.Id > 0)
                .OrderBy(b => b.Id)
                .ToList();

            int numNonSupportBlocks = nonSupportBlocks.Count;
            Dictionary<int, int> blockRowMap = nonSupportBlocks
                .Select((b, idx) => new { b.Id, idx })
                .ToDictionary(x => x.Id, x => x.idx * 3); // Map to row group (3 rows per block)

            // Step 1: Map each (face, vertex) pair to a unique column index
            Dictionary<(int FaceId, int VertexId), int> columnMap = new Dictionary<(int, int), int>();
            int colIdx = 0;
            foreach (var face in geometry.Faces.Values.OrderBy(f => f.Id))
            {
                foreach (int vId in face.VertexIds.OrderBy(v => v))
                {
                    columnMap[(face.Id, vId)] = colIdx++;
                }
            }

            // Step 2: Initialize matrix
            int numRows = numBlocks * 3;
            int numCols = columnMap.Count * 2; // 2 variables (N, T) per vertex
            double[,] matrixA = new double[numRows, numCols];

            // Step 3: Populate matrix
            foreach (var face in geometry.Faces.Values)
            {
                // Skip faces between two support blocks
                bool isJSupport = face.BlockJ <= 0;
                bool isKSupport = face.BlockK <= 0;
                if (isJSupport && isKSupport) continue;

                double[] n = face.Normal;
                double[] t = face.Tangent;

                foreach (int vId in face.VertexIds)
                {
                    ContactPoint vertex = geometry.Vertices[vId];
                    int colN = columnMap[(face.Id, vId)] * 2;
                    int colT = colN + 1;

                    // Handle Block J if it's a non-support block
                    if (!isJSupport && blockRowMap.TryGetValue(face.BlockJ, out int rowJ))
                    {
                        Block blockJ = geometry.Blocks[face.BlockJ];
                        double xRelJ = vertex.X - blockJ.CentroidX;
                        double yRelJ = vertex.Y - blockJ.CentroidY;

                        matrixA[rowJ, colN] = -n[0];         // ΣFx (N) Fx = -N*nx - T*tx
                        matrixA[rowJ, colT] = -t[0];         // ΣFx (T)
                        matrixA[rowJ + 1, colN] = -n[1];     // ΣFy (N)
                        matrixA[rowJ + 1, colT] = -t[1];     // ΣFy (T)
                        matrixA[rowJ + 2, colN] = -xRelJ * n[1] - yRelJ * n[0];  // ΣM (N)
                        matrixA[rowJ + 2, colT] = -xRelJ * t[1] - yRelJ * t[0];   // ΣM (T)
                    }

                    // Handle Block K if it's a non-support block
                    if (!isKSupport && blockRowMap.TryGetValue(face.BlockK, out int rowK))
                    {
                        Block blockK = geometry.Blocks[face.BlockK];
                        double xRelK = vertex.X - blockK.CentroidX;
                        double yRelK = vertex.Y - blockK.CentroidY;

                        matrixA[rowK, colN] = n[0];        // ΣFx (N) reversed
                        matrixA[rowK, colT] = t[0];        // ΣFx (T) reversed
                        matrixA[rowK + 1, colN] = n[1];    // ΣFy (N) reversed
                        matrixA[rowK + 1, colT] = t[1];    // ΣFy (T) reversed
                        matrixA[rowK + 2, colN] = (xRelK * n[1] - yRelK * n[0]);  // ΣM (N) reversed
                        matrixA[rowK + 2, colT] = (xRelK * t[1] - yRelK * t[0]);   // ΣM (T) reversed
                    }
                }
            }
            return matrixA;
        }

        private static int GetColumnIndex(int faceId, int vertexId)
        {
            return _columnMap.TryGetValue((faceId, vertexId), out int idx) ? idx : -1;
        }

    }

    /// Represents a single contact point (vertex).
    public class ContactPoint
    {
        public int Id { get; set; }
        public double X { get; set; }  // Global X coordinate
        public double Y { get; set; }  // Global Y coordinate
    }

    /// Represents a 2D contact face, assumed 2 vertexes, but adaptable  <summary>
    ///   faceID, length, thickness, cohesion, friction, vertex1, 
    ///   vertex2,normalX, normalY, blok j, block j+1
    /// </summary>
    public class Face
    {
        public int Id { get; set; }
        public double Depth { get; set; }
        public double Thickness { get; set; }

        public double CohesionValue { get; set; }
        public double? MuOverride { get; set; }

        public List<int> VertexIds { get; set; } = new List<int>();
        public double[] Normal { get; set; } = new double[2];   // [nx, ny], unit vector
        public double[] Tangent { get; set; } = new double[2];  // [tx, ty], unit vector
        public int BlockJ { get; set; }  // Block on the "from" side of the normal
        public int BlockK { get; set; }  // Block on the "to" side of the normal
    }

    public class Block
    {
        public int Id { get; set; }
        public double CentroidX { get; set; }  // Centroid X
        public double CentroidY { get; set; }  // Centroid Y
    }

    /// A container for all geometry in the problem (faces, vertices, etc.).
    public class GeometryModel
    {
        public Dictionary<int, ContactPoint> Vertices { get; } = new Dictionary<int, ContactPoint>();
        public Dictionary<int, Face> Faces { get; } = new Dictionary<int, Face>();
        public Dictionary<int, Block> Blocks { get; } = new Dictionary<int, Block>();

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
            private GeometryModel _geometry;
            // Gurobi variable array (size = 2 * faceVertexPairs.Count):
            // the even indices are f_n, the odd are f_t, for that pair.
            private GRBVar[] fAll;
            private GRBVar lambda;
            //  storing Gurobi variables (eccVars) mapped by faceId
            private Dictionary<int, GRBVar> faceEccVars = new Dictionary<int, GRBVar>();
            // maps (faceId, vertexId) -> pair index
            private Dictionary<(int face, int vtx), int> pairIndexMap;
          
        private int GetPairIndex(int faceId, int vertexId) =>
            pairIndexMap.TryGetValue((faceId, vertexId), out int idx) ? idx : -1;

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

                pairIndexMap = new Dictionary<(int, int), int>(faceVertexPairs.Count);
                for (int j = 0; j < faceVertexPairs.Count; j++)
                {
                var p = faceVertexPairs[j];
                pairIndexMap[(p.FaceId, p.VertexId)] = j;
                }

            _geometry = geometry;
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
                        DumpColumnMap(data);

                    // 8) Solve
                    model.Optimize();
                        SaveResultsToFile(model, @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\results_parallel.txt");
                    // 9) Print solution
                    PrintSolution(model);

                    }
                }
            }

        private void ValidateForceDirections()
        {
            foreach (var face in _geometry.Faces.Values)
            {
                foreach (int vId in face.VertexIds)
                {
                    int idx = GetPairIndex(face.Id, vId);
                    double fn = fAll[2 * idx].X;
                    double ft = fAll[2 * idx + 1].X;

                    // Convert to global forces
                    double Fx = fn * face.Normal[0] + ft * face.Tangent[0];
                    double Fy = fn * face.Normal[1] + ft * face.Tangent[1];

                    Console.WriteLine($"Face {face.Id} Vtx {vId}: F_global=({Fx:F2}, {Fy:F2})");
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
                if (data.B != null && i < data.B.Length)
                    rhs.AddTerm(data.B[i], lambda);

                model.AddConstr(lhs == rhs, $"Equil_{i}");
                }
            }


        /// For each pair j, friction constraints: -mu * fN_j <= fT_j <= mu * fN_j
        /// plus no tension fN_j >= 0 (already in the variable bound).
        private void AddContactConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            foreach (var kvp in geometry.Faces)        // <‑‑ loop defines kvp
            {
                Face face = kvp.Value;
                double mu = face.MuOverride ?? data.Mu;   // per‑face or default

                // If mu‑value ≤ 0  ⇒ treat as perfectly smooth face
                bool frictionless = mu <= 0.0;

                double area = face.Thickness * face.Depth;
                double cohShare = 0.5 * face.CohesionValue * area;   // kN

                foreach (int vId in face.VertexIds)
                {
                    int idx = GetPairIndex(face.Id, vId);
                    if (idx < 0) continue;       // safety

                    GRBVar fN = fAll[2 * idx];        // normal
                    GRBVar fT = fAll[2 * idx + 1];    // tangential

                    //  No‑tension bound (compression only)
                    model.AddConstr(fN >= 0.0, $"NoTension_{face.Id}_{vId}");

                    if (frictionless)
                    {
                        // shear must vanish
                        model.AddConstr(fT == 0.0, $"NoFric_{face.Id}_{vId}");
                    }
                    else
                    {
                        // Mohr‑Coulomb with cohesion
                        model.AddConstr(fT <= mu * fN + cohShare, $"Fric+_{face.Id}_{vId}");
                        model.AddConstr(fT >= -mu * fN - cohShare, $"Fric-_{face.Id}_{vId}");
                    }
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

                // Get vertex IDs
                int vertexId1 = face.VertexIds[0];
                int vertexId2 = face.VertexIds[1];

                // Get vertex coordinates
                ContactPoint vertex1 = _geometry.Vertices[vertexId1]; // Renamed from v1
                ContactPoint vertex2 = _geometry.Vertices[vertexId2]; // Renamed from v2

                // Get pair indices
                int idx1 = GetPairIndex(face.Id, vertexId1);
                int idx2 = GetPairIndex(face.Id, vertexId2);

                if (idx1 < 0 || idx2 < 0)
                {
                  Console.WriteLine($"Skipping ecc constraints for Face {face.Id}.");
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

                // 3) Moment equilibrium:   
                // Compute moment arm direction based on vertex order
                double sign = 1.0;
                if (vertex1.X > vertex2.X) sign = -1.0; 


                //Moment equilibrium:  fnk * eK = (fn2 - fn1)*(Lk/2)
                GRBQuadExpr mq = new GRBQuadExpr();
                mq.AddTerm(1.0, fnk, eK);
                mq.AddTerm(-sign * Lk / 2.0, fn2); // Adjusted sign
                mq.AddTerm(sign * Lk / 2.0, fn1);
                model.AddQConstr(mq == 0.0, $"MomentEq_face{face.Id}");

                //    fnk * eK = (fn2 - fn1)*(Lk/2) => fnk*eK - (Lk/2)*fn2 + (Lk/2)*fn1 == 0
                //   if sigmaC > 0 and thickness > 0, Strength limit constraints
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

        private void DumpColumnMap(ProblemData data)
        {
            string path = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\mapping.txt";          // choose any writable folder
            using var w = new StreamWriter(path);

            w.WriteLine("col | face vertex | A[0,col]");
            w.WriteLine("-------------------------------");

            for (int j = 0; j < faceVertexPairs.Count; j++)
            {
                var p = faceVertexPairs[j];
                double a0 = data.MatrixA[0, j];            // sample entry from row 0

                w.WriteLine($"{j,3} | {p.FaceId,4} {p.VertexId,6} | {a0,10:F3}");
            }

            Console.WriteLine($"Column map dumped to {path}");
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
                
                
                // 3) If we have face eccentricities, print them too
                foreach (var kvp in faceEccVars)
                {
                    int faceId = kvp.Key;
                    GRBVar eccVar = kvp.Value;
                    double eccVal = eccVar.X;  // or eccVar.Get(GRB.DoubleAttr.X)
                    writer.WriteLine($"Face {faceId}: eccentricity = {eccVal:F3}");

                    // double halfDepth = _geometry.Faces[faceId].Depth / 2.0;
                   // writer.WriteLine($"Face {faceId}: eccentricity ratio = {eccVal / halfDepth:F3}");
                }

                // 4) Compute and save total forces per face
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
                ProblemData data = new ProblemData();
                GeometryModel geometry = new GeometryModel();

                    // Load faces and geometry 
                LoadAllData(@"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\data_parallel_friction_0e4.txt"
                , geometry, data);


                // Build equilibrium matrix programmatically
                data.NumBlocks = geometry.Blocks.Values.Count(b => b.Id > 0); // Exclude supports
                data.MatrixA = ProblemData.BuildEquilibriumMatrix(geometry, data.NumBlocks);

                // add matrix debug output
                if (data.NumRows != data.B.Length)
                        throw new Exception("Matrix rows must match VectorB length.");
                    Console.WriteLine("First 5 rows of Matrix A:");
                    for (int r = 0; r < Math.Min(5, data.NumRows); r++)
                    {
                        string row = "";
                        for (int c = 0; c < Math.Min(10, data.NumCols); c++)
                            row += $"{data.MatrixA[r, c]:F2} ";
                        Console.WriteLine($"Row {r}: {row}");
                    }
          
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
 

                    DumpFaceData(geometry, "interlocking");   // or "plain", etc.
                 
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


        static void LoadAllData(string filePath, GeometryModel geometry, ProblemData data)
        {
            string currentSection = "";
            var vectorB = new List<double>();
            var vectorG = new List<double>();
            int lineNo = 0;

            foreach (var line in File.ReadAllLines(filePath))
            {
                lineNo++;
                string trimmedLine = line.Trim();
                if (string.IsNullOrEmpty(trimmedLine)) continue;

                if (trimmedLine.StartsWith("[") && trimmedLine.EndsWith("]"))
                {
                    currentSection = trimmedLine;
                    continue;
                }

                switch (currentSection)
                {
                    case "[Blocks]":
                        if (trimmedLine.StartsWith("BlockID")) continue; // Skip header
                        var blockParts = trimmedLine.Split(',');
                        //int blockId = int.Parse(blockParts[0]);
                        if (!int.TryParse(blockParts[0], out int blockId))
                        {
                            Console.WriteLine($"Line {lineNo}: Invalid BlockID '{blockParts[0]}'. Skipping block.");
                            continue;
                        }
                        geometry.Blocks.Add(blockId, new Block
                        {
                            Id = blockId,
                            CentroidX = double.Parse(blockParts[1]),
                            CentroidY = double.Parse(blockParts[2])
                        });
                        break;

                    case "[Vertices]":
                        if (trimmedLine.StartsWith("VertexID")) continue;
                        var vertexParts = trimmedLine.Split(',');
                        //int vId = int.Parse(vertexParts[0]);
                        if (!int.TryParse(vertexParts[0], out int vId))
                        {
                            Console.WriteLine($"Line {lineNo}: Invalid BlockID '{vertexParts[0]}'. Skipping block.");
                            continue;
                        }
                        geometry.Vertices.Add(vId, new ContactPoint
                        {
                            Id = vId,
                            X = double.Parse(vertexParts[1]),
                            Y = double.Parse(vertexParts[2])
                        });
                        break;

                    case "[Faces]":
                        if (trimmedLine.StartsWith("FaceID")) continue;
                        var faceParts = trimmedLine.Split(',');

                        // Parse basic face properties
                        Face face = new Face
                        {
                            Id = int.Parse(faceParts[0]),
                            BlockJ = int.Parse(faceParts[1]),
                            BlockK = int.Parse(faceParts[2]),
                            VertexIds = new List<int> { int.Parse(faceParts[3]), int.Parse(faceParts[4]) },
                            Depth = double.Parse(faceParts[9]),
                            Thickness = double.Parse(faceParts[10]),
                            CohesionValue = double.Parse(faceParts[11]),
                            MuOverride = double.Parse(faceParts[12])
                        };

                        // Parse vectors
                        double nx = double.Parse(faceParts[5]);
                        double ny = double.Parse(faceParts[6]);
                        double tx = double.Parse(faceParts[7]);
                        double ty = double.Parse(faceParts[8]);

                        // --- VALIDATION CHECKS ADDED HERE ---
                        // Check normal-tangent orthogonality
                        double dotProduct = nx * tx + ny * ty;
                        if (Math.Abs(dotProduct) > 1e-6)
                            Console.WriteLine($"Line {lineNo}: Face {face.Id} normal/tangent not orthogonal.");

                        // Check normal vector length
                        double normalLength = Math.Sqrt(nx * nx + ny * ny);
                        if (normalLength <= 1e-8)
                        {
                            Console.WriteLine($"Line {lineNo}: Invalid normal vector. Skipping face.");
                            continue;
                        }
                        else if (Math.Abs(normalLength - 1.0) > 1e-6)
                        {
                            Console.WriteLine($"Line {lineNo}: Normalizing face {face.Id} normal.");
                            nx /= normalLength;
                            ny /= normalLength;
                        }

                        face.Normal = new double[] { nx, ny };
                        face.Tangent = new double[] { tx, ty };
                        geometry.Faces.Add(face.Id, face);
                        break;

                    case "[VectorB]":
                        vectorB.AddRange(ParseLoadLine(trimmedLine, "VectorB", lineNo));
                        break;

                    case "[VectorG]":
                        vectorG.AddRange(ParseLoadLine(trimmedLine, "VectorG", lineNo));
                        break;
                }
            }
                // Validate VectorB and VectorG lengths
                int expectedRows = geometry.Blocks.Count * 3;
                ValidateVectorLength(vectorB, expectedRows, "VectorB");
                ValidateVectorLength(vectorG, expectedRows, "VectorG");

                data.B = vectorB.ToArray();
                data.G = vectorG.ToArray();

            // Output the final counts for verification
            Console.WriteLine($"Loaded {geometry.Blocks.Count} blocks.");
            Console.WriteLine($"Loaded {geometry.Vertices.Count} vertices.");
            Console.WriteLine($"Loaded {geometry.Faces.Count} faces.");
            Console.WriteLine($"Loaded {data.B.Length} entries for VectorB."); // Optional: verify vector lengths
            Console.WriteLine($"Loaded {data.G.Length} entries for VectorG.");

            // Cross-check: Number of faces should be equal to number of columns in the matrix divided by 4
            // Each face has 2 vertices, and for each face-vertex pair, there are 2 variables (fN and fT)
            // Total variables = number of face-vertex pairs * 2 = geometry.Faces.Count * 2 * 2
            int expectedNumVariables = geometry.Faces.Count * 2 * 2;

                // Assuming the matrix A will be built later and will have columns equal to the number of variables
                //  placeholder for the actual number of variables
                int actualNumVariables = data.NumCols;
                if (actualNumVariables != expectedNumVariables)
                {
                    Console.WriteLine($"Warning: The number of variables ({actualNumVariables}) does not match the expected number ({expectedNumVariables}). Please check data consistency.");
                }

                // Similarly,  cross-check the number of faces with the number of vertices / 2
                int expectedNumFacesFromVertices = geometry.Vertices.Count / 2;
                if (geometry.Faces.Count != expectedNumFacesFromVertices)
                {
                    Console.WriteLine($"Warning: Number of faces ({geometry.Faces.Count}) does not equal number of vertices / 2 ({expectedNumFacesFromVertices}). Please verify your data.");
                }
          

        }
        





// Helper method to parse a load line
private static IEnumerable<double> ParseLoadLine(string line, string vectorName, int lineNo)
        {
            try
            {
                return line.Split(',')
                           .Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture));
            }
            catch
            {
                throw new FormatException($"Invalid value in {vectorName} at line {lineNo}.");
            }
        }

        // Helper method to validate vector length
        private static void ValidateVectorLength(List<double> vector, int expected, string name)
        {
            if (vector.Count != expected)
                throw new InvalidDataException(
                    $"{name} has {vector.Count} entries. Expected {expected} (3 per block)."
                );
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

        // --- DumpFaceData Method (No changes needed based on analysis) ---
        static void DumpFaceData(GeometryModel g, string tag)
        {
            Console.WriteLine($"\n--- {tag} ---");
            foreach (var f in g.Faces.Values.OrderBy(f => f.Id))
            {
                // Use Null Coalescing Operator: If f.MuOverride is null, use 0.4
                double mu = f.MuOverride; // Assuming MuOverride is double, not double?
                                          // If it IS double?, then use: double mu = f.MuOverride ?? 0.4;
                                          // Based on parsing logic, it should always have a value now.
                // Let's print the vectors too for debugging
                string nxStr = f.Normal != null && f.Normal.Length == 2 ? f.Normal[0].ToString("F6") : "N/A";
                string nyStr = f.Normal != null && f.Normal.Length == 2 ? f.Normal[1].ToString("F6") : "N/A";
                string txStr = f.Tangent != null && f.Tangent.Length == 2 ? f.Tangent[0].ToString("F6") : "N/A";
                string tyStr = f.Tangent != null && f.Tangent.Length == 2 ? f.Tangent[1].ToString("F6") : "N/A";

                Console.WriteLine($"Face {f.Id}: L={f.Depth:F6}  t={f.Thickness:F6}  µ={mu:F6}  c={f.CohesionValue:F6} | N=({nxStr}, {nyStr}) T=({txStr}, {tyStr})");
            }
        }

    } // end of Program

} //end namespace
