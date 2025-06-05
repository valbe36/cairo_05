using System.Globalization;
using Gurobi;
using System.IO; 
using System.Collections.Generic;
using static InterlockingMasonryLocalForces.ProblemData;
using System.Reflection;
using static System.Runtime.InteropServices.JavaScript.JSType;
using System.ComponentModel;
using System.Data.Common;
using System.Diagnostics.Metrics;
using System.Numerics;
using System.Runtime.Intrinsics.X86;
using System.Xml.Linq;
using System;

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
        public int ExpectedNumVariables { get; set; }

        // Add this helper to track column indices
        private static Dictionary<(int FaceId, int VertexId), int> _columnMap;

        public static double[,] BuildEquilibriumMatrix(GeometryModel geometry, int numBlocks)
        {
            // Filter out only blocks with Id>0 in blockRowMap
            var nonSupportBlocks = geometry.Blocks.Values
            .Where(b => b.Id > 0)
            .OrderBy(b => b.Id)
            .ToList();
            Dictionary<int, int> blockRowMap = nonSupportBlocks
            .Select((b, idx) => new { b.Id, idx })
            .ToDictionary(x => x.Id, x => x.idx * 3);
            // Step 1: Create a column map for (faceId, vertexId) → column index
            Dictionary<(int FaceId, int VertexId), int> columnMap = new Dictionary<(int, int), int>();
            int colIdx = 0;
            foreach (var face in geometry.Faces.Values.OrderBy(f => f.Id))
            {
                foreach (int vId in face.VertexIds)
                {
                    columnMap[(face.Id, vId)] = colIdx++;
                }
            }

            int numRows = numBlocks * 3;
            int numCols = columnMap.Count * 2; // 2 variables (N,T) per (face,vertex)
            double[,] matrixA = new double[numRows, numCols];

            // Step 2: Fill matrix row by row
            foreach (var face in geometry.Faces.Values.OrderBy(f => f.Id))
            {
                bool isJSupport = face.BlockJ <= 0;
                bool isKSupport = face.BlockK <= 0;

                // If both J, K are supports, skip entirely
                if (isJSupport && isKSupport)
                    continue;

                double[] n = face.Normal;
                double[] t = face.Tangent;

                foreach (int vId in face.VertexIds)
                {
                    // locate the column for this face‐vertex
                    int colPairIndex = columnMap[(face.Id, vId)];
                    int colN = colPairIndex * 2;
                    int colT = colN + 1;

                    ContactPoint cp = geometry.Vertices[vId];

                    // If J is non‐support, add the “minus” contribution to blockJ
                    if (!isJSupport && blockRowMap.TryGetValue(face.BlockJ, out int jRow))
                    {
                        // local coords of cp relative to blockJ centroid
                        var bJ = geometry.Blocks[face.BlockJ];
                        double xRelJ = cp.X - bJ.CentroidX;
                        double yRelJ = cp.Y - bJ.CentroidY;

                        // Fx
                        matrixA[jRow, colN] -= n[0];
                        matrixA[jRow, colT] -= t[0];
                        // Fy
                        matrixA[jRow + 1, colN] -= n[1];
                        matrixA[jRow + 1, colT] -= t[1];
                        // M
                        matrixA[jRow + 2, colN] -= (xRelJ * n[1] - yRelJ * n[0]);
                        matrixA[jRow + 2, colT] -= (xRelJ * t[1] - yRelJ * t[0]);
                    }

                    // If K is non‐support, add the “plus” contribution to blockK
                    if (!isKSupport && blockRowMap.TryGetValue(face.BlockK, out int kRow))
                    {
                        var bK = geometry.Blocks[face.BlockK];
                        double xRelK = cp.X - bK.CentroidX;
                        double yRelK = cp.Y - bK.CentroidY;

                        matrixA[kRow, colN] += + n[0];
                        matrixA[kRow, colT] += +t[0];
                        matrixA[kRow + 1, colN] += +n[1];
                        matrixA[kRow + 1, colT] += + t[1];
                        matrixA[kRow + 2, colN] += + (xRelK * n[1] - yRelK * n[0]);
                        matrixA[kRow + 2, colT] += + (xRelK * t[1] - yRelK * t[0]);
                    }
                }
            } // end foreach face

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

        public double Thickness { get; set; }
        public double Length { get; set; }
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
            // 1) Collect face-vertex pairs in EXACTLY THE SAME ORDER as matrix construction
            faceVertexPairs = new List<FaceVertexPair>();

            // Use the SAME ordering as BuildEquilibriumMatrix
            foreach (var face in geometry.Faces.Values.OrderBy(f => f.Id))  // ← SAME as matrix!
            {
                foreach (int vId in face.VertexIds)
                {
                    faceVertexPairs.Add(new FaceVertexPair(face.Id, vId));
                }
            }
            pairIndexMap = new Dictionary<(int, int), int>(faceVertexPairs.Count);
            for (int j = 0; j < faceVertexPairs.Count; j++)
            {
                var p = faceVertexPairs[j];
                pairIndexMap[(p.FaceId, p.VertexId)] = j;
            }

            _geometry = geometry;

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

                    // *** ADD THIS LINE FOR VERIFICATION ***
                    VerifyMatrixVariableCorrespondence(geometry, data);

                    model.VerifyEquilibrium(data, faceVertexPairs, fAll, lambda, geometry);

                    // 5) Add friction & no-tension constraints
                    AddContactConstraints(model, geometry, data);

                    // 6)   Add face eccentricity constraints 
                     AddFaceEccConstraintsAroundV1(model, geometry, data);
                    //AddHybridEccConstraints(model, geometry, data);

                    // 7) Objective: maximize lambda
                    GRBLinExpr obj = 0.0;
                        obj.AddTerm(1.0, lambda);
                        model.SetObjective(obj, GRB.MAXIMIZE);
                        model.Write("debugModel.lp");
                        

                    model.Parameters.NumericFocus = 3;      // Maximum precision
                    model.Parameters.FeasibilityTol = 1e-9;  // Less aggressive; 1e-9 Tighter
                    model.Parameters.OptimalityTol = 1e-9;   // Less aggressive; 1e-9 Tighter
                    model.Parameters.BarConvTol = 1e-8;      // Barrier convergence
                    model.Parameters.MarkowitzTol = 0.01;    // Numerical stability
                    model.Parameters.ScaleFlag = 2;          // Automatic scaling
                    model.Parameters.Aggregate = 0;          // Disable presolve aggregation
                    model.Parameters.PreCrush = 1;           // Enable constraint matrix reduction
                    model.Parameters.Quad = 1;               // Use quadratic algorithm if applicable

                    // For debugging - disable some presolve steps
                    model.Parameters.Presolve = 1;           // Enable but conservative

                    // CRITICAL FOR QUADRATIC CONSTRAINTS:
                    model.Parameters.NonConvex = 2;         // Enable non-convex quadratic optimization
                    model.Parameters.Method = 2;            // Use barrier method
                    model.Parameters.Crossover = 0;         // Disable crossover
                    model.Parameters.PreQLinearize = 0;     // Don't linearize quadratic terms

                    // model.Parameters.MIPGap = 0.0005;
                    // model.Parameters.TimeLimit = 100; // Limit to 60 seconds
                    //model.Parameters.Quad = 1;       // Convex quadratic relaxation
            

                    // 8) Solve
                    model.Optimize();

                    //9 checks
                    if (model.Status == GRB.Status.OPTIMAL)
                    {
                        VerifyPostSolution(model, data, data.Mu);
                    }


                    SaveResultsToFile(model, @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\results_cairo.txt", data);
                    // 10) Print solution
                    PrintSolution(model);
                    }
                }
        }

        //  DEBUG  Add Debugging Method to Verify Matrix-Variable Correspondence
        private void VerifyMatrixVariableCorrespondence(GeometryModel geometry, ProblemData data)
        {
            Console.WriteLine("\n=== MATRIX-VARIABLE CORRESPONDENCE CHECK ===");

            // Build the same column map as in BuildEquilibriumMatrix
            Dictionary<(int FaceId, int VertexId), int> matrixColumnMap = new Dictionary<(int, int), int>();
            int matrixColIdx = 0;
            foreach (var face in geometry.Faces.Values.OrderBy(f => f.Id))
            {
                foreach (int vId in face.VertexIds)
                {
                    matrixColumnMap[(face.Id, vId)] = matrixColIdx++;
                }
            }

            // Compare with your variable ordering
            bool orderingMatch = true;
            for (int j = 0; j < faceVertexPairs.Count; j++)
            {
                var pair = faceVertexPairs[j];
                int expectedMatrixCol = matrixColumnMap[(pair.FaceId, pair.VertexId)];

                if (expectedMatrixCol != j)
                {
                    Console.WriteLine($"MISMATCH: Variable pair {j} (Face {pair.FaceId}, Vertex {pair.VertexId}) " +
                                    $"corresponds to matrix column {expectedMatrixCol}");
                    orderingMatch = false;
                }
            }

            if (orderingMatch)
            {
                Console.WriteLine("✓ Matrix columns and variables are correctly aligned");
            }
            else
            {
                Console.WriteLine("✗ CRITICAL ERROR: Matrix-variable ordering mismatch!");
                Console.WriteLine("This will cause incorrect solutions that worsen with friction.");
            }
        }

        // CRITICAL FIX 3: Enhanced Post-Solution Verification
        // CRITICAL FIX 3: Enhanced Post-Solution Verification
        private void VerifyPostSolution(GRBModel model, ProblemData data, double friction)
        {
            if (model.Status != GRB.Status.OPTIMAL) return;

            Console.WriteLine($"\n=== POST-SOLUTION VERIFICATION (μ = {friction}) ===");

            double lambdaVal = lambda.X;
            Console.WriteLine($"Lambda = {lambdaVal:F6}");

            // Manual equilibrium check for first block
            if (faceVertexPairs.Count > 0)
            {
                Console.WriteLine("\nManual equilibrium check for Block 1:");

                double sumFx = 0, sumFy = 0, sumM = 0;
                var block1 = _geometry.Blocks[1];

                // Check all forces affecting block 1
                for (int j = 0; j < faceVertexPairs.Count; j++)
                {
                    var pair = faceVertexPairs[j];
                    var face = _geometry.Faces[pair.FaceId];

                    // Only consider faces that touch block 1
                    if (face.BlockJ != 1 && face.BlockK != 1) continue;

                    double fn = fAll[2 * j].X;     // Normal force
                    double ft = fAll[2 * j + 1].X; // Tangential force

                    // Determine sign based on whether block 1 is J or K
                    double sign = (face.BlockJ == 1) ? -1.0 : +1.0;

                    // Force contributions
                    sumFx += sign * (fn * face.Normal[0] + ft * face.Tangent[0]);
                    sumFy += sign * (fn * face.Normal[1] + ft * face.Tangent[1]);

                    // Moment contributions
                    var vertex = _geometry.Vertices[pair.VertexId];
                    double xRel = vertex.X - block1.CentroidX;
                    double yRel = vertex.Y - block1.CentroidY;
                    double momentArm = xRel * face.Normal[1] - yRel * face.Normal[0];
                    sumM += sign * fn * momentArm;

                    momentArm = xRel * face.Tangent[1] - yRel * face.Tangent[0];
                    sumM += sign * ft * momentArm;
                }

                // Add external loads
                if (data.B != null && data.B.Length > 2)
                {
                    sumFx += data.B[0] * lambdaVal;
                    sumFy += data.B[1] * lambdaVal;
                    sumM += data.B[2] * lambdaVal;
                }

                // Add gravity
                if (data.G != null && data.G.Length > 2)
                {
                    sumFx += data.G[0];
                    sumFy += data.G[1];
                    sumM += data.G[2];
                }

                Console.WriteLine($"  Manual ΣFx = {sumFx:F6} (should ≈ 0)");
                Console.WriteLine($"  Manual ΣFy = {sumFy:F6} (should ≈ 0)");
                Console.WriteLine($"  Manual ΣM = {sumM:F6} (should ≈ 0)");

                if (Math.Abs(sumFx) > 1e-3 || Math.Abs(sumFy) > 1e-3 || Math.Abs(sumM) > 1e-3)
                {
                    Console.WriteLine("  ⚠️  EQUILIBRIUM VIOLATION DETECTED!");
                }
            }
            // Check friction constraint violations manually
            Console.WriteLine("\nManual friction constraint check:");
            foreach (var faceKvp in _geometry.Faces)
            {
                var face = faceKvp.Value;
                double mu = face.MuOverride ?? friction;
                double cohesion = face.CohesionValue;
                double area = face.Length * face.Thickness;

                // Sum forces on this face
                double totalNormal = 0, totalTangential = 0;
                foreach (int vId in face.VertexIds)
                {
                    int idx = GetPairIndex(face.Id, vId);
                    if (idx >= 0)
                    {
                        totalNormal += fAll[2 * idx].X;
                        totalTangential += fAll[2 * idx + 1].X;
                    }
                }

                double frictionLimit = mu * totalNormal + cohesion * area;
                double frictionUsage = Math.Abs(totalTangential) / (frictionLimit + 1e-12);

                Console.WriteLine($"  Face {face.Id}: |T|={Math.Abs(totalTangential):F3}, limit={frictionLimit:F3}, usage={frictionUsage:F3}");

                if (frictionUsage > 1.01) // Allow small tolerance
                {
                    Console.WriteLine($"    ⚠️  FRICTION VIOLATION: {frictionUsage:F3} > 1.0");
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




        // vertex based friction limit
        /*
        private void AddContactConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            foreach (var faceKvp in geometry.Faces)
            {
                Face face = faceKvp.Value;
                double mu = face.MuOverride ?? data.Mu;
                double c = face.CohesionValue;
                // For a 2-vertex face, each vertex gets half of the total face cohesion:
                double area = face.Thickness * face.Length;
                double faceCohesion = c * area; // total face cohesion
                double cPerVertex = faceCohesion / face.VertexIds.Count;  // e.g., 2.0 if face has 2 vertices

                foreach (int vId in face.VertexIds)
                {
                    int idx = GetPairIndex(face.Id, vId);
                    if (idx < 0) continue; // skip if not found

                    GRBVar fN = fAll[2 * idx];     // normal
                    GRBVar fT = fAll[2 * idx + 1]; // tangent

                    // 1) No tension
                    model.AddConstr(fN >= 0.0, $"NoTension_face{face.Id}_v{vId}");

                    // 2) Vertex friction (with partial cohesion)
                    if (mu <= 1e-12 && c <= 1e-12)
                    {
                        // frictionless => restrict tangential to zero
                        model.AddConstr(fT == 0.0, $"NoShear_face{face.Id}_v{vId}");
                    }
                    else
                    {
                        model.AddConstr(fT <= mu * fN + cPerVertex, $"FricPos_face{face.Id}_v{vId}");
                        model.AddConstr(fT >= -mu * fN - cPerVertex, $"FricNeg_face{face.Id}_v{vId}");
                    }
                }
            }
        }*/

        // face based friction limit

        /// For each pair j, friction constraints: -mu * fN_j <= fT_j <= mu * fN_j
        /// plus no tension fN_j >= 0 (already in the variable bound).
        private void AddContactConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
            {
                foreach (var fKvp in geometry.Faces)
                {
                    Face face = fKvp.Value;
                    double mu = face.MuOverride ?? data.Mu;
                    double c = face.CohesionValue;
                    double area = face.Thickness * face.Length;
                    double cohesionTerm = c * area;  // Full face cohesion

                    // 1. First add no-tension constraints for each vertex
                    foreach (int vId in face.VertexIds)
                    {
                        int idx = GetPairIndex(face.Id, vId);
                        if (idx < 0) continue;

                        GRBVar fN = fAll[2 * idx];
                        model.AddConstr(fN >= 0.0, $"NoTension_face{face.Id}_v{vId}");
                    }

                    // 2. Create face-level sum variables for normal and tangential forces
                    List<int> vertexIndices = new List<int>();
                    foreach (int vId in face.VertexIds)
                    {
                        int idx = GetPairIndex(face.Id, vId);
                        if (idx >= 0) vertexIndices.Add(idx);
                    }

                    if (vertexIndices.Count == 0) continue;

                    // Create variables for total normal and tangential force on face
                    GRBVar faceNormal = model.AddVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, $"faceN_{face.Id}");
                    GRBVar faceTangent = model.AddVar(-GRB.INFINITY, GRB.INFINITY, 0.0, GRB.CONTINUOUS, $"faceT_{face.Id}");

                    // Sum up all vertex normal forces
                    GRBLinExpr normalSum = new GRBLinExpr();
                    foreach (int idx in vertexIndices)
                    {
                        normalSum.AddTerm(1.0, fAll[2 * idx]);
                    }
                    model.AddConstr(faceNormal == normalSum, $"SumN_face{face.Id}");

                    // Sum up all vertex tangential forces
                    GRBLinExpr tangentSum = new GRBLinExpr();
                    foreach (int idx in vertexIndices)
                    {
                        tangentSum.AddTerm(1.0, fAll[2 * idx + 1]);
                    }
                    model.AddConstr(faceTangent == tangentSum, $"SumT_face{face.Id}");

                    // 3. Apply face-level friction-cohesion constraints
                    if (mu <= 1e-12 && c <= 1e-12)
                    {
                        // Frictionless, no cohesion => no shear
                        model.AddConstr(faceTangent == 0.0, $"NoShear_face{face.Id}");
                    }
                    else
                    {
                        // Face-level constraints with full cohesion
                        model.AddConstr(faceTangent <= mu * faceNormal + cohesionTerm, $"FricPos_face{face.Id}");
                        model.AddConstr(faceTangent >= -mu * faceNormal - cohesionTerm, $"FricNeg_face{face.Id}");
                    }
                }
            }
          
       
        private void AddFaceEccConstraintsAroundV1(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            double sigmaC = data.SigmaC;

            Console.WriteLine("\n--- Face Eccentricity Constraints (Moment about v1) ---");

            foreach (var fKvp in geometry.Faces)
            {
                Face face = fKvp.Value;
                // We assume each face in 2D has exactly 2 vertices.
                if (face.VertexIds.Count != 2) continue;

                // Get vertex IDs
                int vertexId1 = face.VertexIds[0];
                int vertexId2 = face.VertexIds[1];

                // Get vertex coordinates
                ContactPoint vertex1 = _geometry.Vertices[vertexId1];
                ContactPoint vertex2 = _geometry.Vertices[vertexId2];

                // Get pair indices (for normal forces)
                int idx1 = GetPairIndex(face.Id, vertexId1);
                int idx2 = GetPairIndex(face.Id, vertexId2);

                if (idx1 < 0 || idx2 < 0)
                {
                    Console.WriteLine($"Skipping ecc constraints for Face {face.Id}.");
                    continue;
                }

                // Local normal vars at v1 & v2
                GRBVar fn1 = fAll[2 * idx1];   // normal at (face, v1)
                GRBVar fn2 = fAll[2 * idx2];   // normal at (face, v2)

                // 1) fnk = fn1 + fn2
                // GRBVar fnk = model.AddVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, $"fnk_face{face.Id}");// a problem in quadratic
                double realisticMaxForce = sigmaC * face.Length * face.Thickness;  // Physical limit
                GRBVar fnk = model.AddVar(0.0, realisticMaxForce, 0.0, GRB.CONTINUOUS, $"fnk_face{face.Id}");  // GOOD: Realistic bound
                {
                    GRBLinExpr sumExpr = new GRBLinExpr();
                    sumExpr.AddTerm(1.0, fn1);
                    sumExpr.AddTerm(1.0, fn2);
                    model.AddConstr(fnk == sumExpr, $"Def_fnk_face{face.Id}");
                }

                // 2) Distance from v1 to v2 (Lk). 
                //    We'll just use the Euclidean distance so that Lk > 0.
                //    If it's purely horizontal or vertical, this is still fine.
                double dx = vertex2.X - vertex1.X;
                double dy = vertex2.Y - vertex1.Y;
                double Lk = Math.Sqrt(dx * dx + dy * dy);

                // 3) Eccentricity from v1 (call it eK), bounding it from 0 to Lk 
                //    if we want no tension at either edge. 
                //    If you might allow partial tension, you could do a wider range: 
                //    e.g. eK in [-0.5*Lk, 1.5*Lk], etc.
                GRBVar eK = model.AddVar(0.0, Lk, 0.0, GRB.CONTINUOUS, $"eK_face{face.Id}");

                faceEccVars[face.Id] = eK;  // (Optional) Keep a reference for debugging

                // 4) "Moment about v1":  fnk * eK = fn2 * Lk
                //    i.e. the total normal force times its lever arm about v1 
                //    must match what fn2 * Lk alone would produce at v2.
                {
                    GRBQuadExpr mq = new GRBQuadExpr();
                    //   +1 * (fnk * eK)
                    mq.AddTerm(1.0, fnk, eK);
                    //   - Lk * fn2
                    mq.AddTerm(-Lk, fn2);

                    model.AddQConstr(mq == 0.0, $"MomentEq_face{face.Id}");
                }

                // 5) If sigmaC > 0, add compressive/bending strength constraints.
                double t_k = face.Thickness;
                if (sigmaC > 1e-9 && t_k > 1e-9)
                {
                    // denom = 2 * sigmaC * thickness
                    double denom = 2.0 * sigmaC * t_k;

                    // (a)  eK + fnk / denom <= Lk
                    {
                        GRBLinExpr lhsUp = 0.0;
                        lhsUp.AddTerm(1.0, eK);
                        lhsUp.AddTerm(1.0 / denom, fnk);
                        model.AddConstr(lhsUp <= Lk, $"StrengthUp_face{face.Id}");
                    }

                    // (b)  eK - fnk / denom >= 0
                    {
                        GRBLinExpr lhsLo = 0.0;
                        lhsLo.AddTerm(1.0, eK);
                        lhsLo.AddTerm(-1.0 / denom, fnk);
                        model.AddConstr(lhsLo >= 0.0, $"StrengthLo_face{face.Id}");
                    }
                }
            }
        }

        /// <summary>
         
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

        private void SaveResultsToFile(GRBModel model, string resultsFilePath, ProblemData data)
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


                // 5) Shear limits at the FACE level (not vertex level)
                writer.WriteLine("\n--- Face-level Shear Check ---");

                foreach (var kvp in faceTotalForces)
                {
                    int faceId = kvp.Key;
                    double totalNormal = kvp.Value.fnSum;
                    double totalShear = kvp.Value.ftSum;

                    // Get face data
                    Face face = _geometry.Faces[faceId];
                    double mu = face.MuOverride ?? data.Mu;
                    double cohesion = face.CohesionValue;
                    double area = face.Length * face.Thickness;
                    double shearLimit = mu * totalNormal + cohesion * area;

                    // Shear "usage" ratio
                    double tiny = 1e-9;  // small offset to avoid divide-by-zero
                    double usageShear = Math.Abs(totalShear + tiny) / (shearLimit + tiny);

                    // If usage is near 1.0, you're "on" the friction limit
                    string shearStatus;
                    double usageTolerance = 0.98;  // or pick something else like 0.99
                    if (usageShear >= usageTolerance)
                        shearStatus = "CRITICAL";
                    else
                        shearStatus = "OK";

                    writer.WriteLine($"Face {faceId} Shear: ratio={usageShear:F3}, status={shearStatus} " +
$"(|shear|={Math.Abs(totalShear):F3}, limit={shearLimit:F3})");
                }

                // 6) Bending/Compression Usage Ratio
                writer.WriteLine("\n--- Face Eccentricity Check (Bending/Compression) ---");
                foreach (var kvp in faceEccVars)
                {
                    int faceId = kvp.Key;
                    Face face = _geometry.Faces[faceId];
                    double eK = kvp.Value.X;  // or faceEccVars[faceId].X
                    double totalNormal = faceTotalForces[faceId].fnSum;
                    double length = face.Length;
                    double thickness = face.Thickness;
                    double sigmaC = data.SigmaC;

                    // a) Determine compressed “width” in-plane
                    //    eKeff = eK if eK < length/2, else (length - eK).
                    double eKeff = (eK <= (length / 2.0))
                        ? eK
                        : (length - eK);

                    // b) The compressed cross-section area = t × 2 × eKeff
                    //    (assuming eK is always in [0, length], so 2 eK or 2(length-eK))
                    double areaPrime = 2.0 * thickness * eKeff;

                    // c) Compute the actual stress = (totalNormal / areaPrime),
                    //    with some fallback if areaPrime ~ 0
                    double tiny = 1e-9;
                    double actualStress = 0.0;
                    if (areaPrime > tiny)
                    {
                        actualStress = totalNormal / areaPrime;
                    }
                    else
                    {
                        // If eK=0 or eK=length=> areaPrime=0 => set stress=0 or skip
                        // For example, if eK=0 => entire face is compressed at the vertex
                        // or if eK=length => block is compressed at the other vertex.
                        // It's your call how to handle physically.
                        actualStress = 0.0;
                    }

                    // d) Bending/compression usage ratio
                    //    usageBC = (actualStress / sigmaC), forced never to exceed 1.0 by constraints
                    double usageBC = (actualStress + tiny) / (sigmaC + tiny);

                    // Decide how close we are to the limit:
                    double usageTolerance = 0.98;  // or 0.99, etc.
                    string eccStatus = (usageBC >= usageTolerance) ? "CRITICAL" : "OK";

                    // e) Write it out:
                    writer.WriteLine($"Face {faceId}: eK={eK:F3}, eKeff={eKeff:F3}, totalN={totalNormal:F3}, " +
                                     $"area'={areaPrime:F3}, actualStress={actualStress:F3}, usage={usageBC:F3}, status={eccStatus}");
                }

                //7 save used tangents and vectors
                writer.WriteLine("\n--- Face Normal/Tangent Vectors ---");
                foreach (var fEntry in _geometry.Faces)
                {
                    int faceId = fEntry.Key;
                    Face face = fEntry.Value;
                    double nx = face.Normal[0];
                    double ny = face.Normal[1];
                    double tx = face.Tangent[0];
                    double ty = face.Tangent[1];
                    // We'll store coordinates of "vertex 1" (the second vertex in face.VertexIds).
                    // Default to (0,0) if something is missing.
                    double coordX = 0.0;
                    double coordY = 0.0;

                    // If the face actually has 2 vertices, get vertex #1:
                    if (face.VertexIds.Count == 2)
                    {
                        int vId0 = face.VertexIds[0];
                        if (_geometry.Vertices.ContainsKey(vId0))
                        {
                            ContactPoint cp = _geometry.Vertices[vId0];
                            coordX = cp.X;
                            coordY = cp.Y;
                        }
                    }
                    // Print everything
                    writer.WriteLine(
                        $"Face {faceId}: " +
                        $"Vertex0=({coordX:F3}, {coordY:F3},0),  " +
                        $"Normal=({nx:F3}, {ny:F3}),  " +
                        $"Tangent=({tx:F3}, {ty:F3})"
                    );
                }


            }
            Console.WriteLine($"Results saved to {resultsFilePath}");
        }

    } //end public class LocalOptimizer


    public class EquilibriumMatrixVerifier
    {
        /// <summary>
        /// Comprehensive verification of the equilibrium matrix construction
        /// </summary>
        public static void VerifyEquilibriumMatrix(ProblemData data,
                                                 List<FaceVertexPair> faceVertexPairs,
                                                 GRBVar[] fAll,
                                                 GRBVar lambda,
                                                 GeometryModel geometry)
        {
            Console.WriteLine("=== EQUILIBRIUM MATRIX VERIFICATION ===");

            // Basic dimension checks
            VerifyMatrixDimensions(data, faceVertexPairs, fAll);

            // Structural analysis
            AnalyzeMatrixStructure(data);

            // Force direction verification
            VerifyForceDirections(data, faceVertexPairs, geometry);

            // Load vector verification
            VerifyLoadVectors(data);

            // Detailed matrix content analysis
            AnalyzeMatrixContent(data);

            // Row-wise equilibrium check
            VerifyEquilibriumStructure(data, geometry);

            // Critical checks specific to your implementation
            VerifyNormalTangentVectors(geometry);
            VerifySignConventions(data, faceVertexPairs, geometry);
        }

        /// <summary>
        /// Verify matrix dimensions match expected structure
        /// </summary>
        private static void VerifyMatrixDimensions(ProblemData data,
                                                 List<FaceVertexPair> faceVertexPairs,
                                                 GRBVar[] fAll)
        {
            Console.WriteLine("\n--- DIMENSION VERIFICATION ---");

            int expectedCols = 2 * faceVertexPairs.Count; // 2 forces per face-vertex pair
            int actualCols = data.NumCols;
            int actualRows = data.NumRows;

            Console.WriteLine($"Face-vertex pairs: {faceVertexPairs.Count}");
            Console.WriteLine($"Expected columns: {expectedCols}");
            Console.WriteLine($"Actual columns: {actualCols}");
            Console.WriteLine($"Matrix rows: {actualRows}");
            Console.WriteLine($"Force variables: {fAll.Length}");

            if (expectedCols != actualCols)
            {
                Console.WriteLine($"ERROR: Column mismatch! Expected {expectedCols}, got {actualCols}");
            }

            if (actualCols != fAll.Length)
            {
                Console.WriteLine($"ERROR: Variable count mismatch! Matrix has {actualCols} columns, {fAll.Length} variables");
            }

            // Expected rows should be: 2 * numVertices (Fx, Fy) + numVertices (moments) = 3 * numVertices
            // Or it could be different depending on your formulation
            Console.WriteLine($"Rows should typically be 3 * num_vertices for 2D problems");
        }

        /// <summary>
        /// Analyze the overall structure of matrix A
        /// </summary>
        private static void AnalyzeMatrixStructure(ProblemData data)
        {
            Console.WriteLine("\n--- MATRIX STRUCTURE ANALYSIS ---");

            int rows = data.NumRows;
            int cols = data.NumCols;
            int nonZeros = 0;
            double maxValue = 0;
            double minValue = double.MaxValue;

            // Count non-zeros and analyze values
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    double val = Math.Abs(data.MatrixA[i, j]);
                    if (val > 1e-15)
                    {
                        nonZeros++;
                        maxValue = Math.Max(maxValue, val);
                        minValue = Math.Min(minValue, val);
                    }
                }
            }

            double sparsity = (double)nonZeros / (rows * cols);

            Console.WriteLine($"Non-zero entries: {nonZeros} / {rows * cols}");
            Console.WriteLine($"Sparsity: {sparsity:F4}");
            Console.WriteLine($"Value range: [{minValue:E3}, {maxValue:E3}]");
            Console.WriteLine($"Condition indicator (max/min): {maxValue / minValue:E3}");

            if (maxValue / minValue > 1e12)
            {
                Console.WriteLine("WARNING: Very large condition number - potential numerical issues!");
            }
        }

        /// <summary>
        /// Verify force directions are correctly represented in matrix A
        /// </summary>
        private static void VerifyForceDirections(ProblemData data,
                                                List<FaceVertexPair> faceVertexPairs,
                                                GeometryModel geometry)
        {
            Console.WriteLine("\n--- FORCE DIRECTION VERIFICATION ---");

            // For each face-vertex pair, check if the corresponding columns make sense
            for (int pairIdx = 0; pairIdx < Math.Min(5, faceVertexPairs.Count); pairIdx++) // Check first 5 pairs
            {
                int normalColIdx = 2 * pairIdx;     // f_n column
                int tangentialColIdx = 2 * pairIdx + 1; // f_t column
                var pair = faceVertexPairs[pairIdx];

                Console.WriteLine($"\nPair {pairIdx} (Face {pair.FaceId}, Vertex {pair.VertexId}):");

                // Get the face to access normal/tangent vectors
                if (geometry.Faces.TryGetValue(pair.FaceId, out Face face))
                {
                    Console.WriteLine($"  Face normal: ({face.Normal[0]:F6}, {face.Normal[1]:F6})");
                    Console.WriteLine($"  Face tangent: ({face.Tangent[0]:F6}, {face.Tangent[1]:F6})");

                    // Verify unit vectors
                    double nMag = Math.Sqrt(face.Normal[0] * face.Normal[0] + face.Normal[1] * face.Normal[1]);
                    double tMag = Math.Sqrt(face.Tangent[0] * face.Tangent[0] + face.Tangent[1] * face.Tangent[1]);

                    Console.WriteLine($"  Normal magnitude: {nMag:F6}");
                    Console.WriteLine($"  Tangent magnitude: {tMag:F6}");

                    if (Math.Abs(nMag - 1.0) > 0.01)
                        Console.WriteLine($"  WARNING: Normal vector is not unit length!");
                    if (Math.Abs(tMag - 1.0) > 0.01)
                        Console.WriteLine($"  WARNING: Tangent vector is not unit length!");

                    // Check orthogonality
                    double dot = face.Normal[0] * face.Tangent[0] + face.Normal[1] * face.Tangent[1];
                    Console.WriteLine($"  Normal⋅Tangent: {dot:F6}");
                    if (Math.Abs(dot) > 0.01)
                        Console.WriteLine($"  WARNING: Normal and tangent are not orthogonal!");
                }

                // Analyze normal force column
                AnalyzeForceColumn(data, normalColIdx, "Normal");

                // Analyze tangential force column  
                AnalyzeForceColumn(data, tangentialColIdx, "Tangential");
            }
        }

        /// <summary>
        /// Analyze a single force column in the matrix
        /// </summary>
        private static void AnalyzeForceColumn(ProblemData data, int colIdx, string forceType)
        {
            if (colIdx >= data.NumCols) return;

            Console.WriteLine($"  {forceType} force column {colIdx}:");

            // For equilibrium matrices, each column represents the TOTAL effect of a unit force
            // on ALL equilibrium equations, not a unit vector itself

            int numBlocks = data.NumRows / 3; // Assuming 3 equations per block (Fx, Fy, M)

            double totalFxContrib = 0, totalFyContrib = 0, totalMContrib = 0;

            // Sum contributions to all Fx equations (every 3rd row starting from 0)
            for (int blockIdx = 0; blockIdx < numBlocks; blockIdx++)
            {
                int fxRow = blockIdx * 3;
                if (fxRow < data.NumRows)
                    totalFxContrib += data.MatrixA[fxRow, colIdx];
            }

            // Sum contributions to all Fy equations (every 3rd row starting from 1)
            for (int blockIdx = 0; blockIdx < numBlocks; blockIdx++)
            {
                int fyRow = blockIdx * 3 + 1;
                if (fyRow < data.NumRows)
                    totalFyContrib += data.MatrixA[fyRow, colIdx];
            }

            // Sum contributions to all Moment equations (every 3rd row starting from 2)
            for (int blockIdx = 0; blockIdx < numBlocks; blockIdx++)
            {
                int mRow = blockIdx * 3 + 2;
                if (mRow < data.NumRows)
                    totalMContrib += Math.Abs(data.MatrixA[mRow, colIdx]);
            }

            Console.WriteLine($"    Total Fx contributions: {totalFxContrib:F6}");
            Console.WriteLine($"    Total Fy contributions: {totalFyContrib:F6}");
            Console.WriteLine($"    Total |Moment| contributions: {totalMContrib:F6}");

            // For internal forces, total Fx and Fy should be ≈0 (Newton's 3rd law)
            // For boundary forces, they may be non-zero
            if (Math.Abs(totalFxContrib) < 1e-6 && Math.Abs(totalFyContrib) < 1e-6)
            {
                Console.WriteLine($"    ✓ Internal force (sums to zero)");
            }
            else
            {
                Console.WriteLine($"    → Boundary force or force imbalance");
            }
        }

        /// <summary>
        /// Check if normal and tangential force vectors are orthogonal
        /// </summary>
        private static void CheckOrthogonality(ProblemData data, int normalCol, int tangentialCol)
        {
            if (normalCol >= data.NumCols || tangentialCol >= data.NumCols) return;

            double dotProduct = 0;
            int numVertices = data.NumRows / 3;

            // Compute dot product of force direction vectors (Fx and Fy components only)
            for (int i = 0; i < Math.Min(2 * numVertices, data.NumRows); i++)
            {
                dotProduct += data.MatrixA[i, normalCol] * data.MatrixA[i, tangentialCol];
            }

            Console.WriteLine($"    Normal⋅Tangential dot product: {dotProduct:F6}");

            if (Math.Abs(dotProduct) > 0.1)
            {
                Console.WriteLine($"    WARNING: Normal and tangential vectors are not orthogonal!");
            }
        }

        /// <summary>
        /// Verify load vectors G and B
        /// </summary>
        private static void VerifyLoadVectors(ProblemData data)
        {
            Console.WriteLine("\n--- LOAD VECTOR VERIFICATION ---");

            if (data.G != null)
            {
                Console.WriteLine($"Gravity vector G: length {data.G.Length}");
                double gMagnitude = Math.Sqrt(data.G.Sum(x => x * x));
                Console.WriteLine($"  Magnitude: {gMagnitude:F6}");
                Console.WriteLine($"  Max component: {data.G.Max():F6}");
                Console.WriteLine($"  Min component: {data.G.Min():F6}");
            }

            if (data.B != null)
            {
                Console.WriteLine($"Applied load vector B: length {data.B.Length}");
                double bMagnitude = Math.Sqrt(data.B.Sum(x => x * x));
                Console.WriteLine($"  Magnitude: {bMagnitude:F6}");
                Console.WriteLine($"  Max component: {data.B.Max():F6}");
                Console.WriteLine($"  Min component: {data.B.Min():F6}");
            }

            // Check if G and B have same length as matrix rows
            if (data.G != null && data.G.Length != data.NumRows)
            {
                Console.WriteLine($"ERROR: G vector length {data.G.Length} ≠ matrix rows {data.NumRows}");
            }

            if (data.B != null && data.B.Length != data.NumRows)
            {
                Console.WriteLine($"ERROR: B vector length {data.B.Length} ≠ matrix rows {data.NumRows}");
            }
        }

        /// <summary>
        /// Detailed analysis of matrix content
        /// </summary>
        private static void AnalyzeMatrixContent(ProblemData data)
        {
            Console.WriteLine("\n--- DETAILED MATRIX CONTENT ---");

            // Show first few rows and columns
            int showRows = Math.Min(6, data.NumRows);
            int showCols = Math.Min(8, data.NumCols);

            Console.WriteLine("First few matrix entries:");
            Console.Write("Row\\Col".PadRight(8));
            for (int j = 0; j < showCols; j++)
            {
                Console.Write($"{j,10}");
            }
            Console.WriteLine();

            for (int i = 0; i < showRows; i++)
            {
                Console.Write($"{i,6}: ");
                for (int j = 0; j < showCols; j++)
                {
                    Console.Write($"{data.MatrixA[i, j],10:F3}");
                }
                Console.WriteLine();
            }

            // Check for suspicious patterns
            CheckForSuspiciousPatterns(data);
        }

        /// <summary>
        /// Check for common matrix construction errors
        /// </summary>
        private static void CheckForSuspiciousPatterns(ProblemData data)
        {
            Console.WriteLine("\n--- SUSPICIOUS PATTERN DETECTION ---");

            int zeroRows = 0;
            int zeroCols = 0;

            // Check for zero rows
            for (int i = 0; i < data.NumRows; i++)
            {
                bool isZeroRow = true;
                for (int j = 0; j < data.NumCols; j++)
                {
                    if (Math.Abs(data.MatrixA[i, j]) > 1e-15)
                    {
                        isZeroRow = false;
                        break;
                    }
                }
                if (isZeroRow)
                {
                    Console.WriteLine($"WARNING: Row {i} is all zeros!");
                    zeroRows++;
                }
            }

            // Check for zero columns
            for (int j = 0; j < data.NumCols; j++)
            {
                bool isZeroCol = true;
                for (int i = 0; i < data.NumRows; i++)
                {
                    if (Math.Abs(data.MatrixA[i, j]) > 1e-15)
                    {
                        isZeroCol = false;
                        break;
                    }
                }
                if (isZeroCol)
                {
                    Console.WriteLine($"WARNING: Column {j} is all zeros!");
                    zeroCols++;
                }
            }

            Console.WriteLine($"Zero rows: {zeroRows}, Zero columns: {zeroCols}");

            // Check for identical rows (might indicate duplicated constraints)
            CheckForIdenticalRows(data);
        }

        /// <summary>
        /// Check for identical or nearly identical rows
        /// </summary>
        private static void CheckForIdenticalRows(ProblemData data)
        {
            for (int i = 0; i < data.NumRows - 1; i++)
            {
                for (int j = i + 1; j < data.NumRows; j++)
                {
                    double diff = 0;
                    for (int k = 0; k < data.NumCols; k++)
                    {
                        diff += Math.Abs(data.MatrixA[i, k] - data.MatrixA[j, k]);
                    }

                    if (diff < 1e-12)
                    {
                        Console.WriteLine($"WARNING: Rows {i} and {j} are nearly identical!");
                    }
                }
            }
        }

        /// <summary>
        /// Verify that equilibrium equations have the expected structure
        /// </summary>
        private static void VerifyEquilibriumStructure(ProblemData data, GeometryModel geometry)
        {
            Console.WriteLine("\n--- EQUILIBRIUM STRUCTURE VERIFICATION ---");

            // For a masonry arch, expect:
            // - Force equilibrium equations (ΣFx = 0, ΣFy = 0) 
            // - Moment equilibrium equations (ΣM = 0)

            int numVertices = data.NumRows / 3; // Assuming 3 equations per vertex

            Console.WriteLine($"Assuming {numVertices} vertices with 3 equations each");

            // Check Fx equilibrium equations (should sum all x-components of forces)
            for (int vertexIdx = 0; vertexIdx < Math.Min(3, numVertices); vertexIdx++)
            {
                int rowIdx = vertexIdx; // Fx equation for this vertex

                double sumNormalX = 0;
                double sumTangentialX = 0;

                // Sum contributions from all force pairs
                for (int pairIdx = 0; pairIdx < data.NumCols / 2; pairIdx++)
                {
                    sumNormalX += data.MatrixA[rowIdx, 2 * pairIdx];     // Normal force X component
                    sumTangentialX += data.MatrixA[rowIdx, 2 * pairIdx + 1]; // Tangential force X component
                }

                Console.WriteLine($"Vertex {vertexIdx} Fx equation:");
                Console.WriteLine($"  Sum of normal X components: {sumNormalX:F6}");
                Console.WriteLine($"  Sum of tangential X components: {sumTangentialX:F6}");

                // For global equilibrium, these should sum to specific values depending on the problem
            }
        }

        /// <summary>
        /// Verify normal and tangent vectors in the geometry are properly computed
        /// </summary>
        private static void VerifyNormalTangentVectors(GeometryModel geometry)
        {
            Console.WriteLine("\n--- NORMAL/TANGENT VECTOR VERIFICATION ---");

            foreach (var faceKvp in geometry.Faces)
            {
                var face = faceKvp.Value;
                Console.WriteLine($"Face {face.Id}:");

                // Check if normal and tangent are unit vectors
                double nMag = Math.Sqrt(face.Normal[0] * face.Normal[0] + face.Normal[1] * face.Normal[1]);
                double tMag = Math.Sqrt(face.Tangent[0] * face.Tangent[0] + face.Tangent[1] * face.Tangent[1]);

                Console.WriteLine($"  Normal: ({face.Normal[0]:F6}, {face.Normal[1]:F6}), mag = {nMag:F6}");
                Console.WriteLine($"  Tangent: ({face.Tangent[0]:F6}, {face.Tangent[1]:F6}), mag = {tMag:F6}");

                if (Math.Abs(nMag - 1.0) > 0.01)
                    Console.WriteLine($"  ERROR: Normal vector magnitude {nMag:F6} ≠ 1.0");
                if (Math.Abs(tMag - 1.0) > 0.01)
                    Console.WriteLine($"  ERROR: Tangent vector magnitude {tMag:F6} ≠ 1.0");

                // Check orthogonality
                double dot = face.Normal[0] * face.Tangent[0] + face.Normal[1] * face.Tangent[1];
                Console.WriteLine($"  Dot product: {dot:F6}");
                if (Math.Abs(dot) > 0.01)
                    Console.WriteLine($"  ERROR: Vectors not orthogonal! Dot = {dot:F6}");

                // Verify face geometry makes sense
                if (face.VertexIds.Count == 2)
                {
                    var v1 = geometry.Vertices[face.VertexIds[0]];
                    var v2 = geometry.Vertices[face.VertexIds[1]];

                    // Face vector (v1 to v2)
                    double dx = v2.X - v1.X;
                    double dy = v2.Y - v1.Y;
                    double faceLength = Math.Sqrt(dx * dx + dy * dy);

                    Console.WriteLine($"  Face length: {faceLength:F6} (stored: {face.Length:F6})");

                    if (Math.Abs(faceLength - face.Length) > 0.01)
                        Console.WriteLine($"  WARNING: Computed length {faceLength:F6} ≠ stored length {face.Length:F6}");

                    // Check if tangent aligns with face direction
                    double[] faceDir = { dx / faceLength, dy / faceLength };
                    double tangentDot = Math.Abs(faceDir[0] * face.Tangent[0] + faceDir[1] * face.Tangent[1]);
                    Console.WriteLine($"  Tangent-Face alignment: {tangentDot:F6}");

                    if (tangentDot < 0.99)
                        Console.WriteLine($"  WARNING: Tangent not aligned with face direction!");
                }
            }
        }

        /// <summary>
        /// Verify sign conventions in the equilibrium matrix
        /// </summary>
        private static void VerifySignConventions(ProblemData data, List<FaceVertexPair> faceVertexPairs, GeometryModel geometry)
        {
            Console.WriteLine("\n--- SIGN CONVENTION VERIFICATION ---");

            // Check a few specific cases to verify signs
            for (int pairIdx = 0; pairIdx < Math.Min(3, faceVertexPairs.Count); pairIdx++)
            {
                var pair = faceVertexPairs[pairIdx];
                if (!geometry.Faces.TryGetValue(pair.FaceId, out Face face)) continue;
                if (!geometry.Vertices.TryGetValue(pair.VertexId, out ContactPoint vertex)) continue;

                Console.WriteLine($"\nChecking pair {pairIdx}: Face {pair.FaceId}, Vertex {pair.VertexId}");

                int colN = 2 * pairIdx;
                int colT = 2 * pairIdx + 1;

                // Check blocks J and K
                Console.WriteLine($"  Face connects Block {face.BlockJ} to Block {face.BlockK}");

                // For each block this face touches, verify sign consistency
                CheckBlockContribution(data, colN, colT, face.BlockJ, face, vertex, geometry, "BlockJ", true);
                CheckBlockContribution(data, colN, colT, face.BlockK, face, vertex, geometry, "BlockK", false);
            }
        }

        /// <summary>
        /// Check the contribution of a face-vertex pair to a specific block's equilibrium
        /// </summary>
        private static void CheckBlockContribution(ProblemData data, int colN, int colT, int blockId,
                                                 Face face, ContactPoint vertex, GeometryModel geometry,
                                                 string blockName, bool isBlockJ)
        {
            if (blockId <= 0) return; // Skip support blocks

            if (!geometry.Blocks.TryGetValue(blockId, out Block block)) return;

            // Find this block's row indices (assuming 3 equations per block)
            var nonSupportBlocks = geometry.Blocks.Values
                .Where(b => b.Id > 0)
                .OrderBy(b => b.Id)
                .ToList();

            var blockRowMap = nonSupportBlocks
                .Select((b, idx) => new { b.Id, idx })
                .ToDictionary(x => x.Id, x => x.idx * 3);

            if (!blockRowMap.TryGetValue(blockId, out int baseRow)) return;

            Console.WriteLine($"    {blockName} (ID {blockId}) base row: {baseRow}");

            // Get matrix coefficients for this block
            double fxN = data.MatrixA[baseRow, colN];       // Fx from normal force
            double fxT = data.MatrixA[baseRow, colT];       // Fx from tangential force
            double fyN = data.MatrixA[baseRow + 1, colN];   // Fy from normal force
            double fyT = data.MatrixA[baseRow + 1, colT];   // Fy from tangential force
            double mN = data.MatrixA[baseRow + 2, colN];    // Moment from normal force
            double mT = data.MatrixA[baseRow + 2, colT];    // Moment from tangential force

            Console.WriteLine($"    Matrix coeffs: FxN={fxN:F3}, FxT={fxT:F3}, FyN={fyN:F3}, FyT={fyT:F3}");
            Console.WriteLine($"    Moment coeffs: MN={mN:F3}, MT={mT:F3}");

            // Verify these match expected values based on normal/tangent directions
            double expectedFxN = isBlockJ ? -face.Normal[0] : +face.Normal[0];
            double expectedFxT = isBlockJ ? -face.Tangent[0] : +face.Tangent[0];
            double expectedFyN = isBlockJ ? -face.Normal[1] : +face.Normal[1];
            double expectedFyT = isBlockJ ? -face.Tangent[1] : +face.Tangent[1];

            Console.WriteLine($"    Expected: FxN={expectedFxN:F3}, FxT={expectedFxT:F3}, FyN={expectedFyN:F3}, FyT={expectedFyT:F3}");

            // Check moment arms
            double xRel = vertex.X - block.CentroidX;
            double yRel = vertex.Y - block.CentroidY;
            double expectedMN = isBlockJ ? -(xRel * face.Normal[1] - yRel * face.Normal[0]) : +(xRel * face.Normal[1] - yRel * face.Normal[0]);
            double expectedMT = isBlockJ ? -(xRel * face.Tangent[1] - yRel * face.Tangent[0]) : +(xRel * face.Tangent[1] - yRel * face.Tangent[0]);

            Console.WriteLine($"    Expected moments: MN={expectedMN:F3}, MT={expectedMT:F3}");
            Console.WriteLine($"    Lever arm: ({xRel:F3}, {yRel:F3})");

            // Check for discrepancies
            if (Math.Abs(fxN - expectedFxN) > 1e-6) Console.WriteLine($"    ERROR: FxN mismatch!");
            if (Math.Abs(fxT - expectedFxT) > 1e-6) Console.WriteLine($"    ERROR: FxT mismatch!");
            if (Math.Abs(fyN - expectedFyN) > 1e-6) Console.WriteLine($"    ERROR: FyN mismatch!");
            if (Math.Abs(fyT - expectedFyT) > 1e-6) Console.WriteLine($"    ERROR: FyT mismatch!");
            if (Math.Abs(mN - expectedMN) > 1e-6) Console.WriteLine($"    ERROR: MN mismatch!");
            if (Math.Abs(mT - expectedMT) > 1e-6) Console.WriteLine($"    ERROR: MT mismatch!");
        }
    }

    // Extension method to call the verification
    public static class ModelVerificationExtensions
    {
        public static void VerifyEquilibrium(this GRBModel model, ProblemData data,
                                           List<FaceVertexPair> faceVertexPairs,
                                           GRBVar[] fAll, GRBVar lambda, GeometryModel geometry)
        {
            EquilibriumMatrixVerifier.VerifyEquilibriumMatrix(data, faceVertexPairs, fAll, lambda, geometry);
        }
    }









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
                LoadAllData(@"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\/data_pseudoparallel_friction_0e4.txt"   //data_pseudoparallel_friction_0e4
                , geometry, data);
                ComputeFaceNormalsFromGeometry(geometry);
                // validate data
                ValidateGeometryModel(geometry);

                // Build equilibrium matrix programmatically
                data.NumBlocks = geometry.Blocks.Values.Count(b => b.Id > 0); // Exclude supports
                data.MatrixA = ProblemData.BuildEquilibriumMatrix(geometry, data.NumBlocks);

                // Now check if the actual number of columns matches expected
                if (data.MatrixA != null)
                {
                    int actualCols = data.MatrixA.GetLength(1);
                    if (data.ExpectedNumVariables != actualCols)
                    {
                        Console.WriteLine($"WARNING: Matrix has {actualCols} columns but expected {data.ExpectedNumVariables} variables.");
                        Console.WriteLine("This could indicate an issue with the face-vertex pair counting or matrix building.");
                    }
                    else
                    {
                        Console.WriteLine($"Matrix dimensions verified: {actualCols} columns as expected.");
                    }
                }

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
                // Detect section headers,
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

                        if (!int.TryParse(vertexParts[0], out int vId))
                        {
                            Console.WriteLine($"Line {lineNo}: Invalid VertexID '{vertexParts[0]}'. Skipping vertex.");
                            continue;
                        }

                        ContactPoint vertex = new ContactPoint
                        {
                            Id = vId,
                            X = double.Parse(vertexParts[1]),
                            Y = double.Parse(vertexParts[2])
                        };

                        if (!geometry.Vertices.TryAdd(vId, vertex))
                        {
                            Console.WriteLine($"Line {lineNo}: Duplicate vertex ID {vId} found. Skipping.");
                        }
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
                           
                            Thickness = double.Parse(faceParts[5]),
                            CohesionValue = double.Parse(faceParts[6]),
                            MuOverride = double.Parse(faceParts[7])
                        };

                        // Check that referenced vertices exist
                        bool allVerticesExist = true;
                        foreach (int vertexId in face.VertexIds)
                        {
                            if (!geometry.Vertices.ContainsKey(vertexId))
                            {
                                Console.WriteLine($"Line {lineNo}: Error: Face {face.Id} references non-existent vertex {vertexId}");
                                allVerticesExist = false;
                            }
                        }

                        // Check that referenced blocks exist - with support block handling
                        // For Block J (recognize ≤ 0 as valid support blocks)
                        if (face.BlockJ > 0 && !geometry.Blocks.ContainsKey(face.BlockJ))
                        {
                            Console.WriteLine($"Line {lineNo}: Error: Face {face.Id} references non-existent block J: {face.BlockJ}");
                        }

                        // For Block K (recognize ≤ 0 as valid support blocks)
                        if (face.BlockK > 0 && !geometry.Blocks.ContainsKey(face.BlockK))
                        {
                            Console.WriteLine($"Line {lineNo}: Error: Face {face.Id} references non-existent block K: {face.BlockK}");
                        }

                        // Only continue if all vertices exist
                        if (!allVerticesExist)
                        {
                            Console.WriteLine($"Line {lineNo}: Skipping face {face.Id} due to missing vertices");
                            continue;
                        }
    
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
            int realBlocksCount = geometry.Blocks.Values.Count(b => b.Id > 0);
            int expectedRows = realBlocksCount * 3;
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

            // Cross-check expected #variables = #face-vertex pairs * 2
            int totalFaceVertexPairs = geometry.Faces.Values.Sum(f => f.VertexIds.Count);
            int expectedNumVariables = totalFaceVertexPairs * 2;
            data.ExpectedNumVariables = expectedNumVariables;
            Console.WriteLine($"Expected #variables = {expectedNumVariables} (from {totalFaceVertexPairs} face-vertex pairs)");

            // Optional cross-check: does #faces match half the #vertices? (just a heuristic)
            int expectedNumFacesFromVertices = geometry.Vertices.Count / 2;
            if (geometry.Faces.Count != expectedNumFacesFromVertices)
            {
                Console.WriteLine($"Warning: #faces ({geometry.Faces.Count}) != #vertices/2 ({expectedNumFacesFromVertices}). Check data.");
            }
        }
        private static double DistanceBetween(ContactPoint a, ContactPoint b)
        {
            return Math.Sqrt(Math.Pow(b.X - a.X, 2) + Math.Pow(b.Y - a.Y, 2));
        }

 
        public static void ComputeFaceNormalsFromGeometry(GeometryModel geometry)
        {
            foreach (var face in geometry.Faces.Values)
            {
                if (face.VertexIds.Count != 2)
                {
                    // Skip or warn if not exactly 2 vertices
                    Console.WriteLine($"Face {face.Id} has {face.VertexIds.Count} vertices, expected 2.");
                    continue;
                }
                //  Make sure both blocks exist
                if (!geometry.Blocks.ContainsKey(face.BlockJ) ||
                    !geometry.Blocks.ContainsKey(face.BlockK))
                {
                    Console.WriteLine(
                        $"Face {face.Id}: missing block J={face.BlockJ} or K={face.BlockK} in geometry!");
                    continue;
                }
                // Step 1: pick which block is J or K by centroid comparison
                //         i.e., geometry-based rule: smaller X => J, tie => smaller Y => J

                // We'll store the original IDs in local variables
                int jId = face.BlockJ;
                int kId = face.BlockK;

                // Retrieve the actual Block objects
                var jRef = geometry.Blocks[jId];
                var kRef = geometry.Blocks[kId];

                // Decide if we need to swap them
                bool needSwap = false;
                const double eps = 1e-12;

                // Compare centroidX
                if (Math.Abs(jRef.CentroidX - kRef.CentroidX) < eps)
                {
                    // X are effectively the same, compare Y
                    if (jRef.CentroidY > kRef.CentroidY)
                        needSwap = true;
                }
                else
                {
                    // smaller X => J
                    if (jRef.CentroidX > kRef.CentroidX)
                        needSwap = true;
                }
                // If we decided to swap, swap face.BlockJ and face.BlockK
                if (needSwap)
                {
                    face.BlockJ = kId;
                    face.BlockK = jId;
                }

                // Now re-fetch references after potential swap
                jRef = geometry.Blocks[face.BlockJ];
                kRef = geometry.Blocks[face.BlockK];

                // Step 2: retrieve the two vertices that define the face
                int v1Id = face.VertexIds[0];
                int v2Id = face.VertexIds[1];

                if (!geometry.Vertices.ContainsKey(v1Id) ||
                    !geometry.Vertices.ContainsKey(v2Id))
                {
                    Console.WriteLine($"Face {face.Id}: missing one or both vertices in geometry!");
                    continue;
                }

                var v1 = geometry.Vertices[v1Id];
                var v2 = geometry.Vertices[v2Id];

                // Step 3: compute the tangent from v1 -> v2
                double dx = v2.X - v1.X;
                double dy = v2.Y - v1.Y;
                double length = Math.Sqrt(dx * dx + dy * dy);

                if (length < 1e-12)
                {
                    Console.WriteLine($"Face {face.Id}: zero-length edge. Skipping normal/tangent assignment.");
                    continue;
                }

                face.Length = length; // store face length if needed
                double tx = dx / length;
                double ty = dy / length;

                // Step 4: candidate normal = +90° rotation of tangent
                double nx = +ty;
                double ny = -tx;

                // Step 5: flip normal if dot < 0 => ensure normal points from J->K
                double dCx = (kRef.CentroidX - jRef.CentroidX);
                double dCy = (kRef.CentroidY - jRef.CentroidY);

                double dot = nx * dCx + ny * dCy;
                if (dot < 0)
                {
                    nx = -nx;
                    ny = -ny;
                }

                // Step 6: assign final results
                face.Tangent = new double[] { tx, ty };
                face.Normal = new double[] { nx, ny };
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
        // helper method to check that order of vertices is correct
        static void ValidateGeometryModel(GeometryModel geometry)
        {
            Console.WriteLine("\n=== GEOMETRY VALIDATION ===\n");

            // 1. Check all normal and tangent vectors
            Console.WriteLine("--- Vector Validation ---");
            bool allVectorsValid = true;

            foreach (var face in geometry.Faces.Values)
            {
                bool faceValid = true;
                Console.Write($"Face {face.Id}: ");

                // Check if vectors exist
                if (face.Normal == null || face.Normal.Length != 2 ||
                    face.Tangent == null || face.Tangent.Length != 2)
                {
                    Console.WriteLine("ERROR: Missing vectors");
                    allVectorsValid = false;
                    continue;
                }

                // Check normal vector length
                double normalLength = Math.Sqrt(face.Normal[0] * face.Normal[0] + face.Normal[1] * face.Normal[1]);
                bool normalLengthOk = Math.Abs(normalLength - 1.0) < 1e-6;
                if (!normalLengthOk)
                {
                    Console.Write($"Normal not unit length ({normalLength:F6}) ");
                    faceValid = false;
                }

                // Check tangent vector length
                double tangentLength = Math.Sqrt(face.Tangent[0] * face.Tangent[0] + face.Tangent[1] * face.Tangent[1]);
                bool tangentLengthOk = Math.Abs(tangentLength - 1.0) < 1e-6;
                if (!tangentLengthOk)
                {
                    Console.Write($"Tangent not unit length ({tangentLength:F6}) ");
                    faceValid = false;
                }

                // Check orthogonality
                double dotProduct = face.Normal[0] * face.Tangent[0] + face.Normal[1] * face.Tangent[1];
                bool orthogonalOk = Math.Abs(dotProduct) < 1e-6;
                if (!orthogonalOk)
                {
                    Console.Write($"Not orthogonal (dot={dotProduct:F6}) ");
                    faceValid = false;
                }

                if (faceValid)
                    Console.WriteLine("OK");
                else
                    Console.WriteLine();

                allVectorsValid &= faceValid;
            }

            // 2. Check vertex ordering consistency
            Console.WriteLine("\n--- Vertex Ordering Validation ---");
            bool allOrderingsConsistent = true;
            double? expectedCrossSign = null;

            foreach (var face in geometry.Faces.Values)
            {
                if (face.VertexIds.Count != 2)
                {
                    Console.WriteLine($"Face {face.Id}: SKIP (not exactly 2 vertices)");
                    continue;
                }

                int v1Id = face.VertexIds[0];
                int v2Id = face.VertexIds[1];

                if (!geometry.Vertices.ContainsKey(v1Id) || !geometry.Vertices.ContainsKey(v2Id))
                {
                    Console.WriteLine($"Face {face.Id}: ERROR (missing vertices)");
                    continue;
                }

                var v1 = geometry.Vertices[v1Id];
                var v2 = geometry.Vertices[v2Id];

                // Vector from v1 to v2
                double dx = v2.X - v1.X;
                double dy = v2.Y - v1.Y;

                // Cross product with normal (should be consistent sign for all faces)
                double cross = dx * face.Normal[1] - dy * face.Normal[0];

                // If this is the first face, set the expected sign
                if (!expectedCrossSign.HasValue && Math.Abs(cross) > 1e-6)
                {
                    expectedCrossSign = Math.Sign(cross);
                }

                string status = "OK";
                if (expectedCrossSign.HasValue && Math.Sign(cross) != expectedCrossSign.Value)
                {
                    status = "INCONSISTENT ORDERING";
                    allOrderingsConsistent = false;
                }

                Console.WriteLine($"Face {face.Id}: V{v1Id}→V{v2Id}, Cross={cross:F6}, {status}");
            }

            // Summary
            Console.WriteLine("\n=== VALIDATION SUMMARY ===");
            Console.WriteLine($"All vectors valid: {(allVectorsValid ? "YES" : "NO")}");
            Console.WriteLine($"Consistent vertex ordering: {(allOrderingsConsistent ? "YES" : "NO")}");
            Console.WriteLine("===============================\n");
        }


        
        
        // be informed if some input data is modified by the solver 
        static void LogDataModification(int lineNo, string message)
        {
            string formattedMsg = $"[DATA MODIFIED] Line {lineNo}: {message}";
            Console.ForegroundColor = ConsoleColor.Yellow;  // Highlight modifications
            Console.WriteLine(formattedMsg);
            Console.ResetColor();  // Reset to default color
        }


        // --- DumpFaceData Method (No changes needed based on analysis) ---
        static void DumpFaceData(GeometryModel g, string tag)
        {
            Console.WriteLine($"\n--- {tag} Face Data ---");
            Console.WriteLine("ID | BlockJ | BlockK | Thickness | Cohesion | Friction  | Vertex Order Check");
            Console.WriteLine("-----------------------------------------------------------------------------------------");

            foreach (var f in g.Faces.Values.OrderBy(f => f.Id))
            {
                // Format vectors
                string nxStr = f.Normal != null && f.Normal.Length == 2 ? f.Normal[0].ToString("F6") : "N/A";
                string nyStr = f.Normal != null && f.Normal.Length == 2 ? f.Normal[1].ToString("F6") : "N/A";
                string txStr = f.Tangent != null && f.Tangent.Length == 2 ? f.Tangent[0].ToString("F6") : "N/A";
                string tyStr = f.Tangent != null && f.Tangent.Length == 2 ? f.Tangent[1].ToString("F6") : "N/A";

                // Get friction coefficient
                string muStr = f.MuOverride.HasValue ? f.MuOverride.Value.ToString("F3") : "DEFAULT";

                // Check vertex ordering if we have exactly 2 vertices
                string vertexCheck = "";
                if (f.VertexIds.Count == 2 && g.Vertices.ContainsKey(f.VertexIds[0]) && g.Vertices.ContainsKey(f.VertexIds[1]))
                {
                    var v1 = g.Vertices[f.VertexIds[0]];
                    var v2 = g.Vertices[f.VertexIds[1]];

                    // Vector from v1 to v2
                    double dx = v2.X - v1.X;
                    double dy = v2.Y - v1.Y;

                    // Cross product with normal (to check ordering convention)
                    double cross = dx * f.Normal[1] - dy * f.Normal[0];

                    // Dot product with tangent (should be positive if vertices follow tangent)
                    double dot = dx * f.Tangent[0] + dy * f.Tangent[1];

                    vertexCheck = $"Cross={cross:F3}, Dot={dot:F3}";

                    // Add a warning flag for potentially incorrect vertex ordering
                    if (cross < 0)
                        vertexCheck += " ⚠️ REVERSED";
                }

                Console.WriteLine($"{f.Id,2} | {f.BlockJ,6} | {f.BlockK,6} |  {f.Thickness,9:F3} | " +
                                  $"{f.CohesionValue,8:F3} | {muStr,8} | {vertexCheck}");
            }

            // Check for orthogonality and unit length of vectors
            Console.WriteLine("\n--- Vector Validation ---");
            foreach (var f in g.Faces.Values.OrderBy(f => f.Id))
            {
                // Check normal length
                double normalLength = Math.Sqrt(f.Normal[0] * f.Normal[0] + f.Normal[1] * f.Normal[1]);
                bool normalOk = Math.Abs(normalLength - 1.0) < 1e-6;

                // Check tangent length
                double tangentLength = Math.Sqrt(f.Tangent[0] * f.Tangent[0] + f.Tangent[1] * f.Tangent[1]);
                bool tangentOk = Math.Abs(tangentLength - 1.0) < 1e-6;

                // Check orthogonality
                double dotProduct = f.Normal[0] * f.Tangent[0] + f.Normal[1] * f.Tangent[1];
                bool orthogonalOk = Math.Abs(dotProduct) < 1e-6;

                string statusMsg = "OK";
                if (!normalOk || !tangentOk || !orthogonalOk)
                {
                    statusMsg = "ISSUE: ";
                    if (!normalOk) statusMsg += $"Normal length={normalLength:F6} ";
                    if (!tangentOk) statusMsg += $"Tangent length={tangentLength:F6} ";
                    if (!orthogonalOk) statusMsg += $"Dot product={dotProduct:F6} ";
                }

                Console.WriteLine($"Face {f.Id}: {statusMsg}");
            }
        }










    } // end of Program

} //end namespace
