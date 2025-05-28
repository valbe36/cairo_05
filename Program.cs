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
                        matrixA[jRow, colN] += n[0];
                        matrixA[jRow, colT] += t[0];
                        // Fy
                        matrixA[jRow + 1, colN] += n[1];
                        matrixA[jRow + 1, colT] += t[1];
                        // M
                        matrixA[jRow + 2, colN] += (xRelJ * n[1] - yRelJ * n[0]);
                        matrixA[jRow + 2, colT] += (xRelJ * t[1] - yRelJ * t[0]);
                    }

                    // If K is non‐support, add the “plus” contribution to blockK
                    if (!isKSupport && blockRowMap.TryGetValue(face.BlockK, out int kRow))
                    {
                        var bK = geometry.Blocks[face.BlockK];
                        double xRelK = cp.X - bK.CentroidX;
                        double yRelK = cp.Y - bK.CentroidY;

                        matrixA[kRow, colN] += - n[0];
                        matrixA[kRow, colT] += - t[0];
                        matrixA[kRow + 1, colN] += -n[1];
                        matrixA[kRow + 1, colT] += - t[1];
                        matrixA[kRow + 2, colN] += - (xRelK * n[1] - yRelK * n[0]);
                        matrixA[kRow + 2, colT] += - (xRelK * t[1] - yRelK * t[0]);
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
                    return a.FaceId.CompareTo(b.FaceId);
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
                    AddFaceEccConstraintsAroundV1(model, geometry, data);

                        // 7) Objective: maximize lambda
                        GRBLinExpr obj = 0.0;
                        obj.AddTerm(1.0, lambda);
                        model.SetObjective(obj, GRB.MAXIMIZE);
                        model.Write("debugModel.lp");
                        DumpColumnMap(data);

                    model.Parameters.NumericFocus = 3; // Maximum precision
                    model.Parameters.FeasibilityTol = 1e-9; // Tighter feasibility tolerance
                    model.Parameters.OptimalityTol = 1e-9; // Tighter optimality tolerance
                    model.Parameters.MIPGap = 0.0005;
                    model.Parameters.TimeLimit = 100; // Limit to 60 seconds
                    //model.Parameters.Quad = 1;       // Convex quadratic relaxation
                    // model.Parameters.PreQLinearize = 1; // Linearize quadratic terms

                    // 8) Solve
                    model.Optimize();
                        SaveResultsToFile(model, @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\results_cairo.txt", data);
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
                if (data.B != null && i < data.B.Length)
                    rhs.AddTerm(data.B[i], lambda);

                model.AddConstr(lhs == rhs, $"Equil_{i}");
                }
            }


        /// For each pair j, friction constraints: -mu * fN_j <= fT_j <= mu * fN_j
        /// plus no tension fN_j >= 0 (already in the variable bound).
        private void AddContactConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            foreach (var fKvp in geometry.Faces)
            {
                Face face = fKvp.Value;
                double mu = face.MuOverride ?? data.Mu;
                double c = face.CohesionValue;
                double area = face.Thickness * face.Depth;
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
                GRBVar fnk = model.AddVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, $"fnk_face{face.Id}");
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




       /*
        private void AddFaceEccConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            double sigmaC = data.SigmaC;

            //  debugging
            Console.WriteLine("\n--- Face Eccentricity Sign Determination ---");

            foreach (var fKvp in geometry.Faces)
            {
                Face face = fKvp.Value;
                // We assume each face in 2D has exactly 2 vertices.
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
                // Calculate moment arm direction vector from vertex1 to vertex2
                double dx = vertex2.X - vertex1.X;
                double dy = vertex2.Y - vertex1.Y;
                
                // Get coordinates of the vertices
                var v1 = _geometry.Vertices[vertexId1];
                var v2 = _geometry.Vertices[vertexId2];

                // Calculate the face center
                double faceCenterX = (v1.X + v2.X) / 2.0;
                double faceCenterY = (v1.Y + v2.Y) / 2.0;

                // Distance from vertices to face center
                double dx1 = v1.X - faceCenterX;
                double dy1 = v1.Y - faceCenterY;
                double dx2 = v2.X - faceCenterX;
                double dy2 = v2.Y - faceCenterY;

                // Signed moment arms (positive if vertex is on "positive" side of face)
                // Dot product of vertex-to-center vector with tangent vector
                double sign1 = Math.Sign(dx1 * face.Tangent[0] + dy1 * face.Tangent[1]);
                double sign2 = Math.Sign(dx2 * face.Tangent[0] + dy2 * face.Tangent[1]);

                // Project the distance to the face normal direction
                double momentArm1 = sign1 * Lk / 2.0;  // Simplify to just use half the face length
                double momentArm2 = sign2 * Lk / 2.0;  // This should always be -momentArm1

                Console.WriteLine($"Face {face.Id}: momentArm1={momentArm1:F6}, momentArm2={momentArm2:F6}");

                // Standard moment equation: fnk * eK = fn1 * momentArm1 + fn2 * momentArm2
                GRBQuadExpr mq = new GRBQuadExpr();
                mq.AddTerm(1.0, fnk, eK);
                mq.AddTerm(momentArm1, fn1);
                mq.AddTerm(momentArm2, fn2);
                // mq.AddTerm( Lk / 2.0, fn1);
                //  mq.AddTerm(- Lk / 2.0, fn2);
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
      
        */
        private void DumpColumnMap(ProblemData data)
        {
            string path = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\mapping_cairo.txt";
            using var w = new StreamWriter(path);
            // -----------------------------------------------------------
            // 1) Column mapping
            // -----------------------------------------------------------
            w.WriteLine("=== COLUMN MAPPING ===");
            w.WriteLine("Column | Face-Id, Vertex-Id | Variable");
            w.WriteLine("-------------------------");
            for (int j = 0; j < faceVertexPairs.Count; j++)
            {
                var p = faceVertexPairs[j];
                int faceId = p.FaceId;
                int vertexId = p.VertexId;

                // Two columns per pair: Normal, Tangential
                int colN = 2 * j;
                int colT = colN + 1;

                w.WriteLine($"{colN,5} | F{faceId}, V{vertexId} | Normal");
                w.WriteLine($"{colT,5} | F{faceId}, V{vertexId} | Tangential");
                w.WriteLine("-------------------------");
            }

            // -----------------------------------------------------------
            // 2) Print the full A-matrix
            // -----------------------------------------------------------
            w.WriteLine("\n=== FULL A-MATRIX ===");
            if (data.MatrixA != null)
            {
                int nR = data.NumRows;
                int nC = data.NumCols;
                w.WriteLine($"Matrix dimensions: {nR} x {nC}");
                for (int i = 0; i < nR; i++)
                {
                    w.Write($"Row {i,2}: ");
                    for (int j = 0; j < nC; j++)
                    {
                        w.Write($"{data.MatrixA[i, j],10:F4} ");
                    }
                    w.WriteLine();
                }
            }
            else
            {
                w.WriteLine("No MatrixA defined.");
            }

            // -----------------------------------------------------------
            // 3) Print the full B vector
            // -----------------------------------------------------------
            w.WriteLine("\n=== VECTOR B ===");
            if (data.B != null && data.B.Length > 0)
            {
                for (int i = 0; i < data.B.Length; i++)
                {
                    w.WriteLine($"B[{i}] = {data.B[i]:F6}");
                }
            }
            else
            {
                w.WriteLine("No B vector defined or empty.");
            }

            // -----------------------------------------------------------
            // 4) Print the full G vector
            // -----------------------------------------------------------
            w.WriteLine("\n=== VECTOR G ===");
            if (data.G != null && data.G.Length > 0)
            {
                for (int i = 0; i < data.G.Length; i++)
                {
                    w.WriteLine($"G[{i}] = {data.G[i]:F6}");
                }
            }
            else
            {
                w.WriteLine("No G vector defined or empty.");
            }

            // -----------------------------------------------------------
            // 5) Print block centroid + corner vertices
            // -----------------------------------------------------------
            w.WriteLine("\n=== BLOCK CORNERS & CENTROIDS ===");

            // Build a map: blockId -> set of vertexIds
            var blockCorners = new Dictionary<int, HashSet<int>>();
            // Initialize sets so every block has an entry, even if empty
            foreach (var blk in _geometry.Blocks.Values)
            {
                blockCorners[blk.Id] = new HashSet<int>();
            }

            // Populate sets: scan each face
            foreach (var face in _geometry.Faces.Values)
            {
                // For blockJ
                if (_geometry.Blocks.ContainsKey(face.BlockJ))
                {
                    // Add all face vertices to that block's set
                    foreach (int vid in face.VertexIds)
                    {
                        blockCorners[face.BlockJ].Add(vid);
                    }
                }
                // For blockK
                if (_geometry.Blocks.ContainsKey(face.BlockK))
                {
                    foreach (int vid in face.VertexIds)
                    {
                        blockCorners[face.BlockK].Add(vid);
                    }
                }
            }

            // Now print the corners for each block
            foreach (var blk in _geometry.Blocks.Values.OrderBy(b => b.Id))
            {
                double cx = blk.CentroidX;
                double cy = blk.CentroidY;
                w.WriteLine($"\nBlock ID={blk.Id}, Centroid=({cx:F3}, {cy:F3})");

                // Retrieve all vertex IDs for this block
                var cornerIds = blockCorners[blk.Id].OrderBy(id => id).ToList();
                if (cornerIds.Count == 0)
                {
                    w.WriteLine("  No corners found (maybe it's a support with no faces?).");
                    continue;
                }

                w.WriteLine("  Corner Vertices:");
                foreach (int vid in cornerIds)
                {
                    if (_geometry.Vertices.TryGetValue(vid, out var cpt))
                    {
                        w.WriteLine($"   - vId={vid}, coords=({cpt.X:F3}, {cpt.Y:F3})");
                    }
                    else
                    {
                        w.WriteLine($"   - vId={vid}, (Not found in geometry.Vertices)");
                    }
                }
            }

            w.WriteLine("\n=== End of matrix structure analysis ===");
            w.Flush();

            Console.WriteLine($"Matrix, B/G, block corners, and centroids dumped to {path}");


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


                // 5) Shear limits at the FACE level (not vertex level)
                writer.WriteLine("\n--- Face-level Shear Check ---");
                double tolerance = 1e-6;

                foreach (var kvp in faceTotalForces)
                {
                    int faceId = kvp.Key;
                    double totalNormal = kvp.Value.fnSum;
                    double totalShear = kvp.Value.ftSum;

                    // Get face data
                    Face face = _geometry.Faces[faceId];
                    double mu = face.MuOverride ?? data.Mu;
                    double cohesion = face.CohesionValue;
                    double area = face.Depth * face.Thickness;

                    // Use full cohesion (not half per vertex)
                    double cohesionTerm = cohesion * area;

                    // Calculate shear limit
                    double shearLimit = totalNormal * mu + cohesionTerm;
                    double excessShear = Math.Abs(totalShear) - shearLimit;

                    // Determine status with tolerance
                    string shearStatus;
                    if (excessShear >= -tolerance && excessShear <= tolerance)
                        shearStatus = "FAILURE";
                    else
                        shearStatus = "OK";

                    writer.WriteLine($"Face {faceId}: |total_ft|={Math.Abs(totalShear):F3} vs. limit={shearLimit:F3}, Excess={excessShear:F3},Status={shearStatus}");
                    writer.WriteLine($"   Formula: |f_t| ≤ f_n * mu + c*A\"), Details: |{totalShear:F3}| ≤ {totalNormal:F3} * {mu:F3} + {cohesion:F3} * {area:F3}");

                }

                // 6) Face Eccentricity Check
                writer.WriteLine("\n--- Face Eccentricity Check ---");
                foreach (var kvp in faceEccVars)
                {
                    int faceId = kvp.Key;
                    Face face = _geometry.Faces[faceId];
                    GRBVar eccVar = kvp.Value;
                    double eccVal = eccVar.X;

                    // Use total normal force from our already calculated dictionary
                    double totalNormal = faceTotalForces.ContainsKey(faceId) ?
                                        faceTotalForces[faceId].fnSum : 0.0;

                    double depth = face.Depth;
                    double thickness = face.Thickness;
                    double sigmaC = data.SigmaC;

                    // Calculate eccentricity limit
                    double strengthLimit = depth / 2.0 - totalNormal / (2.0 * sigmaC * thickness);
                    double excessEcc = Math.Abs(eccVal) - Math.Abs(strengthLimit);

                    // Determine status with tolerance
                    string eccStatus;
                    if (excessEcc >= -tolerance && excessEcc <= tolerance)
                        eccStatus = "FAILURE";
                    else
                        eccStatus = "OK";

                    writer.WriteLine($"Face {faceId}: eccentricity = {eccVal:F3} | " +
                                    $"Total Normal = {totalNormal:F3}, Limit = ±{Math.Abs(strengthLimit):F3}, " +
                                    $"Excess =  {eccStatus}");
                    writer.WriteLine($"  Formula: |e| ≤ d/2 - N/(2*σc*t)");
                    writer.WriteLine($"  Details: |{eccVal:F3}| ≤ {depth:F3}/2 - {totalNormal:F3}/(2*{sigmaC:F3}*{thickness:F3})");
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
                LoadAllData(@"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\/data_pseudoparallel_friction_0e4.txt"   //data_pseudoparallel_friction_0e4
                , geometry, data);
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
                            Depth = double.Parse(faceParts[9]),
                            Thickness = double.Parse(faceParts[10]),
                            CohesionValue = double.Parse(faceParts[11]),
                            MuOverride = double.Parse(faceParts[12])
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


                        // Check if face.Depth matches actual distance between vertices
                        double actualDistance = Math.Sqrt(
                            Math.Pow(geometry.Vertices[face.VertexIds[1]].X - geometry.Vertices[face.VertexIds[0]].X, 2) +
                            Math.Pow(geometry.Vertices[face.VertexIds[1]].Y - geometry.Vertices[face.VertexIds[0]].Y, 2)
                        );
                        if (Math.Abs(actualDistance - face.Depth) > 1e-6)
                        {
                            Console.WriteLine($"Line {lineNo}: Warning: Face {face.Id} depth ({face.Depth}) differs from actual vertex distance ({actualDistance})");
                            Console.WriteLine($"Line {lineNo}: Auto-correcting face depth for Face {face.Id}");
                            face.Depth = actualDistance;
                        }
                        if (Math.Abs(actualDistance - face.Depth) > 1e-6)
                        {
                            double originalDepth = face.Depth;
                            face.Depth = actualDistance;
                            LogDataModification(lineNo, $"Face {face.Id} depth corrected: {originalDepth:F6} → {actualDistance:F6}");
                        }

                        // Parse vectors
                        double nx = double.Parse(faceParts[5]);
                        double ny = double.Parse(faceParts[6]);
                        double tx = double.Parse(faceParts[7]);
                        double ty = double.Parse(faceParts[8]);

                        // Check normal-tangent orthogonality
                        double dotProduct = nx * tx + ny * ty;
                        if (Math.Abs(dotProduct) > 1e-6)
                        {
                            Console.WriteLine($"Line {lineNo}: Warning: Face {face.Id} normal/tangent not orthogonal (dot product = {dotProduct})");
                            Console.WriteLine($"Line {lineNo}: Auto-correcting tangent vector for Face {face.Id}");
                            // Make tangent orthogonal to normal (rotate 90 degrees)
                            tx = -ny;
                            ty = nx;
                            // Re-normalize tangent
                            double tangentLength = Math.Sqrt(tx * tx + ty * ty);
                            tx /= tangentLength;
                            ty /= tangentLength;
                        }
                        // Ensure normal vector is normalized
                        double normalLength = Math.Sqrt(nx * nx + ny * ny);
                        if (Math.Abs(normalLength - 1.0) > 1e-6)
                        {
                            Console.WriteLine($"Line {lineNo}: Warning! Face {face.Id} normal not unit length ({normalLength:F6})");
                            nx /= normalLength;
                            ny /= normalLength;
                            Console.WriteLine($"Line {lineNo}: Normalized normal vector for Face {face.Id}: ({nx:F6}, {ny:F6})");
                        }

                        // --- VALIDATION CHECKS ADDED HERE ---


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
            // Cross-check: Count the expected number of variables
            int totalFaceVertexPairs = 0;
            foreach (var face in geometry.Faces.Values)
            {
                totalFaceVertexPairs += face.VertexIds.Count;
            }
            int expectedNumVariables = totalFaceVertexPairs * 2; // 2 variables (fN, fT) per vertex-face pair

            // Store the expected count for later validation after matrix is built
            data.ExpectedNumVariables = expectedNumVariables;

            Console.WriteLine($"Expected number of variables: {expectedNumVariables} (from {totalFaceVertexPairs} face-vertex pairs)");

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
            Console.WriteLine("ID | BlockJ | BlockK | Length | Thickness | Cohesion | Friction | Normal | Tangent | Vertex Order Check");
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

                Console.WriteLine($"{f.Id,2} | {f.BlockJ,6} | {f.BlockK,6} | {f.Depth,6:F3} | {f.Thickness,9:F3} | " +
                                  $"{f.CohesionValue,8:F3} | {muStr,8} | ({nxStr},{nyStr}) | ({txStr},{tyStr}) | {vertexCheck}");
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
