using Gurobi;
using System.Globalization;
using static System.Runtime.InteropServices.JavaScript.JSType;

/*
Mathematical programming model for masonry stability check.
Primary variables: fx[i], fy[i] (global forces at vertex i).
Constraints: Block equilibrium, Vertex No-Tension, Vertex Friction, Face Strength/Eccentricity.
Goal: Maximize load factor lambda.

Input File:
[Vertices]
Header:[Blocks]
Format (per line):  BlockID, CentroidX, CentroidY, AppliedFx, AppliedFy, x_coord (1), y_coord (1)
Header:[Vertices]
Format (per line):vertexID, x_coord, y_coord
Header: [Faces]
Format (per line): faceID, ID_block1, ID_block2, ID_pt1, ID_pt2, Thickness, Cohesion, friction coeff. 
Header: [VectorG]
   
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

            int numRows = numBlocks * 6;
            int numCols = columnMap.Count * 3; // 3 variables (N,,t1, t2) per (face,vertex)
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
                double[] t1 = face.Tangent1;   
                double[] t2 = face.Tangent2;

                foreach (int vId in face.VertexIds)
                {
                    // locate the column for this face‐vertex
                    int colPairIndex = columnMap[(face.Id, vId)];
                    int colN = colPairIndex * 3;      
                    int colT1 = colN + 1;             
                    int colT2 = colN + 2;

                    ContactPoint cp = geometry.Vertices[vId];

                    // If J is non‐support, add the “minus” contribution to blockJ
                    if (!isJSupport && blockRowMap.TryGetValue(face.BlockJ, out int jRow))
                    {
                        // local coords of cp relative to blockJ centroid
                        var bJ = geometry.Blocks[face.BlockJ];
                        double xRelJ = cp.X - bJ.CentroidX;
                        double yRelJ = cp.Y - bJ.CentroidY;
                        double zRelJ = cp.Z - bJ.CentroidZ;
                        // Fx equation
                        matrixA[jRow, colN] -= n[0];      // Normal force X component
                        matrixA[jRow, colT1] -= t1[0];    // Tangent1 force X component
                        matrixA[jRow, colT2] -= t2[0];    

                        // Fy equation  
                        matrixA[jRow + 1, colN] -= n[1];  // Normal force Y component
                        matrixA[jRow + 1, colT1] -= t1[1]; // 
                        matrixA[jRow + 1, colT2] -= t2[1];  

                        // Fz equation (3d only)
                        matrixA[jRow + 2, colN] -= n[2];  // Normal force Z component
                        matrixA[jRow + 2, colT1] -= t1[2]; //  
                        matrixA[jRow + 2, colT2] -= t2[2]; //  

                        // MOMENT EQUILIBRIUM EQUATIONS  
                        // Mx = y*Fz - z*Fy ( 3d only
                        matrixA[jRow + 3, colN] -= (yRelJ * n[2] - zRelJ * n[1]);
                        matrixA[jRow + 3, colT1] -= (yRelJ * t1[2] - zRelJ * t1[1]);
                        matrixA[jRow + 3, colT2] -= (yRelJ * t2[2] - zRelJ * t2[1]);

                        // My = z*Fx - x*Fz  
                        matrixA[jRow + 4, colN] -= (zRelJ * n[0] - xRelJ * n[2]);
                        matrixA[jRow + 4, colT1] -= (zRelJ * t1[0] - xRelJ * t1[2]);
                        matrixA[jRow + 4, colT2] -= (zRelJ * t2[0] - xRelJ * t2[2]);

                        // Mz = x*Fy - y*Fx (SAME as before, but now at row+5)
                        matrixA[jRow + 5, colN] -= (xRelJ * n[1] - yRelJ * n[0]);
                        matrixA[jRow + 5, colT1] -= (xRelJ * t1[1] - yRelJ * t1[0]);
                        matrixA[jRow + 5, colT2] -= (xRelJ * t2[1] - yRelJ * t2[0]);
                    }

                    // If K is non‐support, add the “plus” contribution to blockK
                    if (!isKSupport && blockRowMap.TryGetValue(face.BlockK, out int kRow))
                    {
                        var bK = geometry.Blocks[face.BlockK];
                        double xRelK = cp.X - bK.CentroidX;
                        double yRelK = cp.Y - bK.CentroidY;
                        double zRelK = cp.Z - bK.CentroidZ;

                        // Fx equation
                        matrixA[kRow, colN] += n[0];
                        matrixA[kRow, colT1] += t1[0];
                        matrixA[kRow, colT2] += t2[0];

                        // Fy equation
                        matrixA[kRow + 1, colN] += n[1];
                        matrixA[kRow + 1, colT1] += t1[1];
                        matrixA[kRow + 1, colT2] += t2[1];

                        // Fz equation (NEW)
                        matrixA[kRow + 2, colN] += n[2];
                        matrixA[kRow + 2, colT1] += t1[2];
                        matrixA[kRow + 2, colT2] += t2[2];

                        // MOMENT EQUILIBRIUM EQUATIONS (same structure as J, but with + signs)
                        // Mx = y*Fz - z*Fy
                        matrixA[kRow + 3, colN] += (yRelK * n[2] - zRelK * n[1]);
                        matrixA[kRow + 3, colT1] += (yRelK * t1[2] - zRelK * t1[1]);
                        matrixA[kRow + 3, colT2] += (yRelK * t2[2] - zRelK * t2[1]);

                        // My = z*Fx - x*Fz
                        matrixA[kRow + 4, colN] += (zRelK * n[0] - xRelK * n[2]);
                        matrixA[kRow + 4, colT1] += (zRelK * t1[0] - xRelK * t1[2]);
                        matrixA[kRow + 4, colT2] += (zRelK * t2[0] - xRelK * t2[2]);

                        // Mz = x*Fy - y*Fx
                        matrixA[kRow + 5, colN] += (xRelK * n[1] - yRelK * n[0]);
                        matrixA[kRow + 5, colT1] += (xRelK * t1[1] - yRelK * t1[0]);
                        matrixA[kRow + 5, colT2] += (xRelK * t2[1] - yRelK * t2[0]);
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
        public double Z { get; set; }  // ADD THIS LINE ONLY
    }
    /// Represents a 3D contact face, assumed n vertexes, but adaptable  <summary>
    /// </summary>
    public class Face
    {
        public int Id { get; set; }
        public double Thickness { get; set; }
        public double Length { get; set; }
        public double Area { get; set; }
        public double CohesionValue { get; set; }
        public double? MuOverride { get; set; }

        public List<int> VertexIds { get; set; } = new List<int>();
        public double[] Normal { get; set; } = new double[3];    
        public double[] Tangent { get; set; } = new double[3];   
        public double[] Tangent1 { get; set; } = new double[3];  
        public double[] Tangent2 { get; set; } = new double[3];
        public double[] Centroid { get; set; } = new double[3];
        public int BlockJ { get; set; }  // Block on the "from" side of the normal
        public int BlockK { get; set; }  // Block on the "to" side of the normal
    }

    public class Block
    {
        public int Id { get; set; }
        public double CentroidX { get; set; }  // Centroid X
        public double CentroidY { get; set; }  // Centroid Y
        public double CentroidZ { get; set; }  // ADD THIS LINE ONLY

        // Applied loads
        public double AppliedFx { get; set; } = 0.0;
        public double AppliedFy { get; set; } = 0.0;
        public double AppliedFz { get; set; } = 0.0;  // ADD THIS LINE ONLY

        public double LoadApplicationX { get; set; } = 0.0;
        public double LoadApplicationY { get; set; } = 0.0;
        public double LoadApplicationZ { get; set; } = 0.0;  // ADD THIS LINE ONLY
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



                    model.VerifyEquilibrium(data, faceVertexPairs, fAll, lambda, geometry);

                    // 5) Add friction & no-tension constraints
                    AddContactConstraints(model, geometry, data);

                    // 6)   Add face eccentricity constraints 
            
                    AddBiaxialBending(model, geometry, data);
                    // AddFaceMidpointBendingConstraints(model, geometry, data);

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



        // DEBUG 1 : to be written
        private void VerifyPostSolution(GRBModel model, ProblemData data, double friction)
        {
        }

        /// Create 3 local variables (f_n, f_t1, f_t2) for each face-vertex pair,plus the load factor lambda.
        /// <param name="model"></param>
        /// <param name="data"></param>

        private void CreateVariables(GRBModel model, ProblemData data)
        {
            int m = faceVertexPairs.Count; // number of pairs
            fAll = new GRBVar[3 * m];       

            for (int j = 0; j < m; j++)
            {
                // f_n >= 0
                fAll[3 * j] = model.AddVar(        
                    0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS,
                    $"fN_{j}"
                );

                // f_t1 unbounded (renamed from f_t)
                fAll[3 * j + 1] = model.AddVar(     
                    -GRB.INFINITY, GRB.INFINITY, 0.0, GRB.CONTINUOUS,
                    $"fT1_{j}"                      
                );

                // f_t2 unbounded (NEW VARIABLE)
                fAll[3 * j + 2] = model.AddVar(     
                    -GRB.INFINITY, GRB.INFINITY, 0.0, GRB.CONTINUOUS,
                    $"fT2_{j}"
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

                // Subtract B[i] * lambda (moved to left side)
                if (data.B != null && i < data.B.Length)
                    lhs.AddTerm(data.B[i], lambda);

                // Equilibrium: A*fAll - G - B*lambda = 0
                model.AddConstr(lhs == 0.0, $"Equil_{i}");
                }
            }


        // vertex based friction limit
        
        private void AddContactConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            foreach (var faceKvp in geometry.Faces)
            {
                Face face = faceKvp.Value;
                double mu = face.MuOverride ?? data.Mu;
                double c = face.CohesionValue;

                // Distribute total face cohesion equally among ALL vertices
                // (In reality, cohesion acts only on compressed zone, but we don't know that yet)
                // This is conservative because it distributes cohesion among all vertices, 
                double totalCohesionForce = c * face.Area;
                double cPerVertex = totalCohesionForce / face.VertexIds.Count;
                foreach (int vId in face.VertexIds)
                {
                    int idx = GetPairIndex(face.Id, vId);
                    if (idx < 0) continue; // skip if not found

                    int baseIndex = 3 * idx;  // 3 variables per vertex: fn, ft1, ft2

                    GRBVar fN = fAll[baseIndex];     // normal force
                    GRBVar fT1 = fAll[baseIndex + 1]; // first tangential force
                    GRBVar fT2 = fAll[baseIndex + 2]; // second tangential force

                    // 1) NO-TENSION CONSTRAINT 
                    model.AddConstr(fN >= 0.0, $"NoTension3D_face{face.Id}_v{vId}");

                    // 2) Vertex friction (with partial cohesion)
                    if (mu <= 1e-12 && c <= 1e-12)
                    {
                        // FRICTIONLESS CASE
                        model.AddConstr(fT1 == 0.0, $"NoShear1_face{face.Id}_v{vId}");
                        model.AddConstr(fT2 == 0.0, $"NoShear2_face{face.Id}_v{vId}");
                    }
                    else
                    {
                        // CONSERVATIVE OCTAGONAL PYRAMID APPROXIMATION
                        // Each vertex gets its own friction pyramid with distributed cohesion
                        // This is safer than face-level constraints and conservative because:
                        // 1. Uses inscribed octagon (3-4% under-estimation of true capacity)
                        // 2. Distributes cohesion among all vertices (not just compressed ones)
                        // 3. Per-vertex constraints are more restrictive than face-level

                        GRBLinExpr frictionLimit = mu * fN + cPerVertex;

                        // OCTAGONAL CONSTRAINTS (8 linear inequalities)

                        // 1. Axis-aligned constraints (4 constraints)
                        model.AddConstr(fT1 <= frictionLimit, $"FricOct1_face{face.Id}_v{vId}");   //  ft1 ≤ μfn + c
                        model.AddConstr(-fT1 <= frictionLimit, $"FricOct2_face{face.Id}_v{vId}");  // -ft1 ≤ μfn + c
                        model.AddConstr(fT2 <= frictionLimit, $"FricOct3_face{face.Id}_v{vId}");   //  ft2 ≤ μfn + c
                        model.AddConstr(-fT2 <= frictionLimit, $"FricOct4_face{face.Id}_v{vId}");  // -ft2 ≤ μfn + c

                        // 2. Diagonal constraints (4 constraints)
                        // These cut off the corners of the square to better approximate the circle
                        double sqrt2 = Math.Sqrt(2.0);
                        GRBLinExpr diagonalLimit = sqrt2 * frictionLimit;

                        model.AddConstr(fT1 + fT2 <= diagonalLimit, $"FricOct5_face{face.Id}_v{vId}");   //  (ft1 + ft2) ≤ √2(μfn + c)
                        model.AddConstr(fT1 - fT2 <= diagonalLimit, $"FricOct6_face{face.Id}_v{vId}");   //  (ft1 - ft2) ≤ √2(μfn + c)
                        model.AddConstr(-fT1 + fT2 <= diagonalLimit, $"FricOct7_face{face.Id}_v{vId}");  // (-ft1 + ft2) ≤ √2(μfn + c)
                        model.AddConstr(-fT1 - fT2 <= diagonalLimit, $"FricOct8_face{face.Id}_v{vId}");  // (-ft1 - ft2) ≤ √2(μfn + c)
                    }
                }
            }
        }


        // face based friction limit
        /*
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
          */


        /// <summary>
        /// //// Method  BiaxialBending and compression around midpoint
        /// </summary>
        private void AddBiaxialBendingConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            double sigmaC = data.SigmaC;

            Console.WriteLine("\n--- 3D Biaxial Bending Constraints (Quadratic) ---");

            foreach (var fKvp in geometry.Faces)
            {
                Face face = fKvp.Value;
                if (face.VertexIds.Count < 4)
                {
                    Console.WriteLine($"Face {face.Id}: Need 4 vertices for biaxial bending, found {face.VertexIds.Count}. Skipping.");
                    continue;
                }

                // STEP 1: Compute face dimensions
                var vertices = face.VertexIds.Select(vId => geometry.Vertices[vId]).ToList();
                (double L1, double L2) = ComputeFaceDimensions(vertices, face.Tangent1, face.Tangent2);

                Console.WriteLine($"Face {face.Id}: Dimensions L1={L1:F3} (along T1), L2={L2:F3} (along T2)");

                // STEP 2: Total normal force on face
                GRBVar fnTotal = model.AddVar(0.0, sigmaC * face.Area, 0.0, GRB.CONTINUOUS, $"fnTotal_face{face.Id}");
                {
                    GRBLinExpr sumExpr = new GRBLinExpr();
                    foreach (int vId in face.VertexIds)
                    {
                        int idx = GetPairIndex(face.Id, vId);
                        if (idx >= 0)
                        {
                            sumExpr.AddTerm(1.0, fAll[3 * idx]); // Sum normal forces
                        }
                    }
                    model.AddConstr(fnTotal == sumExpr, $"Def_fnTotal_face{face.Id}");
                }

                // STEP 3: Eccentricity variables
                GRBVar e1 = model.AddVar(-L1 / 2.0, L1 / 2.0, 0.0, GRB.CONTINUOUS, $"e1_face{face.Id}");
                GRBVar e2 = model.AddVar(-L2 / 2.0, L2 / 2.0, 0.0, GRB.CONTINUOUS, $"e2_face{face.Id}");

                // STEP 4: Absolute value variables for eccentricities
                // |e1| and |e2| using linear constraints
                GRBVar e1_abs = model.AddVar(0.0, L1 / 2.0, 0.0, GRB.CONTINUOUS, $"e1_abs_face{face.Id}");
                GRBVar e2_abs = model.AddVar(0.0, L2 / 2.0, 0.0, GRB.CONTINUOUS, $"e2_abs_face{face.Id}");

                // Linear constraints to define absolute values
                model.AddConstr(e1_abs >= e1, $"AbsVal1a_face{face.Id}");
                model.AddConstr(e1_abs >= -e1, $"AbsVal1b_face{face.Id}");
                model.AddConstr(e2_abs >= e2, $"AbsVal2a_face{face.Id}");
                model.AddConstr(e2_abs >= -e2, $"AbsVal2b_face{face.Id}");

                // STEP 5: Moment equilibrium constraints
                AddBiaxialBending(model, face, geometry, fnTotal, e1, e2);

                // STEP 6: BIAXIAL BENDING CONSTRAINT
                // fnTotal ≤ σc * (L1 - 2*e1_abs) * (L2 - 2*e2_abs)
                // This is bilinear, but let Gurobi handle it

                // Effective lengths
                GRBVar L1_eff = model.AddVar(0.0, L1, 0.0, GRB.CONTINUOUS, $"L1_eff_face{face.Id}");
                GRBVar L2_eff = model.AddVar(0.0, L2, 0.0, GRB.CONTINUOUS, $"L2_eff_face{face.Id}");

                model.AddConstr(L1_eff == L1 - 2 * e1_abs, $"L1_eff_def_face{face.Id}");
                model.AddConstr(L2_eff == L2 - 2 * e2_abs, $"L2_eff_def_face{face.Id}");

                // Bilinear constraint: fnTotal ≤ σc * L1_eff * L2_eff
                // Let Gurobi handle this quadratic constraint
                GRBQuadExpr biaxialLimit = new GRBQuadExpr();
                biaxialLimit.AddTerm(sigmaC, L1_eff, L2_eff);

                model.AddQConstr(fnTotal <= biaxialLimit, $"BiaxialBending_face{face.Id}");

                Console.WriteLine($"Face {face.Id}: Added biaxial constraint - let Gurobi handle nonlinearity");

                // Store variables for debugging/analysis
                StoreBiaxialVariables(face.Id, fnTotal, e1, e2, e1_abs, e2_abs, L1_eff, L2_eff);
            }
        }

        /// <summary>
        /// Store biaxial variables for post-solution analysis
        /// </summary>
        private Dictionary<int, BiaxialVariables> biaxialVars = new Dictionary<int, BiaxialVariables>();

        private void StoreBiaxialVariables(int faceId, GRBVar fnTotal, GRBVar e1, GRBVar e2,
                                          GRBVar e1_abs, GRBVar e2_abs, GRBVar L1_eff, GRBVar L2_eff)
        {
            biaxialVars[faceId] = new BiaxialVariables
            {
                FnTotal = fnTotal,
                E1 = e1,
                E2 = e2,
                E1_Abs = e1_abs,
                E2_Abs = e2_abs,
                L1_Eff = L1_eff,
                L2_Eff = L2_eff
            };
        }

        /// <summary>
        /// Helper class to store biaxial variables for analysis
        /// </summary>
        private class BiaxialVariables
        {
            public GRBVar FnTotal { get; set; }
            public GRBVar E1 { get; set; }
            public GRBVar E2 { get; set; }
            public GRBVar E1_Abs { get; set; }
            public GRBVar E2_Abs { get; set; }
            public GRBVar L1_Eff { get; set; }
            public GRBVar L2_Eff { get; set; }
        }

        /// <summary>
        /// Compute face dimensions along tangent1 and tangent2 directions
        /// </summary>
        private (double L1, double L2) ComputeFaceDimensions(List<ContactPoint> vertices, double[] tangent1, double[] tangent2)
        {
            // Project all vertices onto tangent1 and tangent2 directions
            double minT1 = double.MaxValue, maxT1 = double.MinValue;
            double minT2 = double.MaxValue, maxT2 = double.MinValue;

            // Use face centroid as reference point
            double centX = vertices.Average(v => v.X);
            double centY = vertices.Average(v => v.Y);
            double centZ = vertices.Average(v => v.Z);

            foreach (var vertex in vertices)
            {
                // Vector from centroid to vertex
                double dx = vertex.X - centX;
                double dy = vertex.Y - centY;
                double dz = vertex.Z - centZ;

                // Project onto tangent directions
                double projT1 = dx * tangent1[0] + dy * tangent1[1] + dz * tangent1[2];
                double projT2 = dx * tangent2[0] + dy * tangent2[1] + dz * tangent2[2];

                minT1 = Math.Min(minT1, projT1);
                maxT1 = Math.Max(maxT1, projT1);
                minT2 = Math.Min(minT2, projT2);
                maxT2 = Math.Max(maxT2, projT2);
            }

            double L1 = maxT1 - minT1; // Dimension along tangent1
            double L2 = maxT2 - minT2; // Dimension along tangent2

            return (L1, L2);
        }

        /// <summary>
        /// Add moment equilibrium constraints relating eccentricities to forces
        /// </summary>
        private void AddBiaxialBending(GRBModel model, Face face, GeometryModel geometry,
                                                       GRBVar fnTotal, GRBVar e1, GRBVar e2)
        {
            // Moment equilibrium about face centroid
            GRBLinExpr momentAboutT2Axis = new GRBLinExpr(); // Creates e1
            GRBLinExpr momentAboutT1Axis = new GRBLinExpr(); // Creates e2

            foreach (int vId in face.VertexIds)
            {
                int idx = GetPairIndex(face.Id, vId);
                if (idx < 0) continue;

                ContactPoint vertex = geometry.Vertices[vId];

                // Distance from face centroid to vertex
                double dx = vertex.X - face.Centroid[0];
                double dy = vertex.Y - face.Centroid[1];
                double dz = vertex.Z - face.Centroid[2];

                // Project distance onto tangent directions
                double distT1 = dx * face.Tangent1[0] + dy * face.Tangent1[1] + dz * face.Tangent1[2];
                double distT2 = dx * face.Tangent2[0] + dy * face.Tangent2[1] + dz * face.Tangent2[2];

                // Add to moment sums
                momentAboutT2Axis.AddTerm(distT1, fAll[3 * idx]); // Normal force × distance along T1
                momentAboutT1Axis.AddTerm(distT2, fAll[3 * idx]); // Normal force × distance along T2
            }

            // Moment equilibrium: Sum(fi × di) = fnTotal × eccentricity
            GRBQuadExpr momEq1 = new GRBQuadExpr();
            momEq1.Add(momentAboutT2Axis);
            momEq1.AddTerm(-1.0, fnTotal, e1);
            model.AddQConstr(momEq1 == 0.0, $"MomentEq1_face{face.Id}");

            GRBQuadExpr momEq2 = new GRBQuadExpr();
            momEq2.Add(momentAboutT1Axis);
            momEq2.AddTerm(-1.0, fnTotal, e2);
            model.AddQConstr(momEq2 == 0.0, $"MomentEq2_face{face.Id}");
        }

        /// <summary>
        /// Analyze biaxial results after optimization
        /// </summary>
        private void AddBiaxialBendingConstraints(GRBModel model, GeometryModel geometry, ProblemData data)
        {
            if (model.SolCount == 0) return;

            Console.WriteLine("\n=== BIAXIAL BENDING ANALYSIS ===");
            Console.WriteLine("Face | fnTotal |   e1   |   e2   | |e1|  | |e2|  | L1_eff | L2_eff | A_eff | Usage");
            Console.WriteLine("-----|---------|--------|--------|------|------|--------|--------|-------|------");

            foreach (var kvp in biaxialVars)
            {
                int faceId = kvp.Key;
                var vars = kvp.Value;
                double sigmaC = data.SigmaC;
                double fnTotal = vars.FnTotal.X;
                double e1 = vars.E1.X;
                double e2 = vars.E2.X;
                double e1_abs = vars.E1_Abs.X;
                double e2_abs = vars.E2_Abs.X;
                double L1_eff = vars.L1_Eff.X;
                double L2_eff = vars.L2_Eff.X;

                Face face = _geometry.Faces[faceId];
                double A_eff = L1_eff * L2_eff * face.Thickness;
                double maxForce = Data.SigmaC * A_eff;
                double usage = fnTotal / (maxForce + 1e-12);

                Console.WriteLine($"{faceId,4} | {fnTotal,7:F1} | {e1,6:F3} | {e2,6:F3} | " +
                                 $"{e1_abs,4:F3} | {e2_abs,4:F3} | {L1_eff,6:F3} | {L2_eff,6:F3} | " +
                                 $"{A_eff,5:F3} | {usage,5:F3}");

                // Flag critical conditions
                if (usage > 0.95)
                    Console.WriteLine($"      ⚠️  Face {faceId}: High biaxial usage {usage:F3}");
                if (e1_abs > 0.4 * (face.Area / face.Thickness) || e2_abs > 0.4 * (face.Area / face.Thickness))
                    Console.WriteLine($"      ⚠️  Face {faceId}: Large eccentricity detected");
            }
        }

        // =================================================================
        // UNDERSTANDING SOLUTION SPACE
        // =================================================================

        /*
        THIS APPROACH HELPS US UNDERSTAND:

        1. TYPICAL ECCENTRICITY VALUES:
           - Are e1, e2 usually small compared to L1/2, L2/2?
           - Is the bilinear term significant or nearly linear?

        2. CRITICAL FACES:
           - Which faces actually hit biaxial limits?
           - Are most faces uniaxial-dominated or truly biaxial?

        3. LINEARIZATION NEED:
           - If eccentricities are small, maybe linear approximation is sufficient
           - If only few faces are critical, maybe special handling for those only

        4. SOLVER PERFORMANCE:
           - How does Gurobi handle these quadratic constraints?
           - Are solve times acceptable?

        NEXT STEPS AFTER UNDERSTANDING SOLUTION:
        - If quadratic solving is too slow → consider linearization
        - If eccentricities are small → use linear approximation  
        - If only few faces critical → hybrid approach
        - If solution space is well-behaved → keep quadratic

        LET'S RUN THIS AND SEE WHAT THE SOLUTION SPACE LOOKS LIKE!
        */
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

                // 2) Save local (f_n, f_t1, f_t2) for each face–vertex pair
                writer.WriteLine("\n=== 3D Forces: FaceID; Pt(X,Y,Z); fn; ft1; ft2 ===");
                for (int j = 0; j < faceVertexPairs.Count; j++)
                {
                    int indexFn = 3 * j; 
                    int indexFt1 = 3 * j + 1;
                    int indexFt2 = 3 * j + 2;

                    double fnVal, ft1Val, ft2Val;
                    try
                    {
                        fnVal = fAll[indexFn].Get(GRB.DoubleAttr.X);
                        ft1Val = fAll[indexFt1].Get(GRB.DoubleAttr.X);
                        ft2Val = fAll[indexFt2].Get(GRB.DoubleAttr.X);
                    }
                    catch (GRBException e)
                    {
                        writer.WriteLine($"Could not retrieve forces for pair j={j}: {e.Message}");
                        continue;
                    }

                    var pair = faceVertexPairs[j]; // (faceId, vertexId)

                    // Get 3D vertex coordinates
                    double ptX = 0.0, ptY = 0.0, ptZ = 0.0;
                    if (_geometry.Vertices.ContainsKey(pair.VertexId))
                    {
                        ContactPoint vertex = _geometry.Vertices[pair.VertexId];
                        ptX = vertex.X;
                        ptY = vertex.Y;
                        ptZ = vertex.Z; // Added Z coordinate
                    }

                    writer.WriteLine($"{pair.FaceId}; {ptX:F3},{ptY:F3},{ptZ:F3}; {fnVal:F3}; {ft1Val:F3}; {ft2Val:F3}");
                }

                // 3) Save biaxial bending results if available
                if (biaxialVars.Count > 0)
                {
                    writer.WriteLine("\n=== Biaxial Bending Results ===");
                    writer.WriteLine("Face | fnTotal |   e1   |   e2   | |e1|  | |e2|  | L1_eff | L2_eff | A_eff | Usage");
                    writer.WriteLine("-----|---------|--------|--------|------|------|--------|--------|-------|------");

                    foreach (var kvp in biaxialVars)
                    {
                        int faceId = kvp.Key;
                        var vars = kvp.Value;

                        try
                        {
                            double fnTotal = vars.FnTotal.X;
                            double e1 = vars.E1.X;
                            double e2 = vars.E2.X;
                            double e1_abs = vars.E1_Abs.X;
                            double e2_abs = vars.E2_Abs.X;
                            double L1_eff = vars.L1_Eff.X;
                            double L2_eff = vars.L2_Eff.X;

                            Face face = _geometry.Faces[faceId];
                            double A_eff = L1_eff * L2_eff * face.Thickness;
                            double maxForce = data.SigmaC * A_eff;
                            double usage = fnTotal / (maxForce + 1e-12);

                            writer.WriteLine($"{faceId,4} | {fnTotal,7:F1} | {e1,6:F3} | {e2,6:F3} | " +
                                           $"{e1_abs,4:F3} | {e2_abs,4:F3} | {L1_eff,6:F3} | {L2_eff,6:F3} | " +
                                           $"{A_eff,5:F3} | {usage,5:F3}");
                        }
                        catch (GRBException e)
                        {
                            writer.WriteLine($"Face {faceId}: Could not retrieve biaxial variables: {e.Message}");
                        }
                    }
                }

                // 4) Compute and save total forces per face (3D)
                writer.WriteLine("\n=== Face Total Forces: Face, Total Fn, Total Ft1, Total Ft2, |Ft_resultant| ===");

                var faceTotalForces = new Dictionary<int, (double fnSum, double ft1Sum, double ft2Sum)>();

                for (int j = 0; j < faceVertexPairs.Count; j++)
                {
                    int indexFn = 3 * j;
                    int indexFt1 = 3 * j + 1;
                    int indexFt2 = 3 * j + 2;

                    double fnVal = fAll[indexFn].Get(GRB.DoubleAttr.X);
                    double ft1Val = fAll[indexFt1].Get(GRB.DoubleAttr.X);
                    double ft2Val = fAll[indexFt2].Get(GRB.DoubleAttr.X);

                    var pair = faceVertexPairs[j];

                    if (!faceTotalForces.ContainsKey(pair.FaceId))
                        faceTotalForces[pair.FaceId] = (0.0, 0.0, 0.0);

                    var current = faceTotalForces[pair.FaceId];
                    faceTotalForces[pair.FaceId] = (current.fnSum + fnVal, current.ft1Sum + ft1Val, current.ft2Sum + ft2Val);
                }
                // Write totals with resultant tangential force
                foreach (var kvp in faceTotalForces)
                {
                    double ft_resultant = Math.Sqrt(kvp.Value.ft1Sum * kvp.Value.ft1Sum + kvp.Value.ft2Sum * kvp.Value.ft2Sum);
                    writer.WriteLine($"{kvp.Key},{kvp.Value.fnSum:F3},{kvp.Value.ft1Sum:F3},{kvp.Value.ft2Sum:F3},{ft_resultant:F3}");
                }


                // 5) Shear limits at the FACE level (not vertex level)
                foreach (var kvp in faceTotalForces)
                {
                    int faceId = kvp.Key;
                    double totalNormal = kvp.Value.fnSum;
                    double totalTangent1 = kvp.Value.ft1Sum;
                    double totalTangent2 = kvp.Value.ft2Sum;
                    double tangentialResultant = Math.Sqrt(totalTangent1 * totalTangent1 + totalTangent2 * totalTangent2);

                    Face face = _geometry.Faces[faceId];
                    double mu = face.MuOverride ?? data.Mu;
                    double cohesion = face.CohesionValue;
                    double area = face.Area;
                    double frictionLimit = mu * totalNormal + cohesion * area;

                    double usageFriction = tangentialResultant / (frictionLimit + 1e-9);

                    string frictionStatus = usageFriction >= 0.995 ? "CRITICAL" : "OK";

                    writer.WriteLine($"Face {faceId} 3D Friction: ratio={usageFriction:F3}, status={frictionStatus} " +
                                   $"(|ft_res|={tangentialResultant:F3}, limit={frictionLimit:F3})");
                    writer.WriteLine($"  Components: ft1={totalTangent1:F3}, ft2={totalTangent2:F3}, fn={totalNormal:F3}");
                }

                writer.WriteLine("\n--- 3D Face Geometry: ID, Vertices, Centroid, Normal, Tangent1, Tangent2, Area ---");
                foreach (var fEntry in _geometry.Faces)
                {
                    int faceId = fEntry.Key;
                    Face face = fEntry.Value;

                    // Vertex coordinates
                    string vertexCoords = "";
                    foreach (int vId in face.VertexIds)
                    {
                        if (_geometry.Vertices.ContainsKey(vId))
                        {
                            ContactPoint cp = _geometry.Vertices[vId];
                            vertexCoords += $"({cp.X:F3},{cp.Y:F3},{cp.Z:F3}) ";
                        }
                    }

                    // Face vectors
                    string centroid = $"({face.Centroid[0]:F3},{face.Centroid[1]:F3},{face.Centroid[2]:F3})";
                    string normal = $"({face.Normal[0]:F3},{face.Normal[1]:F3},{face.Normal[2]:F3})";
                    string tangent1 = $"({face.Tangent1[0]:F3},{face.Tangent1[1]:F3},{face.Tangent1[2]:F3})";
                    string tangent2 = $"({face.Tangent2[0]:F3},{face.Tangent2[1]:F3},{face.Tangent2[2]:F3})";

                    writer.WriteLine($"Face {faceId}:");
                    writer.WriteLine($"  Vertices: {vertexCoords.Trim()}");
                    writer.WriteLine($"  Centroid: {centroid}");
                    writer.WriteLine($"  Normal:   {normal}");
                    writer.WriteLine($"  Tangent1: {tangent1}");
                    writer.WriteLine($"  Tangent2: {tangent2}");
                    writer.WriteLine($"  Area:     {face.Area:F3}");
                    writer.WriteLine();
                }

                // 7) 3D Block information
                writer.WriteLine("\n--- 3D Block Information: ID, Centroid(X,Y,Z), Applied Forces(Fx,Fy,Fz), Load Point(X,Y,Z) ---");
                foreach (var blockEntry in _geometry.Blocks.Values.Where(b => b.Id > 0).OrderBy(b => b.Id))
                {
                    writer.WriteLine($"Block {blockEntry.Id}:");
                    writer.WriteLine($"  Centroid:     ({blockEntry.CentroidX:F3},{blockEntry.CentroidY:F3},{blockEntry.CentroidZ:F3})");
                    writer.WriteLine($"  Applied F:    ({blockEntry.AppliedFx:F3},{blockEntry.AppliedFy:F3},{blockEntry.AppliedFz:F3})");
                    writer.WriteLine($"  Load Point:   ({blockEntry.LoadApplicationX:F3},{blockEntry.LoadApplicationY:F3},{blockEntry.LoadApplicationZ:F3})");

                    // Show computed moments
                    double xOffset = blockEntry.LoadApplicationX - blockEntry.CentroidX;
                    double yOffset = blockEntry.LoadApplicationY - blockEntry.CentroidY;
                    double zOffset = blockEntry.LoadApplicationZ - blockEntry.CentroidZ;
                    double computedMx = blockEntry.AppliedFy * zOffset - blockEntry.AppliedFz * yOffset;
                    double computedMy = blockEntry.AppliedFz * xOffset - blockEntry.AppliedFx * zOffset;
                    double computedMz = blockEntry.AppliedFx * yOffset - blockEntry.AppliedFy * xOffset;

                    writer.WriteLine($"  Computed M:   ({computedMx:F3},{computedMy:F3},{computedMz:F3})");
                    writer.WriteLine();
                }

            }
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
                LoadAllData(@"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\/data_cairo_friction_0e05.txt"   //data_pseudoparallel_friction_0e4
                , geometry, data);
                ComputeFaceNormalsFromGeometry(geometry);


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


        // Method to automatically construct B vector from block load data
        private static double[] ConstructBVectorFromBlocks(GeometryModel geometry)
        {
            var nonSupportBlocks = geometry.Blocks.Values
                .Where(b => b.Id > 0)
                .OrderBy(b => b.Id)
                .ToList();

            int numRows = nonSupportBlocks.Count * 3;
            double[] vectorB = new double[numRows];

            Console.WriteLine("\n=== AUTOMATIC 3D B VECTOR CONSTRUCTION ===");
            Console.WriteLine("Block | Applied Forces (Fx,Fy,Fz) | Load Point (X,Y,Z) | Centroid (X,Y,Z) | Computed Moments (Mx,My,Mz)");
            Console.WriteLine("------|---------------------------|-------------------|------------------|---------------------------");

            for (int blockIdx = 0; blockIdx < nonSupportBlocks.Count; blockIdx++)
            {
                var block = nonSupportBlocks[blockIdx];

                int baseRow = blockIdx * 6; // CHANGE from * 3 to * 6

                // Direct force assignments (first 3 rows)
                vectorB[baseRow] = block.AppliedFx;     // Fx
                vectorB[baseRow + 1] = block.AppliedFy; // Fy  
                vectorB[baseRow + 2] = block.AppliedFz; // Fz (NEW)

                // AUTOMATIC 3D MOMENT COMPUTATION (last 3 rows)
                // Offset from centroid to load application point
                double xOffset = block.LoadApplicationX - block.CentroidX;
                double yOffset = block.LoadApplicationY - block.CentroidY;
                double zOffset = block.LoadApplicationZ - block.CentroidZ;

                // 3D moment equations from applied forces:
                // Mx = Fy * zOffset - Fz * yOffset
                double computedMx = block.AppliedFy * zOffset - block.AppliedFz * yOffset;
                vectorB[baseRow + 3] = computedMx;

                // My = Fz * xOffset - Fx * zOffset  
                double computedMy = block.AppliedFz * xOffset - block.AppliedFx * zOffset;
                vectorB[baseRow + 4] = computedMy;

                // Mz = Fx * yOffset - Fy * xOffset (SAME as before)
                double computedMz = block.AppliedFx * yOffset - block.AppliedFy * xOffset;
                vectorB[baseRow + 5] = computedMz;

                // Display for verification
                string forces = $"({block.AppliedFx:F2},{block.AppliedFy:F2},{block.AppliedFz:F2})";
                string loadPoint = $"({block.LoadApplicationX:F2},{block.LoadApplicationY:F2},{block.LoadApplicationZ:F2})";
                string centroid = $"({block.CentroidX:F2},{block.CentroidY:F2},{block.CentroidZ:F2})";
                string moments = $"({computedMx:F3},{computedMy:F3},{computedMz:F3})";

                Console.WriteLine($"{block.Id,5} | {forces,25} | {loadPoint,17} | {centroid,16} | {moments}");

                // Show detailed calculation if there are applied forces
                if (Math.Abs(block.AppliedFx) + Math.Abs(block.AppliedFy) + Math.Abs(block.AppliedFz) > 1e-6)
                {
                    Console.WriteLine($"      Moment calculation:");
                    Console.WriteLine($"        Mx = Fy*zOffset - Fz*yOffset = {block.AppliedFy:F3}*{zOffset:F3} - {block.AppliedFz:F3}*{yOffset:F3} = {computedMx:F6}");
                    Console.WriteLine($"        My = Fz*xOffset - Fx*zOffset = {block.AppliedFz:F3}*{xOffset:F3} - {block.AppliedFx:F3}*{zOffset:F3} = {computedMy:F6}");
                    Console.WriteLine($"        Mz = Fx*yOffset - Fy*xOffset = {block.AppliedFx:F3}*{yOffset:F3} - {block.AppliedFy:F3}*{xOffset:F3} = {computedMz:F6}");
                }
            }

            Console.WriteLine($"\nGenerated 3D B vector with {numRows} entries (6 per block)");
            Console.WriteLine("All moments computed automatically from force application points and block centroids.");
            return vectorB;
        }

  



            static void LoadAllData(string filePath, GeometryModel geometry, ProblemData data)
        {
            string currentSection = "";
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

                        if (!int.TryParse(blockParts[0], out int blockId))
                        {
                            Console.WriteLine($"Line {lineNo}: Invalid BlockID '{blockParts[0]}'. Skipping block.");
                            continue;
                        }

                        // Parse block data: ID, CentroidX, CentroidY, CentroidZ, AppliedFx, AppliedFy, AppliedFz, LoadAppX, LoadAppY, LoadAppZ
                        Block block = new Block
                        {
                            Id = blockId,
                            CentroidX = double.Parse(blockParts[1]),
                            CentroidY = double.Parse(blockParts[2]),
                            CentroidZ = double.Parse(blockParts[3])     // ADD this line
                        };

                        // Check if load data is provided (optional - default to zero loads at centroid)
                        // Check if load data is provided (optional - default to zero loads at centroid)
                        if (blockParts.Length >= 10)  // CHANGE from >= 9 to >= 10
                        {
                            block.AppliedFx = double.Parse(blockParts[4]);
                            block.AppliedFy = double.Parse(blockParts[5]);
                            block.AppliedFz = double.Parse(blockParts[6]);     // ADD this line
                            block.LoadApplicationX = double.Parse(blockParts[7]);
                            block.LoadApplicationY = double.Parse(blockParts[8]);
                            block.LoadApplicationZ = double.Parse(blockParts[9]);  // ADD this line
                        }
                        else
                        {
                            // No load data, assume no applied loads
                            block.AppliedFx = 0.0;
                            block.AppliedFy = 0.0;
                            block.AppliedFz = 0.0;    // ADD this line
                            block.LoadApplicationX = block.CentroidX;
                            block.LoadApplicationY = block.CentroidY;
                            block.LoadApplicationZ = block.CentroidZ;  // ADD this line
                        }

                        geometry.Blocks.Add(blockId, block);
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
                            Y = double.Parse(vertexParts[2]),
                            Z = double.Parse(vertexParts[3])    // ADD this line
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

                    case "[VectorG]":
                        vectorG.AddRange(ParseLoadLine(trimmedLine, "VectorG", lineNo));
                        break;
                }
            }

            data.B = ConstructBVectorFromBlocks(geometry);
            data.G = vectorG.ToArray();

            // Validate VectorB and VectorG lengths
            int realBlocksCount = geometry.Blocks.Values.Count(b => b.Id > 0);
            int expectedRows = realBlocksCount * 6;

            if (data.B.Length != expectedRows)
                throw new InvalidDataException($"Generated B vector has {data.B.Length} entries, expected {expectedRows}");

            //ValidateVectorLength(vectorG, expectedRows, "VectorG");

            Console.WriteLine($"✓ Automatically generated B vector with {data.B.Length} entries");
            Console.WriteLine($"✓ Loaded {geometry.Blocks.Count} blocks with load information");
            Console.WriteLine($"✓ Loaded {geometry.Vertices.Count} vertices");
            Console.WriteLine($"✓ Loaded {geometry.Faces.Count} faces");
        }
        private static double DistanceBetween(ContactPoint a, ContactPoint b)
        {
            return Math.Sqrt(Math.Pow(b.X - a.X, 2) + Math.Pow(b.Y - a.Y, 2) + Math.Pow(b.Z - a.Z, 2));  // ADD + Math.Pow(b.Z - a.Z, 2)
        }


        public static void ComputeFaceNormalsFromGeometry(GeometryModel geometry)
        {
            foreach (var face in geometry.Faces.Values)
            {
                    if (face.VertexIds.Count < 3)
                    {
                        Console.WriteLine($"Face {face.Id} has {face.VertexIds.Count} vertices, expected at least 3 for 3D.");
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

                //  store the original IDs in local variables
                int jId = face.BlockJ;
                int kId = face.BlockK;

                // Retrieve the actual Block objects
                var jRef = geometry.Blocks[jId];
                var kRef = geometry.Blocks[kId];

                // Decide if we need to swap them
                bool needSwap = false;
                const double eps = 1e-12;

                // Compare centroidX first
                if (Math.Abs(jRef.CentroidX - kRef.CentroidX) < eps)
                {
                    // X are effectively the same, compare Y
                    if (Math.Abs(jRef.CentroidY - kRef.CentroidY) < eps)
                    {
                        // X and Y are the same, compare Z (3D extension)
                        if (jRef.CentroidZ > kRef.CentroidZ)
                            needSwap = true;
                    }
                    else
                    {
                        if (jRef.CentroidY > kRef.CentroidY)
                            needSwap = true;
                    }
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

                // Re-fetch references after potential swap
                jRef = geometry.Blocks[face.BlockJ];
                kRef = geometry.Blocks[face.BlockK];

                // STEP 2: Check that all vertices exist
                bool allVerticesExist = true;
                foreach (int vertexId in face.VertexIds)
                {
                    if (!geometry.Vertices.ContainsKey(vertexId))
                    {
                        Console.WriteLine($"Face {face.Id}: missing vertex {vertexId} in geometry!");
                        allVerticesExist = false;
                    }
                }
                if (!allVerticesExist) continue;

                // STEP 3: VERTEX ORDERING for consistent matrix computation
                // Ensure vertices are ordered so normal points from J → K
                var vertices = face.VertexIds.Select(vId => geometry.Vertices[vId]).ToList();

                // Use first three vertices to compute trial normal
                var v0 = vertices[0];
                var v1 = vertices[1];
                var v2 = vertices[2];

                // Edge vectors from v0
                double[] edge1 = { v1.X - v0.X, v1.Y - v0.Y, v1.Z - v0.Z };
                double[] edge2 = { v2.X - v0.X, v2.Y - v0.Y, v2.Z - v0.Z };

                // Trial normal = edge1 × edge2
                double[] trialNormal = {
            edge1[1] * edge2[2] - edge1[2] * edge2[1],  // nx
            edge1[2] * edge2[0] - edge1[0] * edge2[2],  // ny
            edge1[0] * edge2[1] - edge1[1] * edge2[0]   // nz
        };

                // Normalize trial normal
                double normalMag = Math.Sqrt(trialNormal[0] * trialNormal[0] +
                                           trialNormal[1] * trialNormal[1] +
                                           trialNormal[2] * trialNormal[2]);

                if (normalMag < 1e-12)
                {
                    Console.WriteLine($"Face {face.Id}: Cannot compute normal, vertices are collinear.");
                    continue;
                }

                trialNormal[0] /= normalMag;
                trialNormal[1] /= normalMag;
                trialNormal[2] /= normalMag;

                // STEP 4: YOUR ORIGINAL DIRECTION CHECK - extended to 3D
                // Ensure normal points from BlockJ toward BlockK
                double dCx = kRef.CentroidX - jRef.CentroidX;
                double dCy = kRef.CentroidY - jRef.CentroidY;
                double dCz = kRef.CentroidZ - jRef.CentroidZ;

                double dot = trialNormal[0] * dCx + trialNormal[1] * dCy + trialNormal[2] * dCz;

                if (dot < 0)
                {
                    // Flip normal and reverse vertex order (extending your original logic)
                    trialNormal[0] = -trialNormal[0];
                    trialNormal[1] = -trialNormal[1];
                    trialNormal[2] = -trialNormal[2];
                    face.VertexIds.Reverse();
                    vertices.Reverse();
                    Console.WriteLine($"Face {face.Id}: Reversed vertex order to ensure normal points J→K");
                }

                // STEP 5: TANGENT COMPUTATION
                // Tangent1: along first edge (normalized)
                vertices = face.VertexIds.Select(vId => geometry.Vertices[vId]).ToList(); // Refresh after potential reversal
                v0 = vertices[0];
                v1 = vertices[1];

                double[] firstEdge = { v1.X - v0.X, v1.Y - v0.Y, v1.Z - v0.Z };
                double firstEdgeMag = Math.Sqrt(firstEdge[0] * firstEdge[0] +
                                               firstEdge[1] * firstEdge[1] +
                                               firstEdge[2] * firstEdge[2]);

                double[] tangent1;
                if (firstEdgeMag > 1e-12)
                {
                    tangent1 = new double[] {
                firstEdge[0] / firstEdgeMag,
                firstEdge[1] / firstEdgeMag,
                firstEdge[2] / firstEdgeMag
            };
                }
                else
                {
                    // Fallback if first edge is degenerate
                    tangent1 = new double[] { 1, 0, 0 };
                }

                // Tangent2: Normal × Tangent1 (in face plane, perpendicular to Tangent1)
                double[] tangent2 = {
            trialNormal[1] * tangent1[2] - trialNormal[2] * tangent1[1],  // tx2
            trialNormal[2] * tangent1[0] - trialNormal[0] * tangent1[2],  // ty2  
            trialNormal[0] * tangent1[1] - trialNormal[1] * tangent1[0]   // tz2
        };

                // Normalize Tangent2
                double t2Mag = Math.Sqrt(tangent2[0] * tangent2[0] + tangent2[1] * tangent2[1] + tangent2[2] * tangent2[2]);
                if (t2Mag > 1e-12)
                {
                    tangent2[0] /= t2Mag;
                    tangent2[1] /= t2Mag;
                    tangent2[2] /= t2Mag;
                }

                // STEP 6: Store results in face
                face.Normal = trialNormal;
                face.Tangent1 = tangent1;
                face.Tangent2 = tangent2;
                face.Tangent = tangent1;  // For compatibility with existing code

                // STEP 7: COMPUTE FACE AREA using simple triangulation from first vertex
                ComputeFaceArea(face, geometry);

                // STEP 8: Compute face centroid
                ComputeFaceCentroid(face, geometry);

                // Verification output
                Console.WriteLine($"3D Face {face.Id}: {vertices.Count} vertices, Area={face.Area:F3}, " +
                                 $"Normal=({trialNormal[0]:F3},{trialNormal[1]:F3},{trialNormal[2]:F3})");

                // Optional: Verify orthogonality
                VerifyFaceVectors(face);
            }
        }

        /// <summary>
        /// Verify that computed vectors are orthogonal and unit length
        /// </summary>
        private static void VerifyFaceVectors(Face face)
        {
            const double tolerance = 1e-6;

            // Check magnitudes
            double nMag = Math.Sqrt(face.Normal[0] * face.Normal[0] + face.Normal[1] * face.Normal[1] + face.Normal[2] * face.Normal[2]);
            double t1Mag = Math.Sqrt(face.Tangent1[0] * face.Tangent1[0] + face.Tangent1[1] * face.Tangent1[1] + face.Tangent1[2] * face.Tangent1[2]);
            double t2Mag = Math.Sqrt(face.Tangent2[0] * face.Tangent2[0] + face.Tangent2[1] * face.Tangent2[1] + face.Tangent2[2] * face.Tangent2[2]);

            // Check orthogonality
            double dot_n_t1 = face.Normal[0] * face.Tangent1[0] + face.Normal[1] * face.Tangent1[1] + face.Normal[2] * face.Tangent1[2];
            double dot_n_t2 = face.Normal[0] * face.Tangent2[0] + face.Normal[1] * face.Tangent2[1] + face.Normal[2] * face.Tangent2[2];
            double dot_t1_t2 = face.Tangent1[0] * face.Tangent2[0] + face.Tangent1[1] * face.Tangent2[1] + face.Tangent1[2] * face.Tangent2[2];

            bool hasIssues = false;
            if (Math.Abs(nMag - 1.0) > tolerance)
            {
                Console.WriteLine($"  Warning: Face {face.Id} normal magnitude = {nMag:F6}, expected 1.0");
                hasIssues = true;
            }
            if (Math.Abs(t1Mag - 1.0) > tolerance)
            {
                Console.WriteLine($"  Warning: Face {face.Id} tangent1 magnitude = {t1Mag:F6}, expected 1.0");
                hasIssues = true;
            }
            if (Math.Abs(t2Mag - 1.0) > tolerance)
            {
                Console.WriteLine($"  Warning: Face {face.Id} tangent2 magnitude = {t2Mag:F6}, expected 1.0");
                hasIssues = true;
            }
            if (Math.Abs(dot_n_t1) > tolerance)
            {
                Console.WriteLine($"  Warning: Face {face.Id} normal·tangent1 = {dot_n_t1:F6}, expected 0.0");
                hasIssues = true;
            }
            if (Math.Abs(dot_n_t2) > tolerance)
            {
                Console.WriteLine($"  Warning: Face {face.Id} normal·tangent2 = {dot_n_t2:F6}, expected 0.0");
                hasIssues = true;
            }
            if (Math.Abs(dot_t1_t2) > tolerance)
            {
                Console.WriteLine($"  Warning: Face {face.Id} tangent1·tangent2 = {dot_t1_t2:F6}, expected 0.0");
                hasIssues = true;
            }

            if (!hasIssues)
            {
                Console.WriteLine($"  ✓ Face {face.Id}: All vectors are orthonormal");
            }
        }

        /// <summary>
        /// Compute face area using simple triangulation from first vertex
        /// </summary>
        private static void ComputeFaceArea(Face face, GeometryModel geometry)
        {
            var vertices = face.VertexIds.Select(vId => geometry.Vertices[vId]).ToList();

            if (vertices.Count < 3)
            {
                face.Area = 0.0;
                face.Length = 0.0;
                return;
            }

            double totalArea = 0.0;
            var v0 = vertices[0];  // First vertex (common vertex for all triangles)

            // Triangulate: v0-v1-v2, v0-v2-v3, v0-v3-v4, etc.
            for (int i = 1; i < vertices.Count - 1; i++)
            {
                var v1 = vertices[i];
                var v2 = vertices[i + 1];

                // Triangle area = 0.5 * |edge1 × edge2|
                double[] edge1 = { v1.X - v0.X, v1.Y - v0.Y, v1.Z - v0.Z };
                double[] edge2 = { v2.X - v0.X, v2.Y - v0.Y, v2.Z - v0.Z };

                // Cross product
                double[] cross = {
            edge1[1] * edge2[2] - edge1[2] * edge2[1],
            edge1[2] * edge2[0] - edge1[0] * edge2[2],
            edge1[0] * edge2[1] - edge1[1] * edge2[0]
        };

                // Magnitude of cross product
                double crossMag = Math.Sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
                totalArea += 0.5 * crossMag;
            }

            face.Area = totalArea;
            face.Length = Math.Sqrt(totalArea);  // Approximate "length" for compatibility

            Console.WriteLine($"Face {face.Id}: Computed area = {totalArea:F6}");
        }

        /// <summary>
        /// Compute face centroid
        /// </summary>
        private static void ComputeFaceCentroid(Face face, GeometryModel geometry)
        {
            var vertices = face.VertexIds.Select(vId => geometry.Vertices[vId]).ToList();

            double sumX = vertices.Sum(v => v.X);
            double sumY = vertices.Sum(v => v.Y);
            double sumZ = vertices.Sum(v => v.Z);

            face.Centroid = new double[] {
        sumX / vertices.Count,
        sumY / vertices.Count,
        sumZ / vertices.Count
    };
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
