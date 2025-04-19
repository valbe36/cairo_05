// MatrixBuilder.cs  –  compile into its own small exe
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using static System.Math;

// MatrixBuilder.cs
// -------------------------------------------------------------------
// build equilibrium matrix A, external load vector B, and self-weight
// vector G from a single txt file with headers: vertices / faces /
// blocks / bvector.  Outputs matrix_A_cairo.txt that your solver
// already reads.
// -------------------------------------------------------------------



namespace MatrixBuilderApp
{
    class Program
    {
        // -----------------------------------------------------------------
        // simple data records
        // -----------------------------------------------------------------
        class Vertex { public int Id; public double X, Y; }
        class Face
        {
            public int Id, Left, Right, V1, V2;
            public double t, c, mu;
        }
        class Block { public int Id; public double Cx, Cy, W; }


        private const string DEFAULT_INPUT = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\input_cairo.txt";
        private const string DEFAULT_OUTPUT = @"C:\Users\vb\OneDrive - Aarhus universitet\Dokumenter 1\work research\54 ICSA\JOURNAL paper\analyses\output_cairo.txt";


        // -----------------------------------------------------------------
        static void Main(string[] args)
        {
            string inFile = args.Length > 0 ? args[0] : DEFAULT_INPUT;
            string outFile = args.Length > 1 ? args[1] : DEFAULT_OUTPUT;

            if (!File.Exists(inFile))
            {
                Console.WriteLine($"Input file '{inFile}' not found."); return;
            }

            // 1) read all sections
            var vertices = new Dictionary<int, Vertex>();
            var faces = new List<Face>();
            var blocks = new List<Block>();
            var Bvector = new List<double>();

            ReadSections(inFile, vertices, faces, blocks, Bvector);

            // --- sanity: every face vertex must exist ---
            var missing = new List<(int faceId, int vId)>();
            foreach (var f in faces)
            {
                if (!vertices.ContainsKey(f.V1)) missing.Add((f.Id, f.V1));
                if (!vertices.ContainsKey(f.V2)) missing.Add((f.Id, f.V2));
            }
            if (missing.Count > 0)
            {
                Console.WriteLine("ERROR: Faces refer to vertex IDs that were not declared:");
                foreach (var (fid, vid) in missing)
                    Console.WriteLine($"   face {fid}  →  vertex {vid}");
                return;                             // stop building – fix the input
            }

            // 2) orient every face so that (v1→v2) is CCW around Left block
            foreach (var f in faces) OrientFace(f, vertices, blocks);

            // 3) build matrices
            BuildMatrices(faces, vertices, blocks, Bvector.ToArray(),
                          out double[,] A, out double[] B, out double[] G);

            // 4) write file in solver format
            WriteOutput(outFile, A, B, G);
            Console.WriteLine($"Matrix written to '{outFile}'.");
        }

        // -----------------------------------------------------------------
        // read sections (vertices / faces / blocks / bvector)
        // -----------------------------------------------------------------
        static void ReadSections(string path,
                                 Dictionary<int, Vertex> V,
                                 List<Face> F,
                                 List<Block> B,
                                 List<double> BV)
       
        {
            string section = "";
            foreach (string raw in File.ReadLines(path))
            {
                string line = raw.Split('#')[0].Trim();
                if (line == "") continue;

                // header lines
                if (line.Equals("vertices", StringComparison.OrdinalIgnoreCase) ||
                    line.Equals("faces", StringComparison.OrdinalIgnoreCase) ||
                    line.Equals("blocks", StringComparison.OrdinalIgnoreCase) ||
                    line.Equals("bvector", StringComparison.OrdinalIgnoreCase))
                { section = line.ToLower(); continue; }

                string[] p = line.Split(',', StringSplitOptions.RemoveEmptyEntries)
                                 .Select(s => s.Trim()).ToArray();

                switch (section)
                {
                    case "vertices":
                        V[int.Parse(p[0])] = new Vertex
                        {
                            Id = int.Parse(p[0]),
                            X = double.Parse(p[1], CultureInfo.InvariantCulture),
                            Y = double.Parse(p[2], CultureInfo.InvariantCulture)
                        };
                        break;

                    case "faces":
                        F.Add(new Face
                        {
                            Id = int.Parse(p[0]),
                            Left = int.Parse(p[1]),
                            Right = int.Parse(p[2]),
                            t = double.Parse(p[3], CultureInfo.InvariantCulture),
                            c = double.Parse(p[4], CultureInfo.InvariantCulture),
                            mu = double.Parse(p[5], CultureInfo.InvariantCulture),
                            V1 = int.Parse(p[6]),
                            V2 = int.Parse(p[7])
                        });
                        break;

                    case "blocks":
                        B.Add(new Block
                        {
                            Id = int.Parse(p[0]),
                            Cx = double.Parse(p[1], CultureInfo.InvariantCulture),
                            Cy = double.Parse(p[2], CultureInfo.InvariantCulture),
                            W = double.Parse(p[3], CultureInfo.InvariantCulture)
                        });
                        break;

                    case "bvector":
                        BV.Add(double.Parse(p[0], CultureInfo.InvariantCulture));
                        break;
                }
            }
        }

        // -----------------------------------------------------------------
        // Ensure (v1→v2) is counter‑clockwise around Left block
        // -----------------------------------------------------------------
        static void OrientFace(Face f,
                               Dictionary<int, Vertex> V,
                               List<Block> B)
        {
            // if Left is a support (0 or -1) but Right is a real block,
            // swap them so that Left always refers to a block in the list
            if (f.Left <= 0 && f.Right > 0)
            {
                // swap block IDs
                int tmp = f.Left; f.Left = f.Right; f.Right = tmp;
                // also swap vertex order to keep edge direction
                int tv = f.V1; f.V1 = f.V2; f.V2 = tv;
            }

            // edge vector e
            double ex = V[f.V2].X - V[f.V1].X;
            double ey = V[f.V2].Y - V[f.V1].Y;
            // candidate outward normal for current ordering
            double nx = ey, ny = -ex;

            var bLeft = B.First(b => b.Id == f.Left);
            var bRight = f.Right <= 0
                         ? new Block { Cx = bLeft.Cx + nx, Cy = bLeft.Cy + ny }
                         : B.First(b => b.Id == f.Right);

            // vector Left → Right
            double dx = bRight.Cx - bLeft.Cx;
            double dy = bRight.Cy - bLeft.Cy;

            // if normal points inside Left block, swap
            if (nx * dx + ny * dy <= 0)
            {
                int tmp = f.V1; f.V1 = f.V2; f.V2 = tmp;
                int t2 = f.Left; f.Left = f.Right; f.Right = t2;
            }
        }

        // -----------------------------------------------------------------
        // Build A, B, G
        // -----------------------------------------------------------------
        static void BuildMatrices(List<Face> faces,
                                  Dictionary<int, Vertex> V,
                                  List<Block> blocks,
                                  double[] BvecInput,
                                  out double[,] A,
                                  out double[] B,
                                  out double[] G)
        {
            int nPairs = faces.Count * 2;
            int nRows = 3 * blocks.Count;
            int nCols = 2 * nPairs;
            A = new double[nRows, nCols];
            B = new double[nRows];
            G = new double[nRows];

            // copy optional B‑vector values
            for (int i = 0; i < Math.Min(BvecInput.Length, B.Length); i++)
                B[i] = BvecInput[i];

            // column map (faceId,vertexId) -> column index
            var colMap = new Dictionary<(int, int), int>();
            int col = 0;
            foreach (var f in faces)
            {
                colMap[(f.Id, f.V1)] = col++;
                colMap[(f.Id, f.V2)] = col++;
            }

            // row index per block
            var rowOf = blocks.Select((b, idx) => (b.Id, idx))
                              .ToDictionary(p => p.Id, p => 3 * p.idx);

            foreach (var f in faces)
            {
                int[] vIds = { f.V1, f.V2 };
                var vPos = vIds.Select(id => V[id]).ToArray();

                // outward normal of LEFT
                double ex = vPos[1].X - vPos[0].X;
                double ey = vPos[1].Y - vPos[0].Y;
                double len = Sqrt(ex * ex + ey * ey);
                double nx = ey / len, ny = -ex / len;
                double tx = -ny, ty = nx;

                foreach ((int blkId, double sign) in new[] { (f.Left, +1.0),
                                                             (f.Right, -1.0) })
                {
                    if (blkId <= 0) continue;                  // support side

                    int r0 = rowOf[blkId];
                    var blk = blocks.First(b => b.Id == blkId);

                    for (int k = 0; k < 2; k++)                 // vertices
                    {
                        int vId = vIds[k];
                        int cIdx = colMap[(f.Id, vId)];
                        int colN = 2 * cIdx;
                        int colT = colN + 1;

                        double rx = vPos[k].X - blk.Cx;
                        double ry = vPos[k].Y - blk.Cy;
                        double zN = rx * ny - ry * nx;
                        double zT = rx * ty - ry * tx;

                        A[r0, colN] += sign * nx;
                        A[r0, colT] += sign * tx;
                        A[r0 + 1, colN] += sign * ny;
                        A[r0 + 1, colT] += sign * ty;
                        A[r0 + 2, colN] += sign * zN;
                        A[r0 + 2, colT] += sign * zT;
                    }
                    // self‑weight (Fx = 0, Fy = −W, M = 0)
                    G[r0 + 1] += blk.W;
                }
            }
        }

        // -----------------------------------------------------------------
        // Write output file in solver format
        // -----------------------------------------------------------------
        static void WriteOutput(string path, double[,] A,
                                double[] B, double[] G)
        {
            int nRows = A.GetLength(0), nCols = A.GetLength(1);
            using var w = new StreamWriter(path);
            w.WriteLine("matrix_A_size");
            w.WriteLine(nRows);
            w.WriteLine(nCols);
            w.WriteLine("matrix_A_values");
            for (int r = 0; r < nRows; r++)
            {
                for (int c = 0; c < nCols; c++)
                    w.Write(A[r, c].ToString("G17",
                        CultureInfo.InvariantCulture) + " ");
                w.WriteLine();
            }
            w.WriteLine("vector_B_values");
            foreach (var x in B)
                w.Write(x.ToString("G17", CultureInfo.InvariantCulture) + " ");
            w.WriteLine();
            w.WriteLine("vector_G_values");
            foreach (var x in G)
                w.Write(x.ToString("G17", CultureInfo.InvariantCulture) + " ");
            w.WriteLine();
        }
    }
}