using EMP_PR3;
using System.Globalization;

CultureInfo.CurrentCulture = new CultureInfo("en-US");
const string file1 = "Area/AreaDescription.txt";

var mesh = new Mesh3D();
mesh.Input(file1);
mesh.BuildMesh();

double[] lambdas = { 1e2, 1e3, 1e4, 1e5, 8 * 1e5 };
double[] omegas = { 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };
double[] sigmas = { 0, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8 };
double[] chis = { 8.81 * 1e-12, 1e-11, 1e-10 };

//var sw = new StreamWriter($"Result/lambda_{mesh.Nodes.Length}.txt");
//for (int i = 0; i < lambdas.Length; i++)
//{
//   mesh.Omega = 1;
//   mesh.Sigma = 1;
//   mesh.Chi = 1;

//   mesh.Lambda = lambdas[i];
//   sw.WriteLine($"Lambda = {lambdas[i]}");

//   FEM fem = new(mesh);
//   fem.SetTest(new Test1(mesh));

//   //sw.WriteLine($"By LU");
//   //fem.SetSolver(new LU());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By LOS");
//   //fem.SetSolver(new LOS());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By BCG");
//   //fem.SetSolver(new BCG());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   sw.WriteLine($"By BCGWithLU");
//   fem.SetSolver(new BCGWithLU());
//   fem.Compute();
//   fem.PrintSolution(sw);
//}
//sw.Close();

//sw = new StreamWriter($"Result/omega_{mesh.Nodes.Length}.txt");
//for (int i = 0; i < omegas.Length; i++)
//{
//   mesh.Lambda = 1;
//   mesh.Sigma = 1;
//   mesh.Chi = 1;

//   mesh.Omega = omegas[i];
//   sw.WriteLine($"Omega = {omegas[i]}");

//   FEM fem = new(mesh);
//   fem.SetTest(new Test1(mesh));

//   //sw.WriteLine($"By LU");
//   //fem.SetSolver(new LU());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By LOS");
//   //fem.SetSolver(new LOS());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By BCG");
//   //fem.SetSolver(new BCG());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   sw.WriteLine($"By BCGWithLU");
//   fem.SetSolver(new BCGWithLU());
//   fem.Compute();
//   fem.PrintSolution(sw);
//}
//sw.Close();

//sw = new StreamWriter($"Result/sigma_{mesh.Nodes.Length}.txt");
//for (int i = 0; i < sigmas.Length; i++)
//{
//   mesh.Lambda = 1;
//   mesh.Omega = 1;
//   mesh.Chi = 1;

//   mesh.Sigma = sigmas[i];
//   sw.WriteLine($"sigma = {sigmas[i]}");

//   FEM fem = new(mesh);
//   fem.SetTest(new Test1(mesh));

//   //sw.WriteLine($"By LU");
//   //fem.SetSolver(new LU());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By LOS");
//   //fem.SetSolver(new LOS());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By BCG");
//   //fem.SetSolver(new BCG());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   sw.WriteLine($"By BCGWithLU");
//   fem.SetSolver(new BCGWithLU());
//   fem.Compute();
//   fem.PrintSolution(sw);
//}
//sw.Close();

//sw = new StreamWriter($"Result/chi_{mesh.Nodes.Length}.txt");
//for (int i = 0; i < chis.Length; i++)
//{
//   mesh.Lambda = 1;
//   mesh.Omega = 1;
//   mesh.Sigma = 1;

//   mesh.Chi = chis[i];
//   sw.WriteLine($"chi = {chis[i]}");

//   FEM fem = new(mesh);
//   fem.SetTest(new Test1(mesh));

//   //sw.WriteLine($"By LU");
//   //fem.SetSolver(new LU());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By LOS");
//   //fem.SetSolver(new LOS());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   //sw.WriteLine($"By BCG");
//   //fem.SetSolver(new BCG());
//   //fem.Compute();
//   //fem.PrintSolution(sw);

//   sw.WriteLine($"By BCGWithLU");
//   fem.SetSolver(new BCGWithLU());
//   fem.Compute();
//   fem.PrintSolution(sw);
//}
//sw.Close();


//var bcg = new BCG();
var bcg = new BCGWithLU();

//var m = new SparseMatrix(5, 10)
//{
//   _au = new double[] { 10, 0, 0, 0, 10, 10, 0, 0, 10, 0 },
//   _di = new double[] { 1, 2, 3, 4, 5 },
//   _al = new double[] { 0, 0, 11, 11, 0, 0, 0, 0, 11, 0 },
//   _ia = new int[] { 0, 0, 1, 3, 6, 10 },
//   _ja = new int[] { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 }
//};

//var v = new Vector(5);
//v[0] = 1;
//v[1] = 2;
//v[2] = 3;
//v[3] = 4;
//v[4] = 5;

var m = new SparseMatrix(5, 28)
{
   _au = new double[] { 10, 0, 0, 0, 10, 10, 0, 0, 10, 0, 0, 10, 10, 0, 0, 0, 10, 0, 0, 10, 0, 0, 0, 10, 0, 10, 0, 0 },
   _di = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 },
   _al = new double[] { 0, 0, 11, 11, 0, 0, 0, 11, 11, 0, 0, 11, 0, 0, 0, 11, 11, 0, 11, 0, 0, 0, 0, 11, 0, 11, 0, 0 },
   _ia = new int[] { 0, 0, 1, 3, 6, 10, 15, 21, 28 },
   _ja = new int[] { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6 }
};

var v = new Vector(8);
v[0] = 1;
v[1] = 2;
v[2] = 3;
v[3] = 4;
v[4] = 5;
v[5] = 6;
v[6] = 7;
v[7] = 8;

var answer = m * v;
//var answer = SparseMatrix.TransposedMatrixMult(m, v);
//for (int i = 0; i < answer.Size; i++)
//{
//   Console.WriteLine(answer[i]);
//}


bcg.SetSLAE(answer, m);
bcg.Solve();
bcg.PrintSolution();