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

var sw = new StreamWriter($"Result/lambda_{mesh.Nodes.Length}.txt");
for (int i = 0; i < lambdas.Length; i++)
{
   mesh.Omega = 1;
   mesh.Sigma = 1;
   mesh.Chi = 1;

   mesh.Lambda = lambdas[i];
   sw.WriteLine($"Lambda = {lambdas[i]}");

   FEM fem = new(mesh);
   fem.SetTest(new Test1(mesh));

   sw.WriteLine($"By LU");
   fem.SetSolver(new LU());
   fem.Compute();
   fem.PrintSolution(sw);

   sw.WriteLine($"By LOS");
   fem.SetSolver(new LOS());
   fem.Compute();
   fem.PrintSolution(sw);
}
sw.Close();

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

//   sw.WriteLine($"By LU");
//   fem.SetSolver(new LU());
//   fem.Compute();
//   fem.PrintSolution(sw);

//   sw.WriteLine($"By LOS");
//   fem.SetSolver(new LOS());
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

//   sw.WriteLine($"By LU");
//   fem.SetSolver(new LU());
//   fem.Compute();
//   fem.PrintSolution(sw);

//   sw.WriteLine($"By LOS");
//   fem.SetSolver(new LOS());
//   fem.Compute();
//   fem.PrintSolution(sw);
//}
//sw.Close();

sw = new StreamWriter($"Result/chi_{mesh.Nodes.Length}.txt");
for (int i = 0; i < chis.Length; i++)
{
   mesh.Lambda = 1;
   mesh.Omega = 1;
   mesh.Sigma = 1;

   mesh.Chi = chis[i];
   sw.WriteLine($"chi = {chis[i]}");

   FEM fem = new(mesh);
   fem.SetTest(new Test1(mesh));

   sw.WriteLine($"By LU");
   fem.SetSolver(new LU());
   fem.Compute();
   fem.PrintSolution(sw);

   sw.WriteLine($"By LOS");
   fem.SetSolver(new LOS());
   fem.Compute();
   fem.PrintSolution(sw);
}
sw.Close();