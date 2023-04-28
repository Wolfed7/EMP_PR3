using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EMP_PR3;

public class FEM
{
   public delegate double Basis(Point3D point);

   private Mesh3D _mesh;
   private SparseMatrix _globalMatrix;
   private Vector _globalVector;
   private Vector _solution;
   private Solver _slae;
   private Vector _localVector, _localVector1, _localVector2;
   private Matrix _stiffnessMatrix;
   private Matrix _massMatrix;

   public static int NodesPerElement => 8;

   public Test TestCase;

   public FEM(Mesh3D mesh)
   {
      _mesh = mesh;
      _stiffnessMatrix = new(NodesPerElement);
      _massMatrix = new(NodesPerElement);
      _localVector = new(2 * NodesPerElement);
      _localVector1 = new(NodesPerElement);
      _localVector2 = new(NodesPerElement);

      _globalMatrix = new SparseMatrix(0, 0);
      _globalVector = new(0);
      _solution = new(0);
      _slae = new LOSWithLU();

   }

   public void SetTest(Test test)
      => TestCase = test;

   public void SetSolver(Solver slae)
      => _slae = slae;


   public void Compute()
   {
      BuildPortrait();

      AssemblySlae();
      AccountFirstConditions();

      _slae.SetSLAE(_globalVector, _globalMatrix);
      _solution = _slae.Solve();
   }

   public void BuildPortrait()
   {
      var list = new HashSet<int>[2 * _mesh.Nodes.Length].Select(_ => new HashSet<int>()).ToList();
      foreach (var element in _mesh.Elements)
      {
         foreach (var position in element)
         {
            foreach (var node in element)
            {
               if (position == node)
               {
                  list[2 * position + 1].Add(2 * position);
               }
               else if (position > node)
               {
                  list[2 * position].Add(2 * node);
                  list[2 * position].Add(2 * node + 1);

                  list[2 * position + 1].Add(2 * node);
                  list[2 * position + 1].Add(2 * node + 1);
               }
            }
         }
      }

      list = list.Select(childlist => childlist.Order().ToHashSet()).ToList();
      int offDiagonalElementsCount = list.Sum(childList => childList.Count);

      _globalMatrix = new(2 * _mesh.Nodes.Length, offDiagonalElementsCount);
      _globalVector = new(2 * _mesh.Nodes.Length);

      _globalMatrix._ia[0] = 0;

      for (int i = 0; i < list.Count; i++)
         _globalMatrix._ia[i + 1] = _globalMatrix._ia[i] + list[i].Count;

      int k = 0;
      foreach (var childList in list)
      {
         foreach (var value in childList)
         {
            _globalMatrix._ja[k++] = value;
         }
      }
   }

   private void AssemblySlae()
   {
      _globalVector.Fill(0);
      _globalMatrix.Clear();

      for (int ielem = 0; ielem < _mesh.Elements.Length; ielem++)
      {
         AssemblyLocalSLAE(ielem);

         for (int i = 0; i < NodesPerElement; i++)
            for (int j = 0; j < NodesPerElement; j++)
            {
               AddElement(2 * _mesh.Elements[ielem][i], 2 * _mesh.Elements[ielem][j], _stiffnessMatrix[i, j]);
               AddElement(2 * _mesh.Elements[ielem][i] + 1, 2 * _mesh.Elements[ielem][j] + 1, _stiffnessMatrix[i, j]);
               AddElement(2 * _mesh.Elements[ielem][i], 2 * _mesh.Elements[ielem][j] + 1, -_massMatrix[i, j]);
               AddElement(2 * _mesh.Elements[ielem][i] + 1, 2 * _mesh.Elements[ielem][j], _massMatrix[i, j]);
            }

         AddElementToVector(ielem);

         _stiffnessMatrix.Clear();
         _massMatrix.Clear();
         _localVector1.Fill(0);
         _localVector2.Fill(0);
      }
   }

   private void AssemblyLocalSLAE(int ielem)
   {
      int Mu(int i) => i % 2;
      int Nu(int i) => i / 2 % 2;
      int Theta(int i) => i / 4;

      double hx = Math.Abs(_mesh.Nodes[_mesh.Elements[ielem][7]].X - _mesh.Nodes[_mesh.Elements[ielem][0]].X);
      double hy = Math.Abs(_mesh.Nodes[_mesh.Elements[ielem][7]].Y - _mesh.Nodes[_mesh.Elements[ielem][0]].Y);
      double hz = Math.Abs(_mesh.Nodes[_mesh.Elements[ielem][7]].Z - _mesh.Nodes[_mesh.Elements[ielem][0]].Z);

      double[,] matrixG = 
      { 
         { 1.0, -1.0 },
         { -1.0, 1.0 } 
      };

      double[,] matrixM = 
      {
         { 2.0 / 6.0, 1.0 / 6.0 }, 
         { 1.0 / 6.0, 2.0 / 6.0 } 
      };

      for (int i = 0; i < NodesPerElement; i++)
      {
         for (int j = 0; j < NodesPerElement; j++)
         {
            _stiffnessMatrix[i, j] = 
               matrixG[Mu(i), Mu(j)] / hx * matrixM[Nu(i), Nu(j)] * hy * matrixM[Theta(i), Theta(j)] * hz +
               matrixM[Mu(i), Mu(j)] * hx * matrixG[Nu(i), Nu(j)] / hy * matrixM[Theta(i), Theta(j)] * hz +
               matrixM[Mu(i), Mu(j)] * hx * matrixM[Nu(i), Nu(j)] * hy * matrixG[Theta(i), Theta(j)] / hz;

            _massMatrix[i, j] = matrixM[Mu(i), Mu(j)] * hx * matrixM[Nu(i), Nu(j)] * hy * matrixM[Theta(i), Theta(j)] * hz;
         }
      }

      for (int i = 0; i < NodesPerElement; i++)
      {
         _localVector1[i] = TestCase.Fs(_mesh.Nodes[_mesh.Elements[ielem][i]]);
         _localVector2[i] = TestCase.Fc(_mesh.Nodes[_mesh.Elements[ielem][i]]);
      }

      _localVector1 = _massMatrix * _localVector1;
      _localVector2 = _massMatrix * _localVector2;

      for (int i = 0; i < NodesPerElement; i++)
      {
         _localVector[2 * i] = _localVector1[i];
         _localVector[2 * i + 1] = _localVector2[i];
      }

      _stiffnessMatrix = _mesh.Lambda * _stiffnessMatrix + (-_mesh.Omega) * _mesh.Omega * _mesh.Chi * _massMatrix;

      _massMatrix = _mesh.Omega * _mesh.Sigma * _massMatrix;
   }

   private void AddElement(int i, int j, double value)
   {
      if (i == j)
      {
         _globalMatrix._di[i] += value;
         return;
      }

      if (i > j)
         for (int icol = _globalMatrix._ia[i]; icol < _globalMatrix._ia[i + 1]; icol++)
         {
            if (_globalMatrix._ja[icol] == j)
            {
               _globalMatrix._al[icol] += value;
               return;
            }
         }
      else
         for (int icol = _globalMatrix._ia[j]; icol < _globalMatrix._ia[j + 1]; icol++)
         {
            if (_globalMatrix._ja[icol] == i)
            {
               _globalMatrix._au[icol] += value;
               return;
            }
         }
   }

   private void AddElementToVector(int ielem)
   {
      for (int i = 0; i < NodesPerElement; i++)
      {
         _globalVector[2 * _mesh.Elements[ielem][i]] += _localVector[2 * i];
         _globalVector[2 * _mesh.Elements[ielem][i] + 1] += _localVector[2 * i + 1];
      }
   }

   public void AccountFirstConditions()
   {
      foreach (var node in _mesh.Boundaries)
      {
         int row = 2 * node;

         _globalMatrix._di[row] = 1;
         _globalVector[row] = TestCase.Us(_mesh.Nodes[node]);

         for (int i = _globalMatrix._ia[row]; i < _globalMatrix._ia[row + 1]; i++)
         {
            _globalMatrix._al[i] = 0;
         }

         for (int col = row + 1; col < _globalMatrix.Size; col++)
         {
            for (int j = _globalMatrix._ia[col]; j < _globalMatrix._ia[col + 1]; j++)
               if (_globalMatrix._ja[j] == row)
               {
                  _globalMatrix._au[j] = 0;
                  break;
               }
         }

         row = 2 * node + 1;

         _globalMatrix._di[row] = 1;
         _globalVector[row] = TestCase.Uc(_mesh.Nodes[node]);

         for (int i = _globalMatrix._ia[row]; i < _globalMatrix._ia[row + 1]; i++)
         {
            _globalMatrix._al[i] = 0;
         }

         for (int col = row + 1; col < _globalMatrix.Size; col++)
         {
            for (int j = _globalMatrix._ia[col]; j < _globalMatrix._ia[col + 1]; j++)
               if (_globalMatrix._ja[j] == row)
               {
                  _globalMatrix._au[j] = 0;
                  break;
               }
         }
      }
   }

   public void PrintSolution(StreamWriter sw)
   {
      Vector exactSolution = new(2 * _mesh.Nodes.Length);
      for (int i = 0; i < exactSolution.Size / 2; i++)
      {
         exactSolution[2 * i] = TestCase.Us(_mesh.Nodes[i]);
         exactSolution[2 * i + 1] = TestCase.Uc(_mesh.Nodes[i]);
      }
      Vector inaccuracySin = new(_mesh.Nodes.Length);
      Vector inaccuracyCos = new(_mesh.Nodes.Length);
      Vector inaccuracy = new(2 * _mesh.Nodes.Length);
      for (int i = 0; i < _mesh.Nodes.Length; i++)
      {
         if (!_mesh.Boundaries.Contains(i / 2))
         {
            inaccuracySin[i] = exactSolution[2 * i] - _solution[2 * i];
            inaccuracyCos[i] = exactSolution[2 * i + 1] - _solution[2 * i + 1];
         }
      }

      for (int i = 0; i < inaccuracy.Size; i++)
      {
         if (!_mesh.Boundaries.Contains(i / 2))
            inaccuracy[i] = exactSolution[i] - _solution[i];
      }

      //sw.WriteLine("Относительная погрешность U");
      //sw.WriteLine($"{inaccuracy.Norm() / exactSolution.Norm()}");
      //sw.WriteLine("Относительная погрешность Us");
      //sw.WriteLine($"{inaccuracySin.Norm() / exactSolution.Norm()}");
      //sw.WriteLine("Относительная погрешность Uc");
      //sw.WriteLine($"{inaccuracyCos.Norm() / exactSolution.Norm()}");
      //sw.WriteLine("Время");
      //sw.WriteLine($"{_slae.SolveTime / 1000.0} сек.");
      //sw.WriteLine();


      sw.Write($"{_slae.SolveTime / 1000.0} ");
      sw.Write($"{inaccuracy.Norm() / exactSolution.Norm()} ");
      sw.Write($"{inaccuracySin.Norm() / exactSolution.Norm()} ");
      sw.Write($"{inaccuracyCos.Norm() / exactSolution.Norm()} ");
      sw.WriteLine();


      Console.WriteLine($"{_slae.SolveTime / 1000.0} сек.");
   }
}

