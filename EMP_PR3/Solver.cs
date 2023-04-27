using System.ComponentModel.Design;
using System.Diagnostics;
using System.Numerics;

namespace EMP_PR3;

public abstract class Solver
{
   protected double _eps = 1e-14;
   protected int _maxIters = 1000;

   protected SparseMatrix _matrix;
   protected Vector _vector;
   protected Vector _solution;
   public double SolveTime { get; protected set; }

   public Solver()
   {
      _matrix = new SparseMatrix(0, 0);
      _vector = new Vector(0);
      _solution = new Vector(0);
   }

   public void SetSLAE(Vector vector, SparseMatrix matrix)
   {
      _vector = vector;
      _matrix = matrix;
   }

   public abstract Vector Solve();

   protected void LU()
   {
      for (int i = 0; i < _matrix.Size; i++)
      {

         for (int j = _matrix._ia[i]; j < _matrix._ia[i + 1]; j++)
         {
            int jCol = _matrix._ja[j];
            int jk = _matrix._ia[jCol];
            int k = _matrix._ia[i];

            int sdvig = _matrix._ja[_matrix._ia[i]] - _matrix._ja[_matrix._ia[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < _matrix._ia[jCol + 1]; k++, jk++)
            {
               sumL += _matrix._al[k] * _matrix._au[jk];
               sumU += _matrix._au[k] * _matrix._al[jk];
            }

            _matrix._al[j] -= sumL;
            _matrix._au[j] -= sumU;
            _matrix._au[j] /= _matrix._di[jCol];
         }

         double sumD = 0.0;
         for (int j = _matrix._ia[i]; j < _matrix._ia[i + 1]; j++)
            sumD += _matrix._al[j] * _matrix._au[j];

         _matrix._di[i] -= sumD;
      }
   }

   protected void ForwardElimination()
   {
      for (int i = 0; i < _matrix.Size; i++)
      {
         for (int j = _matrix._ia[i]; j < _matrix._ia[i + 1]; j++)
         {
            _solution[i] -= _matrix._al[j] * _solution[_matrix._ja[j]];
         }

         _solution[i] /= _matrix._di[i];
      }
   }

   protected void BackwardSubstitution()
   {
      for (int i = _matrix.Size - 1; i >= 0; i--)
      {
         for (int j = _matrix._ia[i + 1] - 1; j >= _matrix._ia[i]; j--)
         {
            _solution[_matrix._ja[j]] -= _matrix._au[j] * _solution[i];
         }
      }
   }

   public void PrintSolution()
   {
      for(int i = 0; i < _solution.Size; i++)
      {
         Console.WriteLine(_solution[i]);
      }
   }
}

public class LOS : Solver
{
   public override Vector Solve()
   {
      _solution = new(_vector.Size);
      Vector.Copy(_vector, _solution);

      Vector r = _vector - _matrix * _solution;
      Vector z = 1 * r;
      Vector p = _matrix * z;
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double discrepancy = r * r;

      for (int i = 1; i <= _maxIters && discrepancy > _eps; i++)
      {
         alpha = (p * r) / (p * p);
         _solution += alpha * z;

         r -= alpha * p;
         tmp = _matrix * r;

         beta = -(p * tmp) / (p * p);

         z = r + beta * z;
         p = tmp + beta * p;

         discrepancy = r * r;
      }

      sw.Stop();
      SolveTime = sw.ElapsedMilliseconds;

      return _solution;
   }
}

public class LU : Solver
{
   public override Vector Solve()
   {
      _solution = new(_vector.Size);
      Vector.Copy(_vector, _solution);
      _matrix = _matrix.ConvertToProfile();

      Stopwatch sw = Stopwatch.StartNew();

      LU();
      ForwardElimination();
      BackwardSubstitution();

      sw.Stop();
      SolveTime = sw.ElapsedMilliseconds;

      return _solution;
   }
}

public class LOSWithLU : Solver
{
   public override Vector Solve() 
   {
      _solution = new(_vector.Size);
      Vector.Copy(_vector, _solution);

      SparseMatrix matrixLU = new(_matrix.Size, _matrix._ja.Length);
      SparseMatrix.Copy(_matrix, matrixLU);

      //matrixLU = matrixLU.ConvertToProfile();

      LU(matrixLU);

      Vector r = DirElim(matrixLU, _vector - _matrix * _solution);
      Vector z = BackSub(matrixLU, r);
      Vector p = DirElim(matrixLU, _matrix * z);
      Vector tmp;
      double alpha;
      double beta;

      Stopwatch sw = Stopwatch.StartNew();

      double discrepancy = r * r;


      for (int i = 1; i <= _maxIters && discrepancy > _eps; i++)
      {
         alpha = (p * r) / (p * p);
         _solution += alpha * z;

         r -= alpha * p;

         tmp = DirElim(matrixLU, _matrix * BackSub(matrixLU, r));
         beta = -(p * tmp) / (p * p);

         z = BackSub(matrixLU, r) + beta * z;
         p = tmp + beta * p;

         discrepancy = r * r;
      }

      sw.Stop();
      SolveTime = sw.ElapsedMilliseconds;

      return _solution;
   }

   protected static void LU(SparseMatrix Matrix)
   {
      for (int i = 0; i < Matrix.Size; i++)
      {

         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
         {
            int jCol = Matrix._ja[j];
            int jk = Matrix._ia[jCol];
            int k = Matrix._ia[i];

            int sdvig = Matrix._ja[Matrix._ia[i]] - Matrix._ja[Matrix._ia[jCol]];

            if (sdvig > 0)
               jk += sdvig;
            else
               k -= sdvig;

            double sumL = 0.0;
            double sumU = 0.0;

            for (; k < j && jk < Matrix._ia[jCol + 1]; k++, jk++)
            {
               sumL += Matrix._al[k] * Matrix._au[jk];
               sumU += Matrix._au[k] * Matrix._al[jk];
            }

            Matrix._al[j] -= sumL;
            Matrix._au[j] -= sumU;
            Matrix._au[j] /= Matrix._di[jCol];
         }

         double sumD = 0.0;
         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
            sumD += Matrix._al[j] * Matrix._au[j];

         Matrix._di[i] -= sumD;
      }
   }

   protected static Vector DirElim(SparseMatrix Matrix, Vector b)
   {
      Vector result = new Vector(b.Size);
      Vector.Copy(b, result);

      for (int i = 0; i < Matrix.Size; i++)
      {
         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
         {
            result[i] -= Matrix._al[j] * result[Matrix._ja[j]];
         }

         result[i] /= Matrix._di[i];
      }

      return result;
   }

   protected static Vector BackSub(SparseMatrix Matrix, Vector b)
   {
      var result = new Vector(b.Size);
      Vector.Copy(b, result);

      for (int i = Matrix.Size - 1; i >= 0; i--)
      {
         for (int j = Matrix._ia[i + 1] - 1; j >= Matrix._ia[i]; j--)
         {
            result[Matrix._ja[j]] -= Matrix._au[j] * result[i];
         }
      }

      return result;
   }

}

//public class BCG : Solver
//{
//   public override Vector Solve()
//   {
//      //_solution = new(_vector.Size);
//      //Vector.Copy(_vector, _solution);

//      //SparseMatrix matrixLU = new(_matrix.Size, _matrix._ja.Length);
//      //SparseMatrix.Copy(_matrix, matrixLU);

//      //LU(matrixLU);

//      //return _solution;


//      _solution = new(_vector.Size);
//      Vector.Copy(_vector, _solution);

//      Vector r = _vector - _matrix * _solution;

//      Vector p = new(r.Size);
//      Vector z = new(r.Size);
//      Vector s = new(r.Size);

//      Vector.Copy(r, p);
//      Vector.Copy(r, z);
//      Vector.Copy(r, s);


//      Stopwatch sw = Stopwatch.StartNew();

//      double discrepancy = r * r;

//      for (int i = 1; i <= _maxIters && discrepancy > _eps; i++)
//      {
//         var Az = (_matrix * z);
//         double alpha = (p * r) / (s * Az);

//         _solution = _solution + alpha * z;
//         r = r - alpha * Az;
//         p = p - 

//      }

//      sw.Stop();
//      SolveTime = sw.ElapsedMilliseconds;

//      return _solution;
//   }

//   protected static void LU(SparseMatrix Matrix)
//   {
//      for (int i = 0; i < Matrix.Size; i++)
//      {

//         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
//         {
//            int jCol = Matrix._ja[j];
//            int jk = Matrix._ia[jCol];
//            int k = Matrix._ia[i];

//            int sdvig = Matrix._ja[Matrix._ia[i]] - Matrix._ja[Matrix._ia[jCol]];

//            if (sdvig > 0)
//               jk += sdvig;
//            else
//               k -= sdvig;

//            double sumL = 0.0;
//            double sumU = 0.0;

//            for (; k < j && jk < Matrix._ia[jCol + 1]; k++, jk++)
//            {
//               sumL += Matrix._al[k] * Matrix._au[jk];
//               sumU += Matrix._au[k] * Matrix._al[jk];
//            }

//            Matrix._al[j] -= sumL;
//            Matrix._au[j] -= sumU;
//            Matrix._au[j] /= Matrix._di[jCol];
//         }

//         double sumD = 0.0;
//         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
//            sumD += Matrix._al[j] * Matrix._au[j];

//         Matrix._di[i] -= sumD;
//      }
//   }

//   protected static Vector DirElim(SparseMatrix Matrix, Vector b)
//   {
//      Vector result = new Vector(b.Size);
//      Vector.Copy(b, result);

//      for (int i = 0; i < Matrix.Size; i++)
//      {
//         for (int j = Matrix._ia[i]; j < Matrix._ia[i + 1]; j++)
//         {
//            result[i] -= Matrix._al[j] * result[Matrix._ja[j]];
//         }

//         result[i] /= Matrix._di[i];
//      }

//      return result;
//   }

//   protected static Vector BackSub(SparseMatrix Matrix, Vector b)
//   {
//      Vector result = new Vector(b.Size);
//      Vector.Copy(b, result);

//      for (int i = Matrix.Size - 1; i >= 0; i--)
//      {
//         for (int j = Matrix._ia[i + 1] - 1; j >= Matrix._ia[i]; j--)
//         {
//            result[Matrix._ja[j]] -= Matrix._au[j] * result[i];
//         }
//      }
//      return result;
//   }

//}