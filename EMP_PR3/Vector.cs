namespace EMP_PR3;

public class Vector
{
   private double[] _container;
   public int Size { get; init; }

   public Vector(int size)
   {
      _container = new double[size];
      Size = size;
   }

   public double this[int index]
   {
      get => _container[index];
      set => _container[index] = value;
   }

   public void Fill(double value)
   {
      for (int i = 0; i < Size; i++)
      {
         _container[i] = value;
      }
   }

   public void Clear(double value)
   {
      for (int i = 0; i < Size; i++)
      {
         _container[i] = 0;
      }
   }

   public static void Copy(Vector source, Vector destination)
   {
      for (int i = 0; i < source.Size; i++)
         destination[i] = source[i];
   }

   public void Norming()
   {
      double norm = Norm();

      for (int i = 0; i < Size; i++)
         _container[i] /= norm;
   }

   public double Norm()
   {
      double result = 0;

      for (int i = 0; i < Size; i++)
         result += _container[i] * _container[i];

      return Math.Sqrt(result);
   }

   public static Vector operator *(Matrix matrix, Vector vector)
   {
      Vector result = new(vector._container.Length);

      for (int i = 0; i < vector.Size; i++)
         for (int j = 0; j < vector.Size; j++)
            result._container[i] += matrix[i, j] * vector._container[j];

      return result;
   }

   public static Vector operator -(Vector fstVector, Vector sndVector)
   {
      Vector result = new(fstVector.Size);

      for (int i = 0; i < fstVector.Size; i++)
         result[i] = fstVector._container[i] - sndVector._container[i];

      return result;
   }

   public static Vector operator +(Vector fstVector, Vector sndVector)
   {
      Vector result = new(fstVector.Size);

      for (int i = 0; i < fstVector.Size; i++)
         result[i] = fstVector._container[i] + sndVector._container[i];

      return result;
   }

   public static Vector operator *(double coef, Vector vector)
   {
      Vector result = new(vector.Size);

      for (int i = 0; i < vector.Size; i++)
         result[i] = vector._container[i] * coef;

      return result;
   }

   public static double operator *(Vector fstVector, Vector sndVector)
   {
      double result = 0;

      for (int i = 0; i < fstVector.Size; i++)
         result += fstVector._container[i] * sndVector._container[i];

      return result;
   }
}