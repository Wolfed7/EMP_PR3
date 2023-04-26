using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EMP_PR3;

public class Point3D : ICloneable
{
   public double X { get; init; }
   public double Y { get; init; }
   public double Z { get; init; }

   public Point3D(double x = 0, double y = 0, double z = 0)
   {
      X = x;
      Y = y;
      Z = z;
   }

   public object Clone()
   {
      return new Point3D(X, Y, Z);
   }

   public override string ToString()
   {
      X.ToString();
      Y.ToString();
      Z.ToString();
      return X + " " + Y + " " + Z;
   }

   public static Point3D Parse(string data)
   {
      var dataList = data.Split().Select(double.Parse).ToList();
      return new Point3D(dataList[0], dataList[1], dataList[2]);
   }

   public static Point3D operator +(Point3D p1, Point3D p2)
   {
      return new Point3D(p1.X + p2.X, p1.Y + p2.Y, p1.Z + p2.Z);
   }

   public static Point3D operator -(Point3D p1, Point3D p2)
   {
      return new Point3D(p1.X - p2.X, p1.Y - p2.Y, p1.Z - p2.Z);
   }
}
