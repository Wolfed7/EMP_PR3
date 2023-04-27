namespace EMP_PR3;

public abstract class Test
{
   protected double omega;
   protected double lambda;
   protected double sigma;
   protected double chi;

   public Test(Mesh3D mesh)
   {
      omega = mesh.Omega;
      lambda = mesh.Lambda;
      sigma = mesh.Sigma;
      chi = mesh.Chi;
   }

   public abstract double Us(Point3D point);
   protected abstract double DivGradUs(Point3D point);
   public abstract double Uc(Point3D point);
   protected abstract double DivGradUc(Point3D point);
   public double Fs(Point3D point) 
      => -lambda * DivGradUs(point) - omega * sigma * Uc(point) - omega * omega * chi * Us(point);
   public double Fc(Point3D point)
      => -lambda * DivGradUc(point) + omega * sigma * Us(point) - omega * omega * chi * Uc(point);
}

public class Test1 : Test
{
   public Test1(Mesh3D mesh) : base(mesh) { }

   public override double Us(Point3D point) 
      => 2 * point.X + point.Y / 2 + 5 * point.Z;
   public override double Uc(Point3D point) 
      => point.Z - 1.5 * point.X - point.Y;
   protected override double DivGradUs(Point3D point)
      => 0;
   protected override double DivGradUc(Point3D point)
      => 0;
}