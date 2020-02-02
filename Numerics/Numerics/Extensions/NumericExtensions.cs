using Numerics;
using Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace System
{
    public static class NumericExtensions
    {
        public const double h = 0.000001;


        public static double Derivate(this Func<double, double> func, double x, int order =1) {

            var h0 = Math.Pow(10, order) * h;
            var d = 0.0;
            var pascalMatrix = new Matrix(new double[order + 1, order + 1]).Pascal();
            for (var i =1;i<=order+1; i++)
            {
       
                d += Sign(i) * pascalMatrix.values[order, i-1] * func(x + (order - i) * h0);

            }

            return d / Math.Pow(h0, order);
        } 

        public static double Derivate(this Func<double[], double> func, double[] variables, int index, int order =1)
        {
            Func<double, double> funcIndex = (double v) => func(variables.Select((c,i)=> { if (i == index) { return v; } return c; }).ToArray());
            return funcIndex.Derivate(variables[index], order);
      
        }

        public static double Derivate(this Func<Vector, double> func, Vector variables, Cartesian cartesian, int order=1)
        {
      
            switch (cartesian)
            {
                case Cartesian.x:
                   Func<double, double> funcX = (double v) => func(new Vector(v, variables.y, variables.z));
                   return funcX.Derivate(variables.x, order);
                 
                case Cartesian.y:
                    Func<double, double> funcY = (double v) => func(new Vector(variables.x, v, variables.z));
                    return funcY.Derivate(variables.y, order);
                case Cartesian.z:
                    Func<double, double> funcZ= (double v) => func(new Vector(variables.x, variables.y, v));
                    return funcZ.Derivate(variables.z, order);
        
            }

            return func(variables);
        }



        public static double Derivate(this Func<(double x,double y), double> func, (double x, double y) variables, Cartesian cartesian, int order = 1)
        {

            switch (cartesian)
            {
                case Cartesian.x:
                    Func<double, double> funcX = (double v) => func((v, variables.y));
                    return funcX.Derivate(variables.x, order);

                case Cartesian.y:
                    Func<double, double> funcY = (double v) => func((variables.x, v));
                    return funcY.Derivate(variables.y, order);

            }

            return func(variables);
        }


        public static ComplexNumber Derivate(this ComplexFunction func, ComplexNumber variables, int order = 1)
        {

            var du= func.u.Derivate((variables.realPart, variables.imaginaryPart),Cartesian.x,order);

            var dv = func.v.Derivate((variables.realPart, variables.imaginaryPart), Cartesian.x, order);
          
            return (new ComplexNumber(du, dv));
  
        }


        public static IDictionary<Vector, Vector> Gradient(this Func<Vector, double> func, double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var vectorField = new Dictionary<Vector, Vector>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                vectorField.Add(new Vector((xmin + i, ymin + i, zmin + i)),func.Gradient((xmin + i, ymin + i, zmin + i)));
            }
            return vectorField;

        }

        public static Vector Gradient(this Func<Vector, double> func, (double,double,double) points)
        {
            var dx = func.Derivate(new Vector(points), Cartesian.x);
            var dy = func.Derivate(new Vector(points), Cartesian.y);
            var dz = func.Derivate(new Vector(points), Cartesian.z);
            return new Vector(dx, dy, dz);

        }

        public static double Laplacian(this Func<Vector, double> func, (double, double, double) points)
        {
            var dx2 = func.Derivate(new Vector(points), Cartesian.x, 2);
            var dy2 = func.Derivate(new Vector(points), Cartesian.y,2);
            var dz2 = func.Derivate(new Vector(points), Cartesian.z,2);
            return dx2+ dy2+ dz2;

        }

        public static double Integrate(this Func<double, double> func, double lowerLimit, double upperLimit)
        {
            var result = 0.0;
            for (var i = lowerLimit; i <= upperLimit; i += h)
            {
                result += h * (func(i) + func(i + h)) / 2;
            }
            return result;
        }


        public static double Integrate(this Func<(double x, double y), double> func, (double lowerLimit, double upperLimit) xlimit, (double lowerLimit, double upperLimit) ylimit)
        {

           var rnd = new Random();

            var result = 0.0;
            var throws = 999999;


            for (var i = 0; i < throws; i ++)
                {

                var x = rnd.NextDouble() *(xlimit.upperLimit - xlimit.lowerLimit) + xlimit.lowerLimit;
                var y = rnd.NextDouble() * (ylimit.upperLimit - ylimit.lowerLimit) + ylimit.lowerLimit; 

                result +=func((x,y));

            }
            return (xlimit.upperLimit - xlimit.lowerLimit) * (ylimit.upperLimit - ylimit.lowerLimit) * result / 999999; 
         
        }



        private static int Sign(int index) => index % 2 == 0 ? -1 : 1;
    }
}

