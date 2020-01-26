using Numerics;
using Numerics.Objects;
using System.Collections.Generic;

namespace System
{
    public static class NumericExtensions
    {
        public const double h = 0.0000001;

        public static double Derivate(this Func<double, double> func, double x) => (func(x + h) - func(x)) / h;

        public static double Derivate(this Func<double[], double> func, double[] variables, int index)
        {
            var variablescopy = (double[])variables.Clone();
            variables[index] += h;
            return (func(variables) - func(variablescopy)) / h;
        }

        public static double Derivate(this Func<Vector, double> func, Vector variables, Cartesian cartesian)
        {
           var delta = variables;
            switch (cartesian)
            {
                case Cartesian.x:
                    delta = new Vector(variables.x+h, variables.y, variables.z);
                    break;
                case Cartesian.y:
                    delta = new Vector(variables.x, variables.y+h, variables.z);
                    break;
                case Cartesian.z:
                    delta = new Vector(variables.x, variables.y, variables.z + h);
                    break;
            }

            return (func(delta) - func(variables)) / h;
        }
        public static IEnumerable<Vector> Gradient(this Func<Vector, double> func, double xMin, double yMin, double zMin, double stepSize, double maxSteps)
        {
            var vectorField = new List<Vector>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                vectorField.Add(func.Gradient((xMin + i, yMin + i, zMin + i)));
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

        public static double Integrate(this Func<double, double> func, double lowerlimit, double upperlimit)
        {
            var result = 0.0;
            for (var i = lowerlimit; i <= upperlimit; i += h)
            {
                result += h * (func(i) + func(i + h)) / 2;
            }
            return result;
        }
    }
}

