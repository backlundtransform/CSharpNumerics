using System;

namespace Numerics
{
    public static class NumericExtensions
    {
        public const double h = 0.0000001;

        public static double Derivate(this Func<double, double> func, double x) => (func(x + h) - func(x)) / h;

        public static double PartialDerivate(this Func<double[], double> func, double[] variables, int index)
        {
            var variablescopy = (double[])variables.Clone();
            variables[index] += h;
            return (func(variables) - func(variablescopy)) / h;
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

