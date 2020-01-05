using System;

namespace Numerics
{
    public static class NumericExtensions
    {

        public static double Derivate(this Func<double, double> func, double variablevalue) => (func(variablevalue + 0.0000001) - func(variablevalue)) / 0.0000001;

    }
}
