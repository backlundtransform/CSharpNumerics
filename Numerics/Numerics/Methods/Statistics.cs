using System;

namespace Numerics.Methods
{
    public class Statistics
    {

        public static Func<double, double> NormalDistribution(double standardDeviation, double mean)
        {
            double func(double t) => 1 / (standardDeviation * Math.Sqrt(2 * Math.PI)) * Math.Exp(-Math.Pow((t-mean) / Math.Sqrt(standardDeviation), 2)/2);
            return func;
        }

    }
}
