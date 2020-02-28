using System.Collections.Generic;

namespace System.Linq
{
    public static class StatisticsExtensions
    {
        public static double Median<T>(this IEnumerable<T> enumerable, Func<T, double> func)
        {
            var count = enumerable.Count();

            if (count == 0)
            {
                return 0;
            }

            var list = enumerable.OrderBy(func).Select(func);
            var midpoint = count / 2;

            if (count % 2 == 0)
            {
                return (Convert.ToDouble(list.ElementAt(midpoint - 1)) + Convert.ToDouble(list.ElementAt(midpoint))) / 2.0;
            }

            return Convert.ToDouble(list.ElementAt(midpoint));
        }

        public static double Variance<T>(this IEnumerable<T> enumerable, Func<T, double> func)
        {
            var count = enumerable.Count();
            var average = enumerable.Average(func);

            var sum = 0.0;

            foreach (var value in enumerable.Select(func))
            {
                sum += Math.Pow(value - average, 2);
            }

            return sum / count;
        }

        public static double StandardDeviation<T>(this IEnumerable<T> enumerable, Func<T, double> func) => Math.Sqrt(enumerable.Variance(func));

        public static double Covariance<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var xy = enumerable.Select(func);

            var meanX = xy.Average(p => p.x);
            var meanY = xy.Average(p => p.y);

            return xy.Sum(p => (p.x - meanX) * (p.y - meanY)) / (xy.Count() - 1);
        }


        public static double GenerateNoise(this Random random, double variance)
        {
            double normalDistribution(double t) =>1/(Math.Sqrt(2* variance * Math.PI))* Math.Exp(-Math.Pow(t / Math.Sqrt(variance),2));
            var foo = normalDistribution(1);

            return normalDistribution(random.Next((int)variance));
        }

        public static double GetRandomDouble(this Random random, double minimum, double maximum)
        {
            return random.NextDouble() * (maximum - minimum) + minimum;
        }

    }
}