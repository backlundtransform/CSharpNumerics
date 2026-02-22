using System.Collections.Generic;


namespace System.Linq
{
    public static class DescriptiveStatisticsExtensions
    {
        /// <summary>
        /// Computes the median value of a sequence after projecting each element to a scalar.
        /// For an even number of elements, returns the average of the two middle values.
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to numeric value.</param>
        /// <returns>The median of the projected values, or 0 if the sequence is empty.</returns>
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

        /// <summary>
        /// Computes the population variance of a sequence after projecting each element to a scalar.
        /// variance = sum((x - mean)^2) / n.
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to numeric value.</param>
        /// <returns>The variance of the projected values.</returns>
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

        /// <summary>
        /// Computes the standard deviation of a sequence after projecting each element to a scalar.
        /// StandardDeviation = sqrt(Variance).
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to numeric value.</param>
        /// <returns>The standard deviation of the projected values.</returns>
        public static double StandardDeviation<T>(this IEnumerable<T> enumerable, Func<T, double> func) => Math.Sqrt(enumerable.Variance(func));

        /// <summary>
        /// Computes the sample covariance between two variables.
        /// cov(X,Y) = sum((x-meanX)(y-meanY)) / (n-1).
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to (x,y).</param>
        /// <returns>The sample covariance.</returns>
        public static double Covariance<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var xy = enumerable.Select(func).ToList();

            var meanX = xy.Average(p => p.x);
            var meanY = xy.Average(p => p.y);

            return xy.Sum(p => (p.x - meanX) * (p.y - meanY)) / (xy.Count() - 1);
        }

        /// <summary>
        /// Computes the coefficient of determination (R²) for paired data.
        /// This implementation computes correlation = cov / (stdDevX * stdDevY) and returns R² = correlation².
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to (x,y).</param>
        /// <returns>R² (in [0,1] for typical data), or 0 if any standard deviation is zero.</returns>
        public static double CoefficientOfDetermination<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var xy = enumerable.Select(func);
            var cov = xy.Covariance((p) => (p.x, p.y));

            var stdDevX = xy.Select(p => p.x).StandardDeviation();
            var stdDevY = xy.Select(p => p.y).StandardDeviation();


            if (stdDevX == 0 || stdDevY == 0)
            {
                return 0.0;
            }


            double correlation = cov / (stdDevY * stdDevX);

            return correlation * correlation;
        }

        /// <summary>
        /// Computes the standard deviation of a sequence of doubles.
        /// If <paramref name="sample"/> is true, uses the sample standard deviation (divisor n-1);
        /// otherwise uses population standard deviation (divisor n).
        /// </summary>
        /// <param name="values">The input values.</param>
        /// <param name="sample">True for sample standard deviation; false for population.</param>
        /// <returns>The standard deviation (0 when undefined under current rules).</returns>
        public static double StandardDeviation(this IEnumerable<double> values, bool sample = true)
        {
            var list = values.ToList();
            int n = list.Count;
            if (n == 0 || (sample && n < 2))
                return 0.0;

            double mean = list.Average();
            double sumOfSquares = list.Sum(x => Math.Pow(x - mean, 2));

            double divisor = sample ? (n - 1) : n;
            return Math.Sqrt(sumOfSquares / divisor);
        }

        /// <summary>
        /// Computes confidence intervals for the mean using an approximate z-score.
        /// Returns (mean - z * σ / sqrt(n), mean + z * σ / sqrt(n)).
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to numeric value.</param>
        /// <param name="confidenceLevel">Confidence level (e.g., 0.95).</param>
        /// <returns>Tuple containing (lowerValue, upperValue).</returns>
        public static (double lowerValue, double upperValue) ConfidenceIntervals<T>(this IEnumerable<T> enumerable, Func<T, double> func, double confidenceLevel) {

            var standardDeviation = enumerable.StandardDeviation(func);

            var mean = enumerable.Average(func);

            var count = enumerable.Count();

            var zIndex = 0.0;

            for (var i = 4.0; i >= 0; i -= 0.001)
            {
                var area = 1 - (1 - Math.Exp(-1.98 * i / Math.Sqrt(2))) * Math.Exp(-Math.Pow(i / Math.Sqrt(2), 2)) / (1.135 * Math.Sqrt(Math.PI) * i / Math.Sqrt(2));

                if (Math.Round(area, 3) ==confidenceLevel)
                {
                    zIndex = i;
                    break;
                }

            }
            var interval = zIndex * standardDeviation / Math.Sqrt(count);

            return (mean-interval, mean+ interval);

        }

        /// <summary>
        /// Computes the cumulative sum of a projected sequence.
        /// Returns a sequence where each element is the sum of all previous projected values.
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to numeric value.</param>
        /// <returns>An enumerable of cumulative sums.</returns>
        public static IEnumerable<double> CumulativeSum<T>(this IEnumerable<T> enumerable, Func<T, double> func)
        {
            double sum = 0;

            var sequence = enumerable.Select(func);

            foreach (var item in sequence)
            {
                sum += item;
                yield return sum;
            }
        }
     

    }
}