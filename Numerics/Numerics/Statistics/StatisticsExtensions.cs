using CSharpNumerics.Statistics.Methods;
using System.Collections.Generic;
using CSharpNumerics.Statistics.Data;

namespace System.Linq
{
    public static class DataExtensions
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

        /// <summary>
        /// Generates noise from a normal distribution with the given variance.
        /// </summary>
        /// <param name="random">Random number generator.</param>
        /// <param name="variance">Variance (σ²).</param>
        /// <returns>A random sample.</returns>
        public static double GenerateNoise(this Random random, double variance)
        {
            var normalDistribution = Statistics.NormalDistribution(Math.Sqrt(variance), 0);

            return normalDistribution(random.Next((int)variance));
        }

        /// <summary>
        /// Generates a uniformly distributed random double in [minimum, maximum).
        /// </summary>
        /// <param name="random">Random number generator.</param>
        /// <param name="minimum">Lower bound.</param>
        /// <param name="maximum">Upper bound.</param>
        /// <returns>A random value in the specified range.</returns>
        public static double GetRandomDouble(this Random random, double minimum, double maximum)
        {
            return random.NextDouble() * (maximum - minimum) + minimum;
        }

        /// <summary>
        /// Computes simple linear regression parameters for paired data.
        /// Returns slope, intercept, and correlation.
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to (x,y).</param>
        /// <returns>(slope, intercept, correlation).</returns>
        public static (double slope, double intercept, double correlation) LinearRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {

            var serie = enumerable.Select(func);
            var sxx = serie.Sum(p=>Math.Pow(p.x,2))-1.0/ serie.Count()* Math.Pow(serie.Sum(p => p.x), 2);
            var syy = serie.Sum(p => Math.Pow(p.y, 2)) - 1.0 / serie.Count() * Math.Pow(serie.Sum(p => p.y), 2);
            var sxy = serie.Sum(p => p.x*p.y) - 1.0 / serie.Count() * (serie.Sum(p => p.x) * serie.Sum(p => p.y));

            var slope = sxy / sxx;
            var intercept = serie.Average(p=>p.y)-slope * serie.Average(p => p.x);
            var correlation = sxy / Math.Sqrt(sxx * syy);

            return (slope, intercept, correlation);
        }


        /// <summary>
        /// Fits an exponential model y = exp(intercept) * exp(slope * x) by applying linear regression to log(y).
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to (x,y).</param>
        /// <returns>A function f(x) representing the fitted exponential curve.</returns>
        public static Func<double,double> ExponentialRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var (slope, intercept, correlation) = enumerable.Select(func).LinearRegression(p => (p.x, Math.Log(p.y)));

            return (double x)=> Math.Exp(intercept)* Math.Exp(slope*x);
        }



        /// <summary>
        /// Applies a logistic function to data points.
        /// Value = 1 / (1 + exp(-(slope * y + intercept))).
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to (x,y).</param>
        /// <param name="slope">Slope parameter.</param>
        /// <param name="intercept">Intercept parameter.</param>
        /// <returns>A sequence of <see cref="Serie"/> values after applying the logistic transform.</returns>
        public static IEnumerable<Serie> LogisticRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func, double slope, double intercept)
        {

            return enumerable.Select(func).Select(p =>new Serie() { Value = 1 / (1 + Math.Exp(-(slope * p.y + intercept))), Index= p.y});
  
        }

        /// <summary>
        /// Applies the rectified linear unit (ReLU) transform: relu(x) = max(0, x).
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The input sequence.</param>
        /// <param name="func">Projection from element to (x,y).</param>
        /// <returns>A sequence of <see cref="Serie"/> values after applying ReLU to x.</returns>
        public static IEnumerable<Serie> Rectifier<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {

            return enumerable.Select(func).Select(p => new Serie() { Value =p.x>0?p.x:0, Index = p.y });

        }
        /// <summary>
        /// Classifies an unknown point using the k-nearest neighbors (k-NN) rule.
        /// Selects the k closest points and returns the most frequent classification.
        /// </summary>
        /// <typeparam name="T">The element type.</typeparam>
        /// <param name="enumerable">The dataset.</param>
        /// <param name="func">Projection from element to (x, y, classification).</param>
        /// <param name="unkown">The unknown point to classify.</param>
        /// <param name="k">Number of neighbors.</param>
        /// <returns>The predicted classification label.</returns>
        public static int KnearestNeighbors<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y, int classification)> func, (double x, double y) unkown, int k)
        {
            var data = enumerable.Select(func);

            var pointInRing = new List<(double x, double y, int classification, double distance)>();

            foreach(var (x, y, classification) in data)
            {
                var distance = Math.Sqrt(Math.Pow(unkown.x - x, 2) - Math.Pow(unkown.y - y, 2));
                pointInRing.Add((x, y, classification, distance));     
            }

            return pointInRing.OrderBy(p => p.distance).Take(k).GroupBy(x => x.classification)
                          .OrderByDescending(x => x.Count())
                          .First().Key;

        }

        /// <summary>
        /// Generates a uniformly distributed random double in [min, max).
        /// </summary>
        /// <param name="rnd">Random number generator.</param>
        /// <param name="min">Lower bound.</param>
        /// <param name="max">Upper bound.</param>
        /// <returns>A random value in the specified range.</returns>
        public static double RandomDouble(this Random rnd, double min, double max)
        {
            return rnd.NextDouble() * (max - min) + min;
        }
    }
}