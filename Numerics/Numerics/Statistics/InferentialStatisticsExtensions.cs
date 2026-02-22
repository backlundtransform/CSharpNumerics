using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Text;

namespace System.Linq

{
    public static class InferentialStatisticsExtensions
    {
        public static (double slope, double intercept, double correlation) LinearRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {

            var serie = enumerable.Select(func);
            var sxx = serie.Sum(p => Math.Pow(p.x, 2)) - 1.0 / serie.Count() * Math.Pow(serie.Sum(p => p.x), 2);
            var syy = serie.Sum(p => Math.Pow(p.y, 2)) - 1.0 / serie.Count() * Math.Pow(serie.Sum(p => p.y), 2);
            var sxy = serie.Sum(p => p.x * p.y) - 1.0 / serie.Count() * (serie.Sum(p => p.x) * serie.Sum(p => p.y));

            var slope = sxy / sxx;
            var intercept = serie.Average(p => p.y) - slope * serie.Average(p => p.x);
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
        public static Func<double, double> ExponentialRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var (slope, intercept, correlation) = enumerable.Select(func).LinearRegression(p => (p.x, Math.Log(p.y)));

            return (double x) => Math.Exp(intercept) * Math.Exp(slope * x);
        }


        }
}
