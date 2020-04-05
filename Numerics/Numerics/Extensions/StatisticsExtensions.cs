using Numerics.Methods;
using Numerics.Models;
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


        public static double GenerateNoise(this Random random, double variance)
        {
            var normalDistribution = Statistics.NormalDistribution(Math.Sqrt(variance), 0);

            return normalDistribution(random.Next((int)variance));
        }

        public static double GetRandomDouble(this Random random, double minimum, double maximum)
        {
            return random.NextDouble() * (maximum - minimum) + minimum;
        }


        public static double LinearInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp) {

            var prev = ts.FirstOrDefault(p => p.TimeStamp < timeStamp);
            var next = ts.LastOrDefault(p => p.TimeStamp > timeStamp);

            if(prev== null || next == null)
            {
                return 0;
            }

            var previousValue= prev.Value;
            var nextValue = next.Value;

            var previousTicks = prev.TimeStamp.Ticks;
            var nextTicks = next.TimeStamp.Ticks;

            var currentTicks = timeStamp.Ticks;

            return previousValue + (nextValue - previousValue) * (currentTicks - previousTicks) / (nextTicks - previousTicks);
        }


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


        public static Func<double,double> ExponentialRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var (slope, intercept, correlation) = enumerable.Select(func).LinearRegression(p => (p.x, Math.Log(p.y)));

            return (double x)=> Math.Exp(intercept)* Math.Exp(slope*x);
        }

    }
}