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
            var xy = enumerable.Select(func).ToList();

            var meanX = xy.Average(p => p.x);
            var meanY = xy.Average(p => p.y);

            return xy.Sum(p => (p.x - meanX) * (p.y - meanY)) / (xy.Count() - 1);
        }

        public static double CoefficientOfDetermination<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {
            var xy = enumerable.Select(func);
            var cov = xy.Covariance((p)=> (p.x, p.y));

            var stdDevX= xy.Select(p => p.x).StandardDeviation();
            var stdDevY= xy.Select(p => p.y).StandardDeviation();


            if (stdDevX == 0 || stdDevY == 0)
            {
                return 0.0;
            }


            double correlation = cov / (stdDevY * stdDevX);


 

            return correlation * correlation;
        }


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

            return ts.LinearInterpolation(p => (p.TimeStamp.Ticks, p.Value), timeStamp.Ticks);
        }


        public static double LinearInterpolation<T>(this IEnumerable<T> ts,  Func<T, (double x, double y)> func, double index)
        {
            var prevs = ts.Select(func);
            var prev = ts.Select(func).LastOrDefault(p => p.x < index);
            var next = ts.Select(func).FirstOrDefault(p => p.x > index);
            var previousValue = prev.y;
            var nextValue = next.y;

            var previousIndex = prev.x;
            var nextIndex = next.x;

            var currentIndex = index;

            return previousValue + (nextValue - previousValue) * (currentIndex - previousIndex) / (nextIndex - previousIndex);
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



        public static IEnumerable<Serie> LogisticRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func, double slope, double intercept)
        {

            return enumerable.Select(func).Select(p =>new Serie() { Value = 1 / (1 + Math.Exp(-(slope * p.y + intercept))), Index= p.y});
  
        }

        public static IEnumerable<Serie> Rectifier<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)
        {

            return enumerable.Select(func).Select(p => new Serie() { Value =p.x>0?p.x:0, Index = p.y });

        }
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

        public static double RandomDouble(this Random rnd, double min, double max)
        {
            return rnd.NextDouble() * (max - min) + min;
        }

    }
}