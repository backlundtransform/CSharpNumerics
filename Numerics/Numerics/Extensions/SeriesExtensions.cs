using CSharpNumerics.Enums;
using Numerics.Models;
using System.Collections.Generic;
using System.Linq;

namespace System
{
    public static class SeriesExtensions
    {
        public static List<Serie> GetSeries(this Func<double, double> func, double minValue, double maxValue, double stepSize)
        {
            var numbers = new List<Serie>();

            var step = (maxValue - minValue) / stepSize;

            for (var i = minValue; i < maxValue; i += step)
            {
                numbers.Add(new Serie() { Value = func(i), Index = i });

            }

            return numbers;

        }

        public static List<TimeSerie> GetTimeSeries(this Func<DateTime, double> func, DateTime startTime, DateTime endTime, double stepSize)
        {
            var timeSeries = new List<TimeSerie>();

            var step = (startTime - endTime).TotalMilliseconds / stepSize;

            for (var i = startTime; i < endTime; i.AddMilliseconds(step))
            {
                timeSeries.Add(new TimeSerie() { Value = func(i), TimeStamp = i });

            }

            return timeSeries;

        }
        public static List<TimeSerie> GenerateTimeSerieWithEquivalentSteps(this List<TimeSerie> timeSeries, int minutes, DateTime startDate, DateTime endDate, GroupOperator groupOperator= GroupOperator.Average, int multiplier=1, bool shouldInterpolate=true)
        {
            var groups = timeSeries.GroupBy(x =>
            {
                var stamp = x.TimeStamp;
                stamp = stamp.AddMinutes(-(stamp.Minute % minutes));
                stamp = stamp.AddMilliseconds(-stamp.Millisecond - 1000 * stamp.Second);
                return stamp;
            }).ToDictionary(g => g.Key, g => g.ToList());

            var newTimeserie = new List<TimeSerie>();

            for (var date = startDate; date <= endDate; date = date.AddMinutes(minutes))
            {
                var value = 0.0;

                if (!groups.ContainsKey(date) && shouldInterpolate)
                {
                    value = timeSeries.LinearInterpolationTimeSerie(date);
                }

                if (groups.ContainsKey(date))
                {

                    switch (groupOperator)
                    {
                        case GroupOperator.Average:
                            value =groups[date].Average(p => p.Value);
                            break;
                        case GroupOperator.Max:
                            value = groups[date].Max(p => p.Value);
                            break;
                        case GroupOperator.Min:
                            value = groups[date].Min(p => p.Value);
                            break;
                        case GroupOperator.Sum:
                            value = groups[date].Sum(p => p.Value);
                            break;
                        case GroupOperator.Median:
                            value = groups[date].Median(p => p.Value);
                            break;
                    }
                }

                   newTimeserie.Add(new TimeSerie { TimeStamp = date, Value = (multiplier) * value });
                
            }
            return newTimeserie;
        }
    }
}