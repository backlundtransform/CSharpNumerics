using Numerics.Models;
using System.Collections.Generic;

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

    }
}
