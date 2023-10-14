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
        public static List<TimeSerie> GenerateTimeSerieWithEquivalentSteps(this List<TimeSerie> timeSeries, 
            int minutes,
            DateTime startDate, 
            DateTime endDate,  
            int multiplier=1)
        {
            var resampledData = new List<TimeSerie>();
            var endTime = endDate;
            var currentTime = startDate;

            while (currentTime <= endTime)
            {
                var lowerFirstValue = timeSeries.LastOrDefault(p => p.TimeStamp <= currentTime);
                var upperFirstValue = timeSeries.FirstOrDefault(p => p.TimeStamp >= currentTime);

                if (lowerFirstValue == null || upperFirstValue == null)
                {
                    currentTime = currentTime.AddMinutes(minutes);
                    continue;
                }

                var currentValue = currentTime.Ticks.Interpolation(lowerFirstValue.TimeStamp.Ticks,
                    upperFirstValue.TimeStamp.Ticks,
                    lowerFirstValue.Value,
                    upperFirstValue.Value);

                resampledData.Add(new TimeSerie { TimeStamp = currentTime, Value = multiplier *currentValue });
             
                currentTime = currentTime.AddMinutes(minutes);
            }

            return resampledData;
      
        }

        public static double Interpolation(this long currentTicks,
            long firstTicks, 
            long lastTicks, 
            double firstValue,
            double lastValue)
        {
            if ((lastTicks - firstTicks) == 0)
            {
                return firstValue;
            }

            return firstValue + (lastValue - firstValue) * (currentTicks - firstTicks) / (lastTicks - firstTicks);
        }
    }
}