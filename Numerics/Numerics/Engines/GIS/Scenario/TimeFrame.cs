using System;

namespace CSharpNumerics.Engines.GIS.Scenario
{
    /// <summary>
    /// Defines a uniform time range for a simulation: start, end (seconds)
    /// and step size. Generates the discrete time values used by the simulator.
    /// </summary>
    public class TimeFrame
    {
        /// <summary>Start time in seconds.</summary>
        public double Start { get; }

        /// <summary>End time in seconds (inclusive).</summary>
        public double End { get; }

        /// <summary>Time step in seconds.</summary>
        public double StepSeconds { get; }

        /// <summary>Number of discrete time steps.</summary>
        public int Count { get; }

        /// <summary>
        /// Creates a time frame.
        /// </summary>
        /// <param name="start">Start time in seconds.</param>
        /// <param name="end">End time in seconds.</param>
        /// <param name="stepSeconds">Step size in seconds.</param>
        public TimeFrame(double start, double end, double stepSeconds)
        {
            if (stepSeconds <= 0) throw new ArgumentException("Step must be positive.", nameof(stepSeconds));
            if (end < start) throw new ArgumentException("End must be >= start.", nameof(end));

            Start = start;
            End = end;
            StepSeconds = stepSeconds;
            Count = Math.Max(1, (int)Math.Floor((end - start) / stepSeconds) + 1);
        }

        /// <summary>
        /// Returns the time value at the given zero-based index.
        /// </summary>
        public double TimeAt(int index)
        {
            if (index < 0 || index >= Count)
                throw new ArgumentOutOfRangeException(nameof(index));
            return Start + index * StepSeconds;
        }

        /// <summary>
        /// Returns all discrete time values as an array.
        /// </summary>
        public double[] ToArray()
        {
            var times = new double[Count];
            for (int i = 0; i < Count; i++)
                times[i] = Start + i * StepSeconds;
            return times;
        }

        /// <summary>
        /// Returns the zero-based index of the time step nearest to <paramref name="timeSeconds"/>.
        /// Clamps to [0, Count-1].
        /// </summary>
        public int NearestIndex(double timeSeconds)
        {
            int idx = (int)Math.Round((timeSeconds - Start) / StepSeconds);
            return idx < 0 ? 0 : idx >= Count ? Count - 1 : idx;
        }
    }
}
