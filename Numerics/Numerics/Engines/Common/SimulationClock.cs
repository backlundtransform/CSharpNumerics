using System;

namespace CSharpNumerics.Engines.Common
{
    /// <summary>
    /// Fixed-timestep clock with accumulator for framerate-independent simulation.
    /// </summary>
    public class SimulationClock
    {
        private double _accumulator;

        /// <summary>Fixed time step in seconds. Default 1/60.</summary>
        public double FixedDt { get; set; } = 1.0 / 60.0;

        /// <summary>Total elapsed simulation time in seconds.</summary>
        public double ElapsedTime { get; private set; }

        /// <summary>Number of fixed steps taken so far.</summary>
        public long TickCount { get; private set; }

        /// <summary>
        /// Feed a variable-length frame delta. Returns the number of fixed steps to execute.
        /// </summary>
        public int Accumulate(double frameDt)
        {
            if (frameDt < 0) throw new ArgumentOutOfRangeException(nameof(frameDt));
            _accumulator += frameDt;
            int steps = 0;
            while (_accumulator >= FixedDt)
            {
                _accumulator -= FixedDt;
                ElapsedTime += FixedDt;
                TickCount++;
                steps++;
            }
            return steps;
        }

        /// <summary>Fraction of a fixed step remaining in the accumulator (0-1), useful for interpolation.</summary>
        public double Alpha => _accumulator / FixedDt;

        /// <summary>Reset all counters to zero.</summary>
        public void Reset()
        {
            _accumulator = 0;
            ElapsedTime = 0;
            TickCount = 0;
        }
    }
}
