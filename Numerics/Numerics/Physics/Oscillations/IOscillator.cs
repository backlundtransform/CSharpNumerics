using CSharpNumerics.Statistics.Data;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Oscillations
{
    /// <summary>
    /// Common interface for one-dimensional oscillators.
    /// Position and Velocity describe the instantaneous state;
    /// <see cref="Step"/> advances the simulation by one time step.
    /// </summary>
    public interface IOscillator
    {
        /// <summary>Current displacement from equilibrium (m).</summary>
        double Position { get; }

        /// <summary>Current velocity (m/s).</summary>
        double Velocity { get; }

        /// <summary>Current simulation time (s).</summary>
        double Time { get; }

        /// <summary>Total mechanical energy at the current state (J).</summary>
        double TotalEnergy { get; }

        /// <summary>Kinetic energy at the current state (J).</summary>
        double KineticEnergy { get; }

        /// <summary>Potential energy at the current state (J).</summary>
        double PotentialEnergy { get; }

        /// <summary>
        /// Advances the oscillator by one time step.
        /// </summary>
        /// <param name="dt">Time step in seconds.</param>
        void Step(double dt);

        /// <summary>
        /// Resets the oscillator to its initial conditions at t = 0.
        /// </summary>
        void Reset();

        /// <summary>
        /// Simulates the oscillator from t = 0 to <paramref name="tEnd"/> and returns the
        /// displacement as a time series.
        /// The oscillator state is restored to its initial conditions after the call.
        /// </summary>
        /// <param name="tEnd">End time in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>Displacement vs time.</returns>
        List<Serie> Trajectory(double tEnd, double dt);

        /// <summary>
        /// Returns the phase-space trajectory (x, v) from t = 0 to <paramref name="tEnd"/>.
        /// The oscillator state is restored to its initial conditions after the call.
        /// </summary>
        /// <param name="tEnd">End time in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>List of (position, velocity) data points.</returns>
        List<Serie> PhasePortrait(double tEnd, double dt);
    }
}
