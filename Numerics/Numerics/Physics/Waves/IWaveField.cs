using CSharpNumerics.Statistics.Data;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Waves
{
    /// <summary>
    /// Common contract for wave-field simulations.
    /// Mirrors <see cref="CSharpNumerics.Physics.Oscillations.IOscillator"/> for spatial fields.
    /// </summary>
    public interface IWaveField
    {
        /// <summary>Current simulation time (s).</summary>
        double Time { get; }

        /// <summary>Total mechanical energy of the wave field (J).</summary>
        double TotalEnergy { get; }

        /// <summary>
        /// Advances the wave field by one time step.
        /// </summary>
        /// <param name="dt">Time step in seconds.</param>
        void Step(double dt);

        /// <summary>
        /// Resets the field to its initial conditions at t = 0.
        /// </summary>
        void Reset();

        /// <summary>
        /// Returns the current displacement field as a list of (position, amplitude) data points.
        /// </summary>
        List<Serie> Snapshot();
    }
}
