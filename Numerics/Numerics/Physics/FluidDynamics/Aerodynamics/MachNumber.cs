using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Mach number computation and classification for compressible flow.
/// M = V / a(h), where a = speed of sound from <see cref="AtmosphereModel"/>.
/// </summary>
public static class MachNumber
{
    /// <summary>
    /// Computes Mach number from velocity and altitude.
    /// </summary>
    /// <param name="velocity">True airspeed in m/s.</param>
    /// <param name="altitude">Altitude in metres.</param>
    /// <returns>Mach number (dimensionless).</returns>
    public static double Compute(double velocity, double altitude)
    {
        double a = AtmosphereModel.SpeedOfSound(altitude);
        return a > 0 ? Math.Abs(velocity) / a : 0;
    }

    /// <summary>True if subsonic (M &lt; 0.8).</summary>
    public static bool IsSubsonic(double mach) => mach < 0.8;

    /// <summary>True if transonic (0.8 ≤ M ≤ 1.2).</summary>
    public static bool IsTransonic(double mach) => mach >= 0.8 && mach <= 1.2;

    /// <summary>True if supersonic (1.2 &lt; M ≤ 5.0).</summary>
    public static bool IsSupersonic(double mach) => mach > 1.2 && mach <= 5.0;

    /// <summary>True if hypersonic (M &gt; 5.0).</summary>
    public static bool IsHypersonic(double mach) => mach > 5.0;
}
