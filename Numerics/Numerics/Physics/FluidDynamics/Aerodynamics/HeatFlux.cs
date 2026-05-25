using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Aerodynamic heating model for rocket ascent and reentry.
/// Computes stagnation-point convective heat flux using the Sutton-Graves approximation:
///   q̇ = k · √(ρ/r_n) · V³
/// where:
///   k    = Sutton-Graves constant (depends on gas composition, ~1.7415e-4 for air/N₂)
///   ρ    = freestream density (kg/m³)
///   r_n  = nose radius (m)
///   V    = freestream velocity (m/s)
/// 
/// Reference: Sutton & Graves, "A General Stagnation-Point Convective-Heating Equation
/// for Arbitrary Gas Mixtures", NASA TR R-376, 1971.
/// </summary>
public static class HeatFlux
{
    /// <summary>Sutton-Graves constant for Earth atmosphere (air), W/(m² · (kg/m³)^0.5 · (m/s)^3).</summary>
    private const double K_Air = 1.7415e-4;

    /// <summary>
    /// Computes stagnation-point convective heat flux using Sutton-Graves.
    /// </summary>
    /// <param name="density">Freestream air density in kg/m³.</param>
    /// <param name="velocity">Freestream velocity in m/s.</param>
    /// <param name="noseRadius">Effective nose radius in metres.</param>
    /// <returns>Heat flux in W/m² (watts per square metre).</returns>
    public static double StagnationPoint(double density, double velocity, double noseRadius)
    {
        if (density <= 0 || velocity <= 0 || noseRadius <= 0) return 0;
        return K_Air * Math.Sqrt(density / noseRadius) * velocity * velocity * velocity;
    }

    /// <summary>
    /// Computes stagnation-point heat flux from altitude and velocity.
    /// Returns 0 above 86 km (edge of ISA atmosphere model).
    /// </summary>
    /// <param name="altitude">Altitude in metres.</param>
    /// <param name="velocity">Freestream velocity in m/s.</param>
    /// <param name="noseRadius">Effective nose radius in metres.</param>
    /// <returns>Heat flux in W/m².</returns>
    public static double FromAltitudeAndVelocity(double altitude, double velocity, double noseRadius)
    {
        if (altitude >= 86000) return 0;
        double rho = AtmosphereModel.Density(altitude);
        return StagnationPoint(rho, velocity, noseRadius);
    }

    /// <summary>
    /// Computes total stagnation enthalpy for reference:
    /// H = V²/2 (kinetic) + h_∞ (static enthalpy of freestream)
    /// At orbital velocities (~7.8 km/s), kinetic dominates.
    /// </summary>
    /// <param name="velocity">Freestream velocity in m/s.</param>
    /// <returns>Specific stagnation enthalpy in J/kg.</returns>
    public static double StagnationEnthalpy(double velocity)
    {
        return 0.5 * velocity * velocity;
    }
}
