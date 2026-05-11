using System;

namespace CSharpNumerics.Physics.FluidDynamics.Buoyancy;

/// <summary>
/// Computes buoyancy forces for thermally-driven flows (smoke, fire, hot gas plumes).
/// Uses the Boussinesq approximation: density variations are neglected except in the
/// gravity term, where they produce a buoyancy acceleration:
///   a_buoy = −g · β · (T − T_ref)
/// 
/// For game smoke/fire simulation, temperature can be replaced by a generic "density"
/// field (hot = low density = rises, cold = high density = sinks).
/// </summary>
public static class BuoyancyForce
{
    /// <summary>
    /// Computes the vertical buoyancy acceleration using the Boussinesq approximation.
    /// Positive result means upward (against gravity).
    /// </summary>
    /// <param name="temperature">Local temperature (or smoke density proxy).</param>
    /// <param name="referenceTemperature">Ambient / reference temperature.</param>
    /// <param name="thermalExpansion">Volumetric thermal expansion coefficient β (1/K). Default 1/300 for air near 300 K.</param>
    /// <param name="gravity">Gravitational acceleration magnitude (m/s²). Default 9.80665.</param>
    /// <returns>Buoyancy acceleration in m/s² (positive = upward).</returns>
    public static double BoussinesqAcceleration(
        double temperature,
        double referenceTemperature,
        double thermalExpansion = 1.0 / 300.0,
        double gravity = 9.80665)
    {
        return gravity * thermalExpansion * (temperature - referenceTemperature);
    }

    /// <summary>
    /// Applies buoyancy as a velocity increment to the vertical component of a velocity field cell.
    /// v_z += g · β · (T − T_ref) · dt
    /// </summary>
    /// <param name="velocityZ">Current vertical velocity component (updated in-place by caller).</param>
    /// <param name="temperature">Local temperature or density proxy.</param>
    /// <param name="referenceTemperature">Ambient temperature.</param>
    /// <param name="dt">Time step in seconds.</param>
    /// <param name="thermalExpansion">Thermal expansion coefficient β.</param>
    /// <param name="gravity">Gravitational acceleration magnitude.</param>
    /// <returns>New vertical velocity after buoyancy impulse.</returns>
    public static double ApplyBuoyancy(
        double velocityZ,
        double temperature,
        double referenceTemperature,
        double dt,
        double thermalExpansion = 1.0 / 300.0,
        double gravity = 9.80665)
    {
        return velocityZ + BoussinesqAcceleration(temperature, referenceTemperature, thermalExpansion, gravity) * dt;
    }

    /// <summary>
    /// Computes buoyancy from a density field (lighter fluid rises).
    /// a_buoy = g · (ρ_ref − ρ) / ρ_ref
    /// </summary>
    /// <param name="density">Local density.</param>
    /// <param name="referenceDensity">Ambient density.</param>
    /// <param name="gravity">Gravitational acceleration magnitude.</param>
    /// <returns>Buoyancy acceleration (positive = upward).</returns>
    public static double DensityBuoyancy(
        double density,
        double referenceDensity,
        double gravity = 9.80665)
    {
        if (referenceDensity <= 0) return 0;
        return gravity * (referenceDensity - density) / referenceDensity;
    }
}
