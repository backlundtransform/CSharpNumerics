using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Gravitation;

/// <summary>
/// Celestial-mechanics relations for gravitating two-body and hierarchical
/// systems: orbital speed (vis-viva), reduced mass, spheres of gravitational
/// influence (Hill sphere, Laplace sphere of influence), the Roche limit, and
/// tidal acceleration. SI units throughout (m, kg, s).
/// </summary>
public static class CelestialMechanics
{
    private const double G = PhysicsConstants.GravitationalConstant;

    /// <summary>
    /// Vis-viva orbital speed v = √(μ·(2/r − 1/a)), where μ is the standard
    /// gravitational parameter (G·M) of the central body.
    /// </summary>
    /// <param name="gravitationalParameter">μ = G·M of the central body (m³/s²).</param>
    /// <param name="radius">Current distance from the central body (m).</param>
    /// <param name="semiMajorAxis">Orbital semi-major axis a (m).</param>
    public static double VisVivaSpeed(double gravitationalParameter, double radius, double semiMajorAxis)
    {
        if (radius <= 0) throw new ArgumentOutOfRangeException(nameof(radius), "Radius must be positive.");
        if (semiMajorAxis == 0) throw new ArgumentOutOfRangeException(nameof(semiMajorAxis), "Semi-major axis must be non-zero.");

        double v2 = gravitationalParameter * (2.0 / radius - 1.0 / semiMajorAxis);
        if (v2 < 0) throw new ArgumentException("No real speed: r and a are inconsistent for a bound orbit.");
        return Math.Sqrt(v2);
    }

    /// <summary>
    /// Vis-viva speed from the central mass directly (μ = G·M).
    /// </summary>
    public static double VisVivaSpeedFromMass(double centralMass, double radius, double semiMajorAxis)
        => VisVivaSpeed(G * centralMass, radius, semiMajorAxis);

    /// <summary>
    /// Reduced mass of a two-body system: μ = m₁·m₂/(m₁ + m₂).
    /// </summary>
    public static double ReducedMass(double mass1, double mass2)
    {
        if (mass1 <= 0 || mass2 <= 0) throw new ArgumentOutOfRangeException(nameof(mass1), "Masses must be positive.");
        return mass1 * mass2 / (mass1 + mass2);
    }

    /// <summary>
    /// Hill-sphere radius of a small body orbiting a larger one:
    /// r_H ≈ a·(1 − e)·(m/(3M))^(1/3). Inside this radius the small body
    /// dominates the gravitational attraction of satellites.
    /// </summary>
    /// <param name="semiMajorAxis">Semi-major axis of the small body's orbit (m).</param>
    /// <param name="eccentricity">Eccentricity of that orbit (0 ≤ e &lt; 1).</param>
    /// <param name="smallMass">Mass of the orbiting body m (kg).</param>
    /// <param name="largeMass">Mass of the primary M (kg).</param>
    public static double HillSphereRadius(double semiMajorAxis, double eccentricity, double smallMass, double largeMass)
    {
        if (semiMajorAxis <= 0) throw new ArgumentOutOfRangeException(nameof(semiMajorAxis), "Semi-major axis must be positive.");
        if (eccentricity < 0 || eccentricity >= 1) throw new ArgumentOutOfRangeException(nameof(eccentricity), "Eccentricity must be in [0, 1).");

        return semiMajorAxis * (1.0 - eccentricity) * Math.Pow(smallMass / (3.0 * largeMass), 1.0 / 3.0);
    }

    /// <summary>
    /// Laplace sphere of influence radius: r_SOI ≈ a·(m/M)^(2/5). The boundary
    /// inside which the small body is treated as the dominant gravitational
    /// source for patched-conic trajectory design.
    /// </summary>
    /// <param name="semiMajorAxis">Semi-major axis of the small body's orbit (m).</param>
    /// <param name="smallMass">Mass of the orbiting body m (kg).</param>
    /// <param name="largeMass">Mass of the primary M (kg).</param>
    public static double SphereOfInfluence(double semiMajorAxis, double smallMass, double largeMass)
    {
        if (semiMajorAxis <= 0) throw new ArgumentOutOfRangeException(nameof(semiMajorAxis), "Semi-major axis must be positive.");
        return semiMajorAxis * Math.Pow(smallMass / largeMass, 2.0 / 5.0);
    }

    /// <summary>
    /// Roche limit for a rigid satellite: d = R·(2·ρ_M/ρ_m)^(1/3), the closest a
    /// rigid body can approach the primary before tidal forces overcome its own
    /// self-gravity.
    /// </summary>
    /// <param name="primaryRadius">Radius of the primary body R (m).</param>
    /// <param name="primaryDensity">Bulk density of the primary ρ_M (kg/m³).</param>
    /// <param name="satelliteDensity">Bulk density of the satellite ρ_m (kg/m³).</param>
    public static double RocheLimitRigid(double primaryRadius, double primaryDensity, double satelliteDensity)
    {
        ValidateRocheInputs(primaryRadius, primaryDensity, satelliteDensity);
        return primaryRadius * Math.Pow(2.0 * primaryDensity / satelliteDensity, 1.0 / 3.0);
    }

    /// <summary>
    /// Roche limit for a fluid (deformable) satellite: d ≈ 2.44·R·(ρ_M/ρ_m)^(1/3).
    /// </summary>
    public static double RocheLimitFluid(double primaryRadius, double primaryDensity, double satelliteDensity)
    {
        ValidateRocheInputs(primaryRadius, primaryDensity, satelliteDensity);
        return 2.44 * primaryRadius * Math.Pow(primaryDensity / satelliteDensity, 1.0 / 3.0);
    }

    /// <summary>
    /// Differential tidal acceleration across a body of radius
    /// <paramref name="bodyRadius"/> located a distance <paramref name="distance"/>
    /// from a mass <paramref name="mass"/>: a_tidal = 2·G·M·r/d³.
    /// </summary>
    /// <param name="mass">Mass raising the tide M (kg).</param>
    /// <param name="distance">Distance between the two bodies' centres d (m).</param>
    /// <param name="bodyRadius">Radius of the affected body r (m).</param>
    public static double TidalAcceleration(double mass, double distance, double bodyRadius)
    {
        if (distance <= 0) throw new ArgumentOutOfRangeException(nameof(distance), "Distance must be positive.");
        return 2.0 * G * mass * bodyRadius / (distance * distance * distance);
    }

    private static void ValidateRocheInputs(double primaryRadius, double primaryDensity, double satelliteDensity)
    {
        if (primaryRadius <= 0) throw new ArgumentOutOfRangeException(nameof(primaryRadius), "Radius must be positive.");
        if (primaryDensity <= 0) throw new ArgumentOutOfRangeException(nameof(primaryDensity), "Density must be positive.");
        if (satelliteDensity <= 0) throw new ArgumentOutOfRangeException(nameof(satelliteDensity), "Density must be positive.");
    }
}
