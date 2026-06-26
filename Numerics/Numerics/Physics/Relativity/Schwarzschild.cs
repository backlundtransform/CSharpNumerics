using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Relativity;

/// <summary>
/// Geometry of the Schwarzschild solution — the spacetime around a static,
/// uncharged, spherically symmetric mass. Provides the characteristic radii
/// and the gravitational time-dilation / redshift experienced by a static
/// observer outside the mass. Masses are in kg, radii in metres.
/// </summary>
public static class Schwarzschild
{
    private const double G = PhysicsConstants.GravitationalConstant;
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double C2 = C * C;

    /// <summary>
    /// Schwarzschild (gravitational) radius r_s = 2GM/c². The event horizon of
    /// a non-rotating black hole of this mass.
    /// </summary>
    public static double Radius(double mass)
    {
        if (mass < 0) throw new ArgumentOutOfRangeException(nameof(mass), "Mass must be non-negative.");
        return 2.0 * G * mass / C2;
    }

    /// <summary>
    /// Photon-sphere radius r = 3GM/c² = 1.5·r_s, where light can orbit on
    /// (unstable) circular null geodesics.
    /// </summary>
    public static double PhotonSphereRadius(double mass)
    {
        return 1.5 * Radius(mass);
    }

    /// <summary>
    /// Innermost stable circular orbit (ISCO) for a massive test particle,
    /// r = 6GM/c² = 3·r_s. Inside this radius no stable circular orbit exists.
    /// </summary>
    public static double Isco(double mass)
    {
        return 3.0 * Radius(mass);
    }

    /// <summary>
    /// Gravitational time-dilation factor dτ/dt = √(1 − r_s/r) for a static
    /// observer at <paramref name="radius"/>: their proper time runs slow
    /// relative to a clock at infinity.
    /// </summary>
    /// <param name="radius">Areal radius (m), must be greater than r_s.</param>
    /// <param name="mass">Central mass (kg).</param>
    /// <exception cref="ArgumentOutOfRangeException">If r ≤ r_s.</exception>
    public static double TimeDilationFactor(double radius, double mass)
    {
        double rs = Radius(mass);
        if (radius <= rs)
            throw new ArgumentOutOfRangeException(nameof(radius),
                "Radius must be outside the Schwarzschild radius.");

        return Math.Sqrt(1.0 - rs / radius);
    }

    /// <summary>
    /// Gravitational redshift z = 1/√(1 − r_s/r) − 1 of light emitted by a
    /// static source at <paramref name="radius"/> and received at infinity.
    /// </summary>
    public static double Redshift(double radius, double mass)
    {
        return 1.0 / TimeDilationFactor(radius, mass) - 1.0;
    }

    /// <summary>
    /// Ratio of clock rates between two static radii, dτ(r1)/dτ(r2)
    /// = √((1 − r_s/r1)/(1 − r_s/r2)). Values &lt; 1 mean the clock at
    /// <paramref name="radiusLower"/> runs slow relative to the upper one.
    /// </summary>
    public static double ClockRateRatio(double radiusLower, double radiusUpper, double mass)
    {
        return TimeDilationFactor(radiusLower, mass) / TimeDilationFactor(radiusUpper, mass);
    }

    /// <summary>
    /// Newtonian escape velocity √(2GM/r). In the Schwarzschild geometry this
    /// equals c exactly at the event horizon (r = r_s).
    /// </summary>
    public static double EscapeVelocity(double radius, double mass)
    {
        if (radius <= 0) throw new ArgumentOutOfRangeException(nameof(radius), "Radius must be positive.");
        return Math.Sqrt(2.0 * G * mass / radius);
    }
}
