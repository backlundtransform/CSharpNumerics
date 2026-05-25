using System;

namespace CSharpNumerics.Physics.Mechanics;

/// <summary>
/// WGS84 Earth ellipsoid parameters and derived constants.
/// </summary>
public static class EarthModel
{
    /// <summary>Semi-major axis (equatorial radius) in meters.</summary>
    public const double SemiMajorAxis = 6378137.0;

    /// <summary>Flattening factor.</summary>
    public const double Flattening = 1.0 / 298.257223563;

    /// <summary>Semi-minor axis (polar radius) in meters.</summary>
    public static readonly double SemiMinorAxis = SemiMajorAxis * (1.0 - Flattening);

    /// <summary>First eccentricity squared: e² = 2f - f².</summary>
    public static readonly double EccentricitySquared = 2.0 * Flattening - Flattening * Flattening;

    /// <summary>Second eccentricity squared: e'² = (a²-b²)/b².</summary>
    public static readonly double SecondEccentricitySquared =
        (SemiMajorAxis * SemiMajorAxis - SemiMinorAxis * SemiMinorAxis) / (SemiMinorAxis * SemiMinorAxis);

    /// <summary>Earth's rotation rate in rad/s.</summary>
    public const double RotationRate = 7.2921150e-5;

    /// <summary>Standard gravitational parameter GM (m³/s²).</summary>
    public const double GM = 3.986004418e14;

    /// <summary>J2 zonal harmonic coefficient (oblateness).</summary>
    public const double J2 = 1.08263e-3;

    /// <summary>Mean radius in meters (volumetric).</summary>
    public const double MeanRadius = 6371000.0;

    /// <summary>Standard gravitational acceleration at sea level (m/s²).</summary>
    public const double G0 = 9.80665;

    /// <summary>
    /// Radius of curvature in the prime vertical at a given geodetic latitude.
    /// N(φ) = a / sqrt(1 - e²·sin²φ)
    /// </summary>
    public static double PrimeVerticalRadius(double latitude)
    {
        double sinLat = Math.Sin(latitude);
        return SemiMajorAxis / Math.Sqrt(1.0 - EccentricitySquared * sinLat * sinLat);
    }

    /// <summary>
    /// Surface velocity due to Earth's rotation at a given latitude (m/s).
    /// v = ω · a · cos(φ) (approximate, ignoring ellipsoid correction).
    /// </summary>
    public static double SurfaceVelocity(double latitude)
    {
        double N = PrimeVerticalRadius(latitude);
        return RotationRate * (N + 0.0) * Math.Cos(latitude);
    }
}
