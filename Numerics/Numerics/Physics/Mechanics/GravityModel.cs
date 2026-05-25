using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Mechanics;

/// <summary>
/// Computes gravitational acceleration including J2 oblateness perturbation.
/// Uses the WGS84 Earth model for all parameters.
/// </summary>
public static class GravityModel
{
    /// <summary>
    /// Computes gravitational acceleration at a given ECEF or ECI position vector.
    /// Includes J2 zonal harmonic for oblateness effects.
    /// Returns acceleration vector in the same frame as the input position.
    /// </summary>
    /// <param name="position">Position vector from Earth's center (meters).</param>
    /// <returns>Gravitational acceleration vector (m/s²), pointing toward Earth.</returns>
    public static Vector Acceleration(Vector position)
    {
        double x = position.x;
        double y = position.y;
        double z = position.z;

        double r2 = x * x + y * y + z * z;
        double r = Math.Sqrt(r2);
        double r5 = r2 * r2 * r;

        double mu = EarthModel.GM;
        double j2 = EarthModel.J2;
        double re = EarthModel.SemiMajorAxis;
        double re2 = re * re;

        double z2_r2 = z * z / r2;
        double coeff = -mu / (r2 * r); // -μ/r³
        double j2Coeff = 1.5 * j2 * re2 / r2;

        // J2 perturbed acceleration components
        double ax = coeff * x * (1.0 + j2Coeff * (1.0 - 5.0 * z2_r2));
        double ay = coeff * y * (1.0 + j2Coeff * (1.0 - 5.0 * z2_r2));
        double az = coeff * z * (1.0 + j2Coeff * (3.0 - 5.0 * z2_r2));

        return new Vector(ax, ay, az);
    }

    /// <summary>
    /// Computes gravitational acceleration from an ECIPosition.
    /// </summary>
    public static Vector Acceleration(ECIPosition position)
    {
        return Acceleration(new Vector(position.X, position.Y, position.Z));
    }

    /// <summary>
    /// Computes gravitational acceleration from an ECEFPosition.
    /// </summary>
    public static Vector Acceleration(ECEFPosition position)
    {
        return Acceleration(new Vector(position.X, position.Y, position.Z));
    }

    /// <summary>
    /// Computes the magnitude of gravitational acceleration at a given radius and latitude.
    /// g(r,φ) ≈ μ/r² · [1 + 1.5·J2·(Re/r)²·(1 - 3sin²φ)]
    /// Useful for scalar gravity comparisons.
    /// </summary>
    /// <param name="radius">Distance from Earth's center (meters).</param>
    /// <param name="latitude">Geodetic latitude (radians).</param>
    public static double Magnitude(double radius, double latitude)
    {
        double r2 = radius * radius;
        double re_r = EarthModel.SemiMajorAxis / radius;
        double sinLat = Math.Sin(latitude);

        return EarthModel.GM / r2 * (1.0 + 1.5 * EarthModel.J2 * re_r * re_r * (1.0 - 3.0 * sinLat * sinLat));
    }

    /// <summary>
    /// Simple point-mass gravity (no J2). g = -μ/r³ · r_vec
    /// </summary>
    public static Vector PointMassAcceleration(Vector position)
    {
        double r2 = position.x * position.x + position.y * position.y + position.z * position.z;
        double r = Math.Sqrt(r2);
        double coeff = -EarthModel.GM / (r2 * r);

        return new Vector(coeff * position.x, coeff * position.y, coeff * position.z);
    }
}
