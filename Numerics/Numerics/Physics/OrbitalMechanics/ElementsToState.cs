using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.OrbitalMechanics;

/// <summary>
/// Converts classical Keplerian orbital elements to an ECI state vector (position, velocity).
/// </summary>
public static class ElementsToState
{
    /// <summary>
    /// Converts orbital elements to position and velocity vectors in the ECI frame.
    /// </summary>
    /// <param name="elements">Classical Keplerian orbital elements.</param>
    /// <returns>Tuple of (position, velocity) vectors in ECI frame (meters, m/s).</returns>
    public static (Vector Position, Vector Velocity) ToStateVector(OrbitalElements elements)
    {
        double a = elements.SemiMajorAxis;
        double e = elements.Eccentricity;
        double i = elements.Inclination;
        double raan = elements.RAAN;
        double omega = elements.ArgumentOfPeriapsis;
        double nu = elements.TrueAnomaly;
        double mu = elements.Mu;

        // Semi-latus rectum
        double p = a * (1.0 - e * e);

        // Radius at true anomaly
        double r = p / (1.0 + e * Math.Cos(nu));

        // Position in perifocal (PQW) frame
        double xPQW = r * Math.Cos(nu);
        double yPQW = r * Math.Sin(nu);

        // Velocity in perifocal frame
        double sqrtMuP = Math.Sqrt(mu / p);
        double vxPQW = -sqrtMuP * Math.Sin(nu);
        double vyPQW = sqrtMuP * (e + Math.Cos(nu));

        // Rotation matrix elements (perifocal → ECI)
        double cosO = Math.Cos(raan);
        double sinO = Math.Sin(raan);
        double cosI = Math.Cos(i);
        double sinI = Math.Sin(i);
        double cosW = Math.Cos(omega);
        double sinW = Math.Sin(omega);

        // First column of rotation matrix
        double r11 = cosO * cosW - sinO * sinW * cosI;
        double r21 = sinO * cosW + cosO * sinW * cosI;
        double r31 = sinW * sinI;

        // Second column of rotation matrix
        double r12 = -cosO * sinW - sinO * cosW * cosI;
        double r22 = -sinO * sinW + cosO * cosW * cosI;
        double r32 = cosW * sinI;

        // Transform position to ECI
        double x = r11 * xPQW + r12 * yPQW;
        double y = r21 * xPQW + r22 * yPQW;
        double z = r31 * xPQW + r32 * yPQW;

        // Transform velocity to ECI
        double vx = r11 * vxPQW + r12 * vyPQW;
        double vy = r21 * vxPQW + r22 * vyPQW;
        double vz = r31 * vxPQW + r32 * vyPQW;

        return (new Vector(x, y, z), new Vector(vx, vy, vz));
    }
}
