using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.OrbitalMechanics;

/// <summary>
/// Converts between state vectors (position, velocity) and classical Keplerian orbital elements.
/// All methods assume an inertial (ECI) reference frame.
/// </summary>
public static class StateToElements
{
    /// <summary>
    /// Converts an ECI state vector (position, velocity) to classical orbital elements.
    /// </summary>
    /// <param name="position">Position vector in ECI frame (meters).</param>
    /// <param name="velocity">Velocity vector in ECI frame (m/s).</param>
    /// <param name="mu">Gravitational parameter μ = GM (m³/s²).</param>
    /// <returns>Classical Keplerian orbital elements.</returns>
    public static OrbitalElements FromStateVector(Vector position, Vector velocity, double mu)
    {
        double r = position.GetMagnitude();
        double v = velocity.GetMagnitude();

        // Specific angular momentum h = r × v
        Vector h = Cross(position, velocity);
        double hMag = h.GetMagnitude();

        // Node vector n = K × h (K = [0,0,1])
        Vector n = new Vector(-h.y, h.x, 0);
        double nMag = n.GetMagnitude();

        // Eccentricity vector e = (v×h)/μ - r̂
        Vector vCrossH = Cross(velocity, h);
        Vector eVec = new Vector(
            vCrossH.x / mu - position.x / r,
            vCrossH.y / mu - position.y / r,
            vCrossH.z / mu - position.z / r);
        double e = eVec.GetMagnitude();

        // Specific orbital energy
        double energy = 0.5 * v * v - mu / r;

        // Semi-major axis
        double a;
        if (Math.Abs(1.0 - e) > 1e-10)
            a = -mu / (2.0 * energy);
        else
            a = hMag * hMag / mu; // Parabolic: use semi-latus rectum

        // Inclination
        double i = Math.Acos(Clamp(h.z / hMag, -1.0, 1.0));

        // RAAN (Ω)
        double raan;
        if (nMag > 1e-10)
            raan = Math.Acos(Clamp(n.x / nMag, -1.0, 1.0));
        else
            raan = 0; // Equatorial orbit

        if (n.y < 0) raan = 2.0 * Math.PI - raan;

        // Argument of periapsis (ω)
        double omega;
        if (nMag > 1e-10 && e > 1e-10)
        {
            omega = Math.Acos(Clamp(Dot(n, eVec) / (nMag * e), -1.0, 1.0));
            if (eVec.z < 0) omega = 2.0 * Math.PI - omega;
        }
        else if (e > 1e-10)
        {
            // Equatorial orbit: measure from X-axis
            omega = Math.Atan2(eVec.y, eVec.x);
            if (omega < 0) omega += 2.0 * Math.PI;
        }
        else
        {
            omega = 0; // Circular orbit
        }

        // True anomaly (ν)
        double nu;
        if (e > 1e-10)
        {
            nu = Math.Acos(Clamp(Dot(eVec, position) / (e * r), -1.0, 1.0));
            if (Dot(position, velocity) < 0) nu = 2.0 * Math.PI - nu;
        }
        else if (nMag > 1e-10)
        {
            // Circular non-equatorial: measure from ascending node
            nu = Math.Acos(Clamp(Dot(n, position) / (nMag * r), -1.0, 1.0));
            if (position.z < 0) nu = 2.0 * Math.PI - nu;
        }
        else
        {
            // Circular equatorial: measure from X-axis
            nu = Math.Atan2(position.y, position.x);
            if (nu < 0) nu += 2.0 * Math.PI;
        }

        return new OrbitalElements(a, e, i, raan, omega, nu, mu);
    }

    private static Vector Cross(Vector a, Vector b)
    {
        return new Vector(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x);
    }

    private static double Dot(Vector a, Vector b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    private static double Clamp(double value, double min, double max)
    {
        if (value < min) return min;
        if (value > max) return max;
        return value;
    }
}
