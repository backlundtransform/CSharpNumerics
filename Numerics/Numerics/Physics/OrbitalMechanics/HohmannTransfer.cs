using System;
using CSharpNumerics.Physics.Mechanics;

namespace CSharpNumerics.Physics.OrbitalMechanics;

/// <summary>
/// Computes Hohmann transfer orbit parameters between two circular orbits.
/// A Hohmann transfer is the minimum-energy two-impulse transfer between coplanar circular orbits.
/// </summary>
public static class HohmannTransfer
{
    /// <summary>
    /// Computes the Hohmann transfer parameters between two circular orbit radii.
    /// </summary>
    /// <param name="r1">Radius of initial circular orbit (meters).</param>
    /// <param name="r2">Radius of target circular orbit (meters).</param>
    /// <param name="mu">Gravitational parameter (m³/s²). Default: Earth.</param>
    /// <returns>Transfer parameters.</returns>
    public static HohmannResult Compute(double r1, double r2, double mu = 0)
    {
        if (mu == 0) mu = EarthModel.GM;

        // Circular velocities
        double v1 = Math.Sqrt(mu / r1);
        double v2 = Math.Sqrt(mu / r2);

        // Transfer orbit semi-major axis
        double aTransfer = (r1 + r2) / 2.0;

        // Velocities on transfer orbit at departure and arrival
        double vTransfer1 = Math.Sqrt(mu * (2.0 / r1 - 1.0 / aTransfer));
        double vTransfer2 = Math.Sqrt(mu * (2.0 / r2 - 1.0 / aTransfer));

        // Delta-V for each burn
        double dv1 = vTransfer1 - v1; // Positive for raising orbit
        double dv2 = v2 - vTransfer2; // Positive for raising orbit

        // Transfer time (half of transfer orbit period)
        double transferTime = Math.PI * Math.Sqrt(aTransfer * aTransfer * aTransfer / mu);

        return new HohmannResult
        {
            DeltaV1 = dv1,
            DeltaV2 = dv2,
            TotalDeltaV = Math.Abs(dv1) + Math.Abs(dv2),
            TransferTime = transferTime,
            TransferSemiMajorAxis = aTransfer,
            InitialVelocity = v1,
            FinalVelocity = v2
        };
    }

    /// <summary>
    /// Computes Hohmann transfer between two circular orbits defined by altitude above Earth's surface.
    /// </summary>
    /// <param name="altitude1">Initial orbit altitude (meters above surface).</param>
    /// <param name="altitude2">Target orbit altitude (meters above surface).</param>
    /// <returns>Transfer parameters.</returns>
    public static HohmannResult FromAltitudes(double altitude1, double altitude2)
    {
        double r1 = EarthModel.SemiMajorAxis + altitude1;
        double r2 = EarthModel.SemiMajorAxis + altitude2;
        return Compute(r1, r2, EarthModel.GM);
    }
}

/// <summary>
/// Result of a Hohmann transfer computation.
/// </summary>
public struct HohmannResult
{
    /// <summary>ΔV for the first burn (departure). Positive = prograde.</summary>
    public double DeltaV1;

    /// <summary>ΔV for the second burn (arrival/circularization). Positive = prograde.</summary>
    public double DeltaV2;

    /// <summary>Total |ΔV| for the transfer.</summary>
    public double TotalDeltaV;

    /// <summary>Transfer time in seconds (half the transfer orbit period).</summary>
    public double TransferTime;

    /// <summary>Semi-major axis of the transfer ellipse (meters).</summary>
    public double TransferSemiMajorAxis;

    /// <summary>Initial circular orbit velocity (m/s).</summary>
    public double InitialVelocity;

    /// <summary>Final circular orbit velocity (m/s).</summary>
    public double FinalVelocity;
}
