using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Mechanics;

/// <summary>
/// Computes the Coriolis acceleration for objects moving in the Earth-rotating (ECEF) frame.
/// a_coriolis = -2(ω × v) where ω is Earth's rotation vector.
/// </summary>
public static class CoriolisForce
{
    /// <summary>
    /// Computes the Coriolis acceleration given velocity in the ECEF frame.
    /// a = -2(ω × v), where ω = (0, 0, ω_earth).
    /// </summary>
    /// <param name="velocityECEF">Velocity in the ECEF frame (m/s).</param>
    /// <returns>Coriolis acceleration in the ECEF frame (m/s²).</returns>
    public static Vector Acceleration(Vector velocityECEF)
    {
        double omega = EarthModel.RotationRate;

        // ω × v where ω = (0, 0, ω)
        // ω × v = (ω·vy, -ω·vx, 0) ... actually:
        // ω × v = |i  j  k |
        //         |0  0  ω |
        //         |vx vy vz|
        // = (0*vz - ω*vy, ω*vx - 0*vz, 0*vy - 0*vx)
        // = (-ω·vy, ω·vx, 0)

        // Coriolis: -2(ω × v)
        double ax = 2.0 * omega * velocityECEF.y;
        double ay = -2.0 * omega * velocityECEF.x;
        double az = 0.0;

        return new Vector(ax, ay, az);
    }

    /// <summary>
    /// Computes the centrifugal acceleration at a position in the ECEF frame.
    /// a_centrifugal = -ω × (ω × r) = ω²·(x, y, 0) (outward from rotation axis).
    /// </summary>
    /// <param name="positionECEF">Position in the ECEF frame (meters).</param>
    /// <returns>Centrifugal acceleration in the ECEF frame (m/s²).</returns>
    public static Vector CentrifugalAcceleration(Vector positionECEF)
    {
        double omega2 = EarthModel.RotationRate * EarthModel.RotationRate;

        return new Vector(omega2 * positionECEF.x, omega2 * positionECEF.y, 0.0);
    }

    /// <summary>
    /// Computes the total fictitious acceleration (Coriolis + centrifugal) for a body in the ECEF frame.
    /// </summary>
    public static Vector TotalFictitiousAcceleration(Vector positionECEF, Vector velocityECEF)
    {
        Vector coriolis = Acceleration(velocityECEF);
        Vector centrifugal = CentrifugalAcceleration(positionECEF);
        return new Vector(
            coriolis.x + centrifugal.x,
            coriolis.y + centrifugal.y,
            coriolis.z + centrifugal.z);
    }
}
