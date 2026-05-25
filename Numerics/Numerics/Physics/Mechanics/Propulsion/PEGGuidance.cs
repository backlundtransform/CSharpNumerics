using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Powered Explicit Guidance (PEG) — an iterative guidance algorithm that steers
/// the rocket to achieve target orbital parameters. Based on the linear tangent
/// steering law used by the Space Shuttle and modern launch vehicles.
/// 
/// PEG iteratively computes the thrust direction to null out the velocity-to-be-gained
/// (Vgo) at burnout, targeting a specific orbit.
/// </summary>
public class PEGGuidance
{
    /// <summary>Target semi-major axis (meters).</summary>
    public double TargetSemiMajorAxis { get; set; }

    /// <summary>Target eccentricity. Default 0 (circular).</summary>
    public double TargetEccentricity { get; set; } = 0;

    /// <summary>Target inclination (radians).</summary>
    public double TargetInclination { get; set; }

    /// <summary>Gravitational parameter μ (m³/s²). Default: Earth.</summary>
    public double Mu { get; set; } = EarthModel.GM;

    /// <summary>Number of PEG iterations per guidance cycle. Default 3.</summary>
    public int Iterations { get; set; } = 3;

    /// <summary>Whether the guidance has converged.</summary>
    public bool HasConverged { get; private set; }

    /// <summary>Velocity-to-be-gained magnitude (m/s).</summary>
    public double VelocityToGain { get; private set; }

    /// <summary>Estimated time to cutoff (seconds).</summary>
    public double TimeToCutoff { get; private set; }

    private Vector _thrustDirection;
    private double _prevVgo = double.MaxValue;

    /// <summary>
    /// Computes the desired thrust direction given current state and engine parameters.
    /// </summary>
    /// <param name="position">Current ECI position (meters).</param>
    /// <param name="velocity">Current ECI velocity (m/s).</param>
    /// <param name="thrustAccel">Current thrust acceleration magnitude (m/s²).</param>
    /// <param name="exhaustVelocity">Effective exhaust velocity Ve = Isp*g0 (m/s).</param>
    /// <returns>Unit vector in ECI frame indicating desired thrust direction.</returns>
    public Vector ComputeThrustDirection(Vector position, Vector velocity,
        double thrustAccel, double exhaustVelocity)
    {
        double r = position.GetMagnitude();
        double v = velocity.GetMagnitude();

        // Target velocity at current radius for desired orbit
        double rTarget = TargetSemiMajorAxis * (1.0 - TargetEccentricity);
        double vTarget = Math.Sqrt(Mu * (2.0 / r - 1.0 / TargetSemiMajorAxis));

        // Desired velocity direction: tangent to target orbit
        // For circular orbit, perpendicular to radius
        Vector rUnit = (1.0 / r) * position;
        Vector hDesired = ComputeDesiredAngularMomentum(position, velocity);
        double hMag = hDesired.GetMagnitude();
        Vector hUnit = hMag > 1e-10 ? (1.0 / hMag) * hDesired : new Vector(0, 0, 1);

        // Desired velocity direction: h × r (perpendicular to radius, in orbital plane)
        Vector vDesiredDir = Cross(hUnit, rUnit);
        double vDesiredDirMag = vDesiredDir.GetMagnitude();
        if (vDesiredDirMag > 1e-10)
            vDesiredDir = (1.0 / vDesiredDirMag) * vDesiredDir;

        // Velocity to be gained
        Vector vDesired = vTarget * vDesiredDir;
        Vector vgo = new Vector(
            vDesired.x - velocity.x,
            vDesired.y - velocity.y,
            vDesired.z - velocity.z);

        VelocityToGain = vgo.GetMagnitude();

        // Iterative convergence check
        if (Math.Abs(VelocityToGain - _prevVgo) / (VelocityToGain + 1e-10) < 0.01)
            HasConverged = true;
        _prevVgo = VelocityToGain;

        // Time to cutoff estimate (Tsiolkovsky)
        if (thrustAccel > 0.01)
            TimeToCutoff = exhaustVelocity / thrustAccel * (1.0 - Math.Exp(-VelocityToGain / exhaustVelocity));
        else
            TimeToCutoff = double.PositiveInfinity;

        // Thrust direction: linear tangent steering
        // Apply gravity bias: steer slightly above Vgo to compensate for gravity loss during burn
        double tgo = TimeToCutoff;
        Vector gravBias = new Vector(0, 0, 0);
        if (tgo > 0 && tgo < 10000)
        {
            // Approximate gravity loss direction (toward center)
            Vector gDir = (-1.0 / r) * position;
            double gMag = Mu / (r * r);
            gravBias = (0.5 * gMag * tgo) * gDir;
        }

        Vector thrustDir = new Vector(
            vgo.x + gravBias.x,
            vgo.y + gravBias.y,
            vgo.z + gravBias.z);

        double thrustDirMag = thrustDir.GetMagnitude();
        if (thrustDirMag < 1e-10) return (1.0 / v) * velocity; // fallback: prograde

        _thrustDirection = (1.0 / thrustDirMag) * thrustDir;
        return _thrustDirection;
    }

    /// <summary>Resets convergence state for a new guidance phase.</summary>
    public void Reset()
    {
        HasConverged = false;
        _prevVgo = double.MaxValue;
        VelocityToGain = 0;
        TimeToCutoff = 0;
    }

    private Vector ComputeDesiredAngularMomentum(Vector position, Vector velocity)
    {
        // Current angular momentum direction (defines orbital plane)
        Vector hCurrent = Cross(position, velocity);
        double hMag = hCurrent.GetMagnitude();

        if (hMag < 1e-10)
            return new Vector(0, 0, 1); // Default to equatorial

        // For now, maintain current orbital plane (inclination change requires plane change maneuver)
        return hCurrent;
    }

    private static Vector Cross(Vector a, Vector b)
    {
        return new Vector(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x);
    }
}
