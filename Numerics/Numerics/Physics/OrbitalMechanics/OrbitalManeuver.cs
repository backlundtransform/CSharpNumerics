using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;

namespace CSharpNumerics.Physics.OrbitalMechanics;

/// <summary>
/// Represents a generalized impulsive orbital maneuver (instantaneous velocity change).
/// </summary>
public class OrbitalManeuver
{
    /// <summary>ΔV magnitude in m/s.</summary>
    public double DeltaV { get; }

    /// <summary>Direction of the burn in ECI frame (unit vector).</summary>
    public Vector Direction { get; }

    /// <summary>Time at which the maneuver occurs (seconds from epoch).</summary>
    public double Time { get; }

    /// <summary>Description/label for the maneuver.</summary>
    public string Label { get; }

    public OrbitalManeuver(double deltaV, Vector direction, double time, string label = "")
    {
        DeltaV = deltaV;
        Direction = (1.0 / direction.GetMagnitude()) * direction;
        Time = time;
        Label = label;
    }

    /// <summary>
    /// Creates a prograde maneuver (in the direction of current velocity).
    /// </summary>
    public static OrbitalManeuver Prograde(double deltaV, Vector currentVelocity, double time)
    {
        Vector dir = (1.0 / currentVelocity.GetMagnitude()) * currentVelocity;
        return new OrbitalManeuver(deltaV, dir, time, "Prograde");
    }

    /// <summary>
    /// Creates a retrograde maneuver (opposite to current velocity).
    /// </summary>
    public static OrbitalManeuver Retrograde(double deltaV, Vector currentVelocity, double time)
    {
        Vector dir = (-1.0 / currentVelocity.GetMagnitude()) * currentVelocity;
        return new OrbitalManeuver(deltaV, dir, time, "Retrograde");
    }

    /// <summary>
    /// Applies this maneuver to a velocity vector, returning the new velocity.
    /// </summary>
    public Vector Apply(Vector currentVelocity)
    {
        return new Vector(
            currentVelocity.x + DeltaV * Direction.x,
            currentVelocity.y + DeltaV * Direction.y,
            currentVelocity.z + DeltaV * Direction.z);
    }
}

/// <summary>
/// Computes the ΔV required to circularize an orbit at a given point.
/// </summary>
public static class CircularizationBurn
{
    /// <summary>
    /// Computes the ΔV needed to circularize at the current orbital radius.
    /// Positive ΔV = prograde burn (speed up), negative = retrograde (slow down).
    /// </summary>
    /// <param name="position">Current position in ECI (meters).</param>
    /// <param name="velocity">Current velocity in ECI (m/s).</param>
    /// <param name="mu">Gravitational parameter (m³/s²). Default: Earth.</param>
    /// <returns>ΔV needed (m/s). Positive = prograde.</returns>
    public static double Compute(Vector position, Vector velocity, double mu = 0)
    {
        if (mu == 0) mu = EarthModel.GM;

        double r = position.GetMagnitude();
        double vCurrent = velocity.GetMagnitude();
        double vCircular = Math.Sqrt(mu / r);

        return vCircular - vCurrent;
    }

    /// <summary>
    /// Computes the ΔV to circularize at apoapsis of the current orbit.
    /// </summary>
    /// <param name="elements">Current orbital elements.</param>
    /// <returns>ΔV needed at apoapsis (m/s). Typically positive for Hohmann transfer completion.</returns>
    public static double AtApoapsis(OrbitalElements elements)
    {
        double ra = elements.Apoapsis;
        double mu = elements.Mu;

        // Velocity at apoapsis on current orbit: v_a = h/r_a
        double h = elements.SpecificAngularMomentum;
        double vApoapsis = h / ra;

        // Circular velocity at apoapsis radius
        double vCircular = Math.Sqrt(mu / ra);

        return vCircular - vApoapsis;
    }

    /// <summary>
    /// Computes the ΔV to circularize at periapsis of the current orbit.
    /// </summary>
    /// <param name="elements">Current orbital elements.</param>
    /// <returns>ΔV needed at periapsis (m/s). Typically negative (retrograde).</returns>
    public static double AtPeriapsis(OrbitalElements elements)
    {
        double rp = elements.Periapsis;
        double mu = elements.Mu;

        // Velocity at periapsis on current orbit: v_p = h/r_p
        double h = elements.SpecificAngularMomentum;
        double vPeriapsis = h / rp;

        // Circular velocity at periapsis radius
        double vCircular = Math.Sqrt(mu / rp);

        return vCircular - vPeriapsis;
    }
}
