using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Implements gravity turn guidance: an initial pitch-over (kick) followed by
/// a gravity-assisted turn where the rocket follows its velocity vector.
/// This produces a near-optimal trajectory for atmospheric ascent.
/// </summary>
public class GravityTurnGuidance
{
    /// <summary>Altitude at which the pitch-over kick begins (meters). Default 500 m.</summary>
    public double KickAltitude { get; set; } = 500.0;

    /// <summary>Initial kick angle from vertical (radians). Typical: 2-5°. Default 3°.</summary>
    public double KickAngle { get; set; } = 3.0 * Math.PI / 180.0;

    /// <summary>Duration of the pitch-over maneuver (seconds). Default 10 s.</summary>
    public double KickDuration { get; set; } = 10.0;

    /// <summary>Target orbit inclination (radians). Determines azimuth heading.</summary>
    public double TargetInclination { get; set; } = 0.0;

    /// <summary>Launch site latitude (radians). Used to compute launch azimuth.</summary>
    public double LaunchLatitude { get; set; } = 0.0;

    private double _kickStartTime = -1;
    private bool _kickComplete;

    /// <summary>Current guidance phase.</summary>
    public GravityTurnPhase Phase { get; private set; } = GravityTurnPhase.VerticalAscent;

    /// <summary>
    /// Computes the commanded attitude quaternion given current flight state.
    /// </summary>
    /// <param name="altitude">Current altitude (meters).</param>
    /// <param name="velocity">Current velocity vector (world frame).</param>
    /// <param name="time">Current simulation time (seconds).</param>
    /// <returns>Commanded attitude quaternion.</returns>
    public Quaternion ComputeAttitude(double altitude, Vector velocity, double time)
    {
        double speed = velocity.GetMagnitude();

        if (altitude < KickAltitude)
        {
            // Phase 1: Vertical ascent (straight up)
            Phase = GravityTurnPhase.VerticalAscent;
            return VerticalAttitude();
        }

        if (!_kickComplete)
        {
            // Phase 2: Pitch-over kick
            Phase = GravityTurnPhase.PitchOver;
            if (_kickStartTime < 0) _kickStartTime = time;

            double elapsed = time - _kickStartTime;
            double fraction = Math.Min(elapsed / KickDuration, 1.0);
            double currentKickAngle = fraction * KickAngle;

            if (fraction >= 1.0) _kickComplete = true;

            return PitchedAttitude(currentKickAngle);
        }

        // Phase 3: Gravity turn — follow velocity vector (prograde)
        Phase = GravityTurnPhase.GravityTurn;
        if (speed < 1.0) return VerticalAttitude();

        return ProgradeAttitude(velocity);
    }

    /// <summary>
    /// Computes the required launch azimuth to achieve target inclination.
    /// azimuth = asin(cos(i) / cos(lat))
    /// </summary>
    public double LaunchAzimuth()
    {
        double cosI = Math.Cos(TargetInclination);
        double cosLat = Math.Cos(LaunchLatitude);
        if (Math.Abs(cosLat) < 1e-10) return 0;

        double sinAz = cosI / cosLat;
        if (sinAz > 1.0) sinAz = 1.0;
        if (sinAz < -1.0) sinAz = -1.0;
        return Math.Asin(sinAz);
    }

    /// <summary>Resets guidance state for a new ascent.</summary>
    public void Reset()
    {
        _kickStartTime = -1;
        _kickComplete = false;
        Phase = GravityTurnPhase.VerticalAscent;
    }

    private static Quaternion VerticalAttitude()
    {
        // Nose-up in NED: body X → world -Z (up)
        return Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2);
    }

    private static Quaternion PitchedAttitude(double pitchFromVertical)
    {
        // Pitch over from vertical by the given angle (toward +X/north)
        double angle = Math.PI / 2 - pitchFromVertical;
        return Quaternion.FromAxisAngle(new Vector(0, 1, 0), angle);
    }

    private static Quaternion ProgradeAttitude(Vector velocity)
    {
        // Align body X-axis with velocity vector
        Vector velUnit = (1.0 / velocity.GetMagnitude()) * velocity;

        // Construct quaternion that rotates (1,0,0) to velUnit
        Vector bodyX = new Vector(1, 0, 0);
        Vector cross = new Vector(
            bodyX.y * velUnit.z - bodyX.z * velUnit.y,
            bodyX.z * velUnit.x - bodyX.x * velUnit.z,
            bodyX.x * velUnit.y - bodyX.y * velUnit.x);
        double dot = bodyX.x * velUnit.x + bodyX.y * velUnit.y + bodyX.z * velUnit.z;

        double crossMag = cross.GetMagnitude();
        if (crossMag < 1e-10)
        {
            if (dot > 0) return Quaternion.Identity;
            return Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI);
        }

        double angle = Math.Atan2(crossMag, dot);
        Vector axis = (1.0 / crossMag) * cross;
        return Quaternion.FromAxisAngle(axis, angle);
    }
}

/// <summary>
/// Phases of a gravity turn ascent.
/// </summary>
public enum GravityTurnPhase
{
    VerticalAscent,
    PitchOver,
    GravityTurn
}
