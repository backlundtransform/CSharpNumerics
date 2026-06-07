using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Game.Flight;

/// <summary>
/// Full 12-state vector for a 6DOF aircraft:
///   Position     (x, y, z)           — metres in world frame
///   Velocity     (vx, vy, vz)        — m/s in world frame
///   Attitude     quaternion (w,x,y,z) — body orientation relative to world
///   AngularRate  (p, q, r)           — rad/s in body frame (roll, pitch, yaw rates)
/// 
/// Conventions:
///   World frame: x = North, y = East, z = Down (NED)
///   Body frame:  x = Forward, y = Right, z = Down
/// </summary>
public class AircraftState
{
    /// <summary>Position in world frame (m).</summary>
    public Vector Position { get; set; }

    /// <summary>Velocity in world frame (m/s).</summary>
    public Vector Velocity { get; set; }

    /// <summary>Attitude quaternion (body → world rotation).</summary>
    public Quaternion Attitude { get; set; }

    /// <summary>Angular velocity in body frame (rad/s): (p = roll rate, q = pitch rate, r = yaw rate).</summary>
    public Vector AngularRate { get; set; }

    /// <summary>
    /// Creates an aircraft state with all values at zero / identity.
    /// </summary>
    public AircraftState()
    {
        Position = new Vector(0, 0, 0);
        Velocity = new Vector(0, 0, 0);
        Attitude = Quaternion.Identity;
        AngularRate = new Vector(0, 0, 0);
    }

    /// <summary>
    /// Creates an aircraft state from explicit values.
    /// </summary>
    public AircraftState(Vector position, Vector velocity, Quaternion attitude, Vector angularRate)
    {
        Position = position;
        Velocity = velocity;
        Attitude = attitude;
        AngularRate = angularRate;
    }

    /// <summary>
    /// True airspeed (magnitude of world-frame velocity) in m/s.
    /// </summary>
    public double Airspeed => Velocity.GetMagnitude();

    /// <summary>
    /// Altitude (negative of z in NED convention, so positive = up).
    /// </summary>
    public double Altitude => -Position.z;

    /// <summary>
    /// Velocity expressed in body frame.
    /// </summary>
    public Vector BodyVelocity => FrameTransforms.WorldToBody(Attitude, Velocity);

    /// <summary>
    /// Angle of attack computed from body-frame velocity.
    /// </summary>
    public double AngleOfAttack => FrameTransforms.AngleOfAttack(BodyVelocity);

    /// <summary>
    /// Sideslip angle computed from body-frame velocity.
    /// </summary>
    public double SideslipAngle => FrameTransforms.SideslipAngle(BodyVelocity);

    /// <summary>
    /// Euler angles (roll, pitch, yaw) in radians extracted from the attitude quaternion.
    /// </summary>
    public (double roll, double pitch, double yaw) EulerAngles => Attitude.ToEulerAngles();

    /// <summary>
    /// Creates a deep copy of this state.
    /// </summary>
    public AircraftState Clone()
    {
        return new AircraftState(
            new Vector(Position.x, Position.y, Position.z),
            new Vector(Velocity.x, Velocity.y, Velocity.z),
            new Quaternion(Attitude.w, Attitude.x, Attitude.y, Attitude.z),
            new Vector(AngularRate.x, AngularRate.y, AngularRate.z));
    }
}
