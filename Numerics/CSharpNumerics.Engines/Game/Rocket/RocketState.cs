using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Full state vector for 6DOF rocket simulation with variable mass.
/// 14-state: position(3), velocity(3), attitude quaternion(4), angular rates(3), mass(1).
/// 
/// Conventions:
///   World frame: x = North, y = East, z = Down (NED) during atmospheric flight,
///   or ECI (x = vernal equinox, z = north pole) for orbital operations.
///   Body frame:  x = Forward (nose), y = Right (starboard), z = Down (belly).
/// </summary>
public class RocketState
{
    /// <summary>Position in world frame (m).</summary>
    public Vector Position { get; set; }

    /// <summary>Velocity in world frame (m/s).</summary>
    public Vector Velocity { get; set; }

    /// <summary>Attitude quaternion (body → world rotation).</summary>
    public Quaternion Attitude { get; set; }

    /// <summary>Angular velocity in body frame (rad/s): (p = roll, q = pitch, r = yaw).</summary>
    public Vector AngularRate { get; set; }

    /// <summary>Current total vehicle mass including propellant (kg).</summary>
    public double Mass { get; set; }

    /// <summary>
    /// Creates a rocket state with defaults: at origin, vertical, at rest.
    /// </summary>
    public RocketState()
    {
        Position = new Vector(0, 0, 0);
        Velocity = new Vector(0, 0, 0);
        Attitude = Quaternion.Identity;
        AngularRate = new Vector(0, 0, 0);
        Mass = 0;
    }

    /// <summary>
    /// Creates a rocket state from explicit values.
    /// </summary>
    public RocketState(Vector position, Vector velocity, Quaternion attitude, Vector angularRate, double mass)
    {
        Position = position;
        Velocity = velocity;
        Attitude = attitude;
        AngularRate = angularRate;
        Mass = mass;
    }

    /// <summary>Speed in m/s (magnitude of velocity).</summary>
    public double Speed => Velocity.GetMagnitude();

    /// <summary>Altitude (negative of z in NED convention, so positive = up).</summary>
    public double Altitude => -Position.z;

    /// <summary>
    /// Creates a deep copy of the current state.
    /// </summary>
    public RocketState Clone()
    {
        return new RocketState(Position, Velocity, Attitude, AngularRate, Mass);
    }

    /// <summary>
    /// Packs state into a VectorN for numerical integration (14 elements).
    /// [x, y, z, vx, vy, vz, qw, qx, qy, qz, p, q, r, mass]
    /// </summary>
    public VectorN ToVector()
    {
        return new VectorN(new double[]
        {
            Position.x, Position.y, Position.z,
            Velocity.x, Velocity.y, Velocity.z,
            Attitude.w, Attitude.x, Attitude.y, Attitude.z,
            AngularRate.x, AngularRate.y, AngularRate.z,
            Mass
        });
    }

    /// <summary>
    /// Unpacks state from a VectorN (14 elements).
    /// </summary>
    public static RocketState FromVector(VectorN v)
    {
        return new RocketState(
            position: new Vector(v[0], v[1], v[2]),
            velocity: new Vector(v[3], v[4], v[5]),
            attitude: new Quaternion(v[6], v[7], v[8], v[9]).Normalize(),
            angularRate: new Vector(v[10], v[11], v[12]),
            mass: v[13]);
    }
}
