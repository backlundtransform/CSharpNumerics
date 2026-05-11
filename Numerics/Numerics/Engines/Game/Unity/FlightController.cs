using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.Numerics.Objects;
using System;
using static CSharpNumerics.Engines.Game.Unity.UnityAdapter;

namespace CSharpNumerics.Engines.Game.Unity;

/// <summary>
/// Unity-facing flight controller interface.
///
/// Bridges Unity input (float axes) to <see cref="ControlInput"/> and
/// reads <see cref="AircraftState"/> to produce Unity-compatible transform data.
///
/// Usage in Unity:
///   1. Create FlightController with a FlightDynamicsEngine
///   2. Each FixedUpdate: SetInputAxis(), StepSimulation()
///   3. Each Update: read GetPosition(), GetRotation() for rendering
///
/// No Unity dependency — uses surrogate types.
/// </summary>
public class FlightController
{
    private readonly FlightDynamicsEngine _engine;
    private readonly ControlInput _input = new();

    /// <summary>The underlying flight dynamics engine.</summary>
    public FlightDynamicsEngine Engine => _engine;

    /// <summary>Current control input state.</summary>
    public ControlInput Input => _input;

    /// <summary>Current aircraft state.</summary>
    public AircraftState State => _engine.State;

    /// <summary>Simulation time in seconds.</summary>
    public double SimTime => _engine.Time;

    /// <summary>
    /// Creates a flight controller for a flight dynamics engine.
    /// </summary>
    public FlightController(FlightDynamicsEngine engine)
    {
        _engine = engine ?? throw new ArgumentNullException(nameof(engine));
    }

    /// <summary>
    /// Set a control axis value. Maps to ControlInput properties.
    /// </summary>
    /// <param name="axis">Axis name: "throttle", "pitch", "roll", "yaw", "flaps".</param>
    /// <param name="value">Normalized value (0-1 for throttle/flaps, -1 to +1 for pitch/roll/yaw).</param>
    public void SetInputAxis(string axis, double value)
    {
        switch (axis?.ToLowerInvariant())
        {
            case "throttle": _input.Throttle = Math.Max(0, Math.Min(1, value)); break;
            case "pitch": _input.Pitch = Math.Max(-1, Math.Min(1, value)); break;
            case "roll": _input.Roll = Math.Max(-1, Math.Min(1, value)); break;
            case "yaw": _input.Yaw = Math.Max(-1, Math.Min(1, value)); break;
            case "flaps": _input.Flaps = Math.Max(0, Math.Min(1, value)); break;
        }
    }

    /// <summary>Set gear state.</summary>
    public void SetGear(bool down) => _input.GearDown = down;

    /// <summary>
    /// Advance the simulation by dt seconds.
    /// </summary>
    public void StepSimulation(double dt)
    {
        _engine.SetInput(_input);
        _engine.Step(dt);
    }

    /// <summary>
    /// Get the aircraft position in Unity coordinates (Y-up).
    /// NED → Unity: CSN(North, East, Down) → Unity(North, -Down, East) = (x, -z, y)
    /// </summary>
    public UnityVector3 GetPosition()
    {
        return ToUnityVector3(State.Position);
    }

    /// <summary>
    /// Get the aircraft rotation as a Unity quaternion.
    /// </summary>
    public UnityQuaternion GetRotation()
    {
        var q = State.Attitude;
        // Convert quaternion from NED body frame to Unity Y-up
        // Swap Y↔Z components for coordinate change
        return new UnityQuaternion((float)q.x, (float)q.z, (float)q.y, (float)q.w);
    }

    /// <summary>
    /// Get the aircraft velocity in Unity coordinates.
    /// </summary>
    public UnityVector3 GetVelocity()
    {
        return ToUnityVector3(State.Velocity);
    }

    /// <summary>
    /// Get airspeed in m/s.
    /// </summary>
    public double GetAirspeed() => State.Velocity.GetMagnitude();

    /// <summary>
    /// Get altitude (negative of NED Z position).
    /// </summary>
    public double GetAltitude() => -State.Position.z;

    /// <summary>
    /// Get heading in degrees (0 = North, clockwise).
    /// </summary>
    public double GetHeading()
    {
        double heading = Math.Atan2(State.Velocity.y, State.Velocity.x) * 180.0 / Math.PI;
        if (heading < 0) heading += 360;
        return heading;
    }

    /// <summary>
    /// Get a complete HUD data snapshot for UI display.
    /// </summary>
    public HUDData GetHUDData()
    {
        return new HUDData
        {
            Airspeed = GetAirspeed(),
            Altitude = GetAltitude(),
            Heading = GetHeading(),
            Throttle = _input.Throttle,
            GearDown = _input.GearDown,
            Pitch = _input.Pitch,
            Roll = _input.Roll
        };
    }

    /// <summary>
    /// HUD display data container.
    /// </summary>
    public struct HUDData
    {
        public double Airspeed;
        public double Altitude;
        public double Heading;
        public double Throttle;
        public bool GearDown;
        public double Pitch;
        public double Roll;
    }
}
