using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics.Propulsion;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Selects the appropriate guidance law per mission phase and computes the
/// commanded attitude each simulation step.
/// </summary>
public class GuidanceComputer
{
    private readonly GravityTurnGuidance _gravityTurn = new GravityTurnGuidance();
    private readonly PEGGuidance _peg = new PEGGuidance();
    private readonly AttitudeController _attitudeController = new AttitudeController();

    /// <summary>Gravity turn guidance instance (configurable).</summary>
    public GravityTurnGuidance GravityTurn => _gravityTurn;

    /// <summary>PEG guidance instance (configurable).</summary>
    public PEGGuidance PEG => _peg;

    /// <summary>Attitude controller instance (configurable).</summary>
    public AttitudeController AttitudeController => _attitudeController;

    /// <summary>Active guidance mode.</summary>
    public GuidanceMode Mode { get; set; } = GuidanceMode.GravityTurn;

    /// <summary>Altitude at which guidance switches from gravity turn to PEG (meters). Default 80 km.</summary>
    public double PEGActivationAltitude { get; set; } = 80000.0;

    /// <summary>Last computed commanded attitude.</summary>
    public Quaternion CommandedAttitude { get; private set; } = Quaternion.Identity;

    /// <summary>Last computed torque command from attitude controller.</summary>
    public Vector TorqueCommand { get; private set; }

    /// <summary>
    /// Updates guidance each step: selects mode, computes commanded attitude, and
    /// generates torque commands.
    /// </summary>
    /// <param name="state">Current rocket state.</param>
    /// <param name="time">Current simulation time (s).</param>
    /// <param name="altitude">Current altitude (m).</param>
    /// <param name="thrustAccel">Current thrust acceleration (m/s²).</param>
    /// <param name="exhaustVelocity">Effective exhaust velocity (m/s).</param>
    /// <param name="dt">Time step (s).</param>
    public void Update(RocketState state, double time, double altitude,
        double thrustAccel, double exhaustVelocity, double dt)
    {
        // Auto mode selection based on altitude
        if (Mode == GuidanceMode.Auto)
        {
            if (altitude < PEGActivationAltitude)
                Mode = GuidanceMode.GravityTurn;
            else
                Mode = GuidanceMode.PEG;
        }

        switch (Mode)
        {
            case GuidanceMode.GravityTurn:
                CommandedAttitude = _gravityTurn.ComputeAttitude(altitude, state.Velocity, time);
                if (altitude >= PEGActivationAltitude)
                    Mode = GuidanceMode.PEG;
                break;

            case GuidanceMode.PEG:
                Vector thrustDir = _peg.ComputeThrustDirection(
                    state.Position, state.Velocity, thrustAccel, exhaustVelocity);
                CommandedAttitude = DirectionToQuaternion(thrustDir);
                break;

            case GuidanceMode.Prograde:
                double speed = state.Speed;
                if (speed > 1.0)
                    CommandedAttitude = DirectionToQuaternion((1.0 / speed) * state.Velocity);
                break;

            case GuidanceMode.Hold:
                // Keep current commanded attitude (no change)
                break;
        }

        // Compute torque to track commanded attitude
        TorqueCommand = _attitudeController.ComputeTorque(
            state.Attitude, CommandedAttitude, state.AngularRate, dt);
    }

    /// <summary>Resets all guidance subsystems.</summary>
    public void Reset()
    {
        _gravityTurn.Reset();
        _peg.Reset();
        _attitudeController.Reset();
        Mode = GuidanceMode.GravityTurn;
        CommandedAttitude = Quaternion.Identity;
        TorqueCommand = new Vector(0, 0, 0);
    }

    private static Quaternion DirectionToQuaternion(Vector direction)
    {
        double mag = direction.GetMagnitude();
        if (mag < 1e-10) return Quaternion.Identity;

        Vector dir = (1.0 / mag) * direction;
        Vector bodyX = new Vector(1, 0, 0);

        Vector cross = new Vector(
            bodyX.y * dir.z - bodyX.z * dir.y,
            bodyX.z * dir.x - bodyX.x * dir.z,
            bodyX.x * dir.y - bodyX.y * dir.x);
        double dot = bodyX.x * dir.x + bodyX.y * dir.y + bodyX.z * dir.z;
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
/// Guidance modes for the guidance computer.
/// </summary>
public enum GuidanceMode
{
    /// <summary>Automatically selects based on altitude.</summary>
    Auto,
    /// <summary>Gravity turn (atmospheric ascent).</summary>
    GravityTurn,
    /// <summary>Powered Explicit Guidance (vacuum/upper stage).</summary>
    PEG,
    /// <summary>Track velocity vector (prograde).</summary>
    Prograde,
    /// <summary>Hold current attitude.</summary>
    Hold
}
