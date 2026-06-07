using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.AI;

/// <summary>
/// Multi-agent formation controller for coordinating wingmen.
/// Combines RL-learned station-keeping with formation geometry rules.
///
/// Each wingman maintains a desired offset from the leader in the leader's
/// body frame. The controller generates control inputs to maintain spacing
/// while responding to wind perturbation and leader maneuvers.
/// </summary>
public class FormationController
{
    /// <summary>Formation position definitions.</summary>
    public class FormationSlot
    {
        /// <summary>Desired offset from leader in leader's body frame (m). x=back, y=right, z=down.</summary>
        public Vector Offset { get; }

        /// <summary>Name of this slot.</summary>
        public string Name { get; }

        /// <summary>The flight engine for this wingman.</summary>
        public FlightDynamicsEngine Engine { get; set; }

        /// <summary>Optional RL agent for fine-tuning station-keeping.</summary>
        public IAgent RLAgent { get; set; }

        public FormationSlot(string name, Vector offset)
        {
            Name = name;
            Offset = offset;
        }
    }

    private readonly FlightDynamicsEngine _leader;
    private readonly List<FormationSlot> _slots = new();

    /// <summary>Proportional gain for position error → velocity command.</summary>
    public double PositionGain { get; set; } = 0.5;

    /// <summary>Proportional gain for velocity error → control input.</summary>
    public double VelocityGain { get; set; } = 0.3;

    /// <summary>Maximum allowed position error before aggressive correction (m).</summary>
    public double MaxPositionError { get; set; } = 100.0;

    /// <summary>Desired spacing tolerance (m). Within this, no correction needed.</summary>
    public double SpacingTolerance { get; set; } = 10.0;

    /// <summary>Formation slots.</summary>
    public IReadOnlyList<FormationSlot> Slots => _slots;

    /// <summary>
    /// Creates a formation controller for the given leader.
    /// </summary>
    public FormationController(FlightDynamicsEngine leader)
    {
        _leader = leader ?? throw new ArgumentNullException(nameof(leader));
    }

    /// <summary>
    /// Add a wingman to the formation.
    /// </summary>
    /// <param name="name">Slot name (e.g. "Wing 2").</param>
    /// <param name="offset">Desired offset in leader body frame.</param>
    /// <param name="engine">Wingman's flight engine.</param>
    /// <param name="rlAgent">Optional RL agent for fine control.</param>
    public FormationSlot AddWingman(string name, Vector offset, FlightDynamicsEngine engine, IAgent rlAgent = null)
    {
        var slot = new FormationSlot(name, offset) { Engine = engine, RLAgent = rlAgent };
        _slots.Add(slot);
        return slot;
    }

    /// <summary>
    /// Update all wingmen for one timestep.
    /// Computes desired position from leader state, generates control inputs,
    /// and steps each wingman's flight engine.
    /// </summary>
    /// <param name="dt">Timestep in seconds.</param>
    public void Step(double dt)
    {
        var leaderState = _leader.State;

        foreach (var slot in _slots)
        {
            var desiredWorldPos = ComputeDesiredPosition(leaderState, slot.Offset);
            var input = ComputeControlInput(slot, leaderState, desiredWorldPos);
            slot.Engine.SetInput(input);
            slot.Engine.Step(dt);
        }
    }

    /// <summary>
    /// Get the position error for each wingman (distance from desired position).
    /// </summary>
    public double[] GetPositionErrors()
    {
        var errors = new double[_slots.Count];
        var leaderState = _leader.State;

        for (int i = 0; i < _slots.Count; i++)
        {
            var desired = ComputeDesiredPosition(leaderState, _slots[i].Offset);
            var actual = _slots[i].Engine.State.Position;
            var diff = desired - actual;
            errors[i] = diff.GetMagnitude();
        }

        return errors;
    }

    /// <summary>
    /// Compute the desired world position from leader state and body-frame offset.
    /// </summary>
    private Vector ComputeDesiredPosition(AircraftState leaderState, Vector bodyOffset)
    {
        var worldOffset = FrameTransforms.BodyToWorld(leaderState.Attitude, bodyOffset);
        return leaderState.Position + worldOffset;
    }

    /// <summary>
    /// Generate control inputs for a wingman using proportional control.
    /// If an RL agent is attached, it refines the PD output.
    /// </summary>
    private ControlInput ComputeControlInput(FormationSlot slot, AircraftState leaderState, Vector desiredPos)
    {
        var wingState = slot.Engine.State;

        // Position error in world frame
        var posError = desiredPos - wingState.Position;
        double posErrorMag = posError.GetMagnitude();

        // Skip correction if within tolerance
        if (posErrorMag < SpacingTolerance)
        {
            // Match leader's throttle and attitude approximately
            return new ControlInput(0.5, 0, 0, 0);
        }

        // Desired velocity: proportional to position error, clamped
        double errorScale = Math.Min(posErrorMag / MaxPositionError, 1.0);
        double gain = PositionGain * errorScale * Math.Min(posErrorMag, 50.0) / Math.Max(posErrorMag, 0.01);
        var desiredVel = gain * posError;

        // Velocity error
        var velError = desiredVel - (wingState.Velocity - leaderState.Velocity);

        // Convert velocity error to body frame for control
        var velErrorBody = FrameTransforms.WorldToBody(wingState.Attitude, velError);

        // Simple proportional control mapping
        double throttle = Math.Clamp(0.5 + VelocityGain * velErrorBody.x / 50.0, 0, 1);
        double pitchCmd = Math.Clamp(-VelocityGain * velErrorBody.z / 20.0, -1, 1);
        double rollCmd = Math.Clamp(VelocityGain * velErrorBody.y / 20.0, -1, 1);
        double yawCmd = 0; // simplified — rudder for sideslip only

        // If RL agent is attached, let it refine
        if (slot.RLAgent != null)
        {
            var obs = BuildObservation(wingState, leaderState, posError, velError);
            var rlAction = slot.RLAgent.SelectContinuousAction(obs);

            // Blend: 70% PD + 30% RL
            throttle = Math.Clamp(throttle * 0.7 + rlAction[0] * 0.3, 0, 1);
            pitchCmd = Math.Clamp(pitchCmd * 0.7 + rlAction[1] * 0.3, -1, 1);
            rollCmd = Math.Clamp(rollCmd * 0.7 + rlAction[2] * 0.3, -1, 1);
            if (rlAction.Length > 3)
                yawCmd = Math.Clamp(rlAction[3] * 0.3, -1, 1);
        }

        return new ControlInput(throttle, pitchCmd, rollCmd, yawCmd);
    }

    private VectorN BuildObservation(AircraftState wingState, AircraftState leaderState,
        Vector posError, Vector velError)
    {
        var euler = wingState.EulerAngles;
        return new VectorN(new[]
        {
            posError.x / 100.0,
            posError.y / 100.0,
            posError.z / 100.0,
            velError.x / 50.0,
            velError.y / 50.0,
            velError.z / 50.0,
            euler.roll,
            euler.pitch,
            wingState.Airspeed / 100.0,
            wingState.AngularRate.x,
            wingState.AngularRate.y,
            wingState.AngularRate.z
        });
    }
}
