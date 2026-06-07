using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.RL;

/// <summary>
/// Multi-agent pursuit-evasion RL environment using the flight dynamics engine.
///
/// Two aircraft — a pursuer and an evader — compete in a 3D airspace.
/// The environment is formulated from the pursuer's perspective.
///
/// <b>Observation (18D):</b>
///   Pursuer state (6): [altitude, airspeed, roll, pitch, yaw, alpha]
///   Evader relative (6): [rel_x, rel_y, rel_z, rel_vx, rel_vy, rel_vz] (normalized)
///   Pursuer rates (3): [p, q, r]
///   Engagement geometry (3): [range, aspect_angle, closing_speed]
///
/// <b>Action (4D continuous):</b>
///   [throttle, pitch, roll, yaw]
///
/// <b>Reward:</b>
///   Positive for reducing range, maintaining firing solution (aspect angle),
///   large bonus for "kill" (range &lt; kill radius and aspect &lt; threshold).
///   Penalty for losing too much energy or crashing.
///
/// The evader follows a simple scripted evasion policy (random turns + altitude changes).
/// </summary>
public class DogfightEnv : IEnvironment
{
    public int ObservationSize => 18;
    public int ActionSize => 4;
    public bool IsDiscrete => false;

    private readonly FlightDynamicsEngine _pursuer;
    private readonly FlightDynamicsEngine _evader;
    private readonly AircraftConfig _config;
    private readonly double _dt;
    private readonly int _maxSteps;
    private readonly double _killRadius;
    private readonly double _killAspect;
    private readonly double _initSeparation;

    private int _step;
    private Random _rng;
    private double _prevRange;

    // Evader maneuver state
    private double _evaderRoll;
    private double _evaderPitch;
    private int _evaderManeuverTimer;

    /// <summary>
    /// Creates a dogfight environment.
    /// </summary>
    /// <param name="config">Aircraft config for both aircraft.</param>
    /// <param name="dt">Simulation timestep.</param>
    /// <param name="maxSteps">Maximum steps per episode.</param>
    /// <param name="killRadius">Range for a successful intercept (m).</param>
    /// <param name="killAspect">Maximum aspect angle for a kill (rad). ~30 deg default.</param>
    /// <param name="initSeparation">Initial separation distance (m).</param>
    public DogfightEnv(
        AircraftConfig config = null,
        double dt = 0.05,
        int maxSteps = 3000,
        double killRadius = 300.0,
        double killAspect = 0.52, // ~30 degrees
        double initSeparation = 3000.0)
    {
        _config = config ?? AircraftConfig.GenericJet();
        _pursuer = new FlightDynamicsEngine(_config);
        _evader = new FlightDynamicsEngine(_config);
        _dt = dt;
        _maxSteps = maxSteps;
        _killRadius = killRadius;
        _killAspect = killAspect;
        _initSeparation = initSeparation;
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        _step = 0;

        double alt = 3000.0;
        double speed = 150.0;

        // Pursuer starts behind the evader
        _pursuer.Reset();
        _pursuer.Init();
        _pursuer.SetState(new AircraftState(
            new Vector(0, 0, -alt),
            new Vector(speed, 0, 0),
            Quaternion.Identity,
            new Vector(0, 0, 0)));

        // Evader starts ahead with slight offset
        double offsetY = (_rng.NextDouble() - 0.5) * 1000;
        _evader.Reset();
        _evader.Init();
        _evader.SetState(new AircraftState(
            new Vector(_initSeparation, offsetY, -alt),
            new Vector(speed, 0, 0),
            Quaternion.Identity,
            new Vector(0, 0, 0)));

        _evaderRoll = 0;
        _evaderPitch = 0;
        _evaderManeuverTimer = 0;

        _prevRange = GetRange();

        return (GetObservation(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        // Pursuer action
        double throttle = Math.Clamp(action[0], 0, 1);
        double pitch = Math.Clamp(action[1], -1, 1);
        double roll = Math.Clamp(action[2], -1, 1);
        double yaw = Math.Clamp(action[3], -1, 1);

        _pursuer.SetInput(new ControlInput(throttle, pitch, roll, yaw));
        _pursuer.Step(_dt);

        // Evader: scripted evasive maneuvers
        StepEvader();

        _step++;

        // Compute engagement geometry
        double range = GetRange();
        double aspectAngle = GetAspectAngle();
        double closingSpeed = _prevRange - range;

        // Reward
        double reward = 0;

        // Reward closing speed
        reward += closingSpeed * 0.01;

        // Reward being within firing envelope
        if (range < _killRadius * 3 && aspectAngle < _killAspect * 2)
            reward += 0.5;

        // Kill achieved
        bool kill = range < _killRadius && aspectAngle < _killAspect;
        if (kill) reward += 50.0;

        // Energy management: penalize very low airspeed
        if (_pursuer.State.Airspeed < 50) reward -= 1.0;

        // Crash penalties
        bool pursuerCrashed = _pursuer.State.Altitude < 0;
        bool evaderCrashed = _evader.State.Altitude < 0;
        if (pursuerCrashed) reward -= 100.0;

        // Range too large (lost contact)
        bool lost = range > _initSeparation * 3;
        if (lost) reward -= 10.0;

        _prevRange = range;

        bool done = kill || pursuerCrashed || evaderCrashed || lost || _step >= _maxSteps;

        var info = new Dictionary<string, object>
        {
            { "range", range },
            { "aspect", aspectAngle },
            { "kill", kill },
            { "pursuer_alt", _pursuer.State.Altitude },
            { "evader_alt", _evader.State.Altitude }
        };

        return (GetObservation(), reward, done, info);
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        // Map discrete to continuous
        double throttle = (action & 1) != 0 ? 0.9 : 0.5;
        double pitch = ((action >> 1) & 3) switch { 1 => 0.5, 2 => -0.5, _ => 0 };
        double roll = ((action >> 3) & 3) switch { 1 => 0.7, 2 => -0.7, _ => 0 };
        return Step(new VectorN(new[] { throttle, pitch, roll, 0.0 }));
    }

    private VectorN GetObservation()
    {
        var ps = _pursuer.State;
        var es = _evader.State;
        var pEuler = ps.EulerAngles;

        // Relative position and velocity
        var relPos = es.Position - ps.Position;
        var relVel = es.Velocity - ps.Velocity;

        double range = relPos.GetMagnitude();
        double aspect = GetAspectAngle();
        double closingSpeed = range > 0.1 ? -(relVel.x * relPos.x + relVel.y * relPos.y + relVel.z * relPos.z) / range : 0;

        return new VectorN(new[]
        {
            // Pursuer state (6)
            ps.Altitude / 5000.0,
            ps.Airspeed / 300.0,
            pEuler.roll,
            pEuler.pitch,
            pEuler.yaw,
            ps.AngleOfAttack,
            // Relative state (6)
            relPos.x / 5000.0,
            relPos.y / 5000.0,
            relPos.z / 5000.0,
            relVel.x / 200.0,
            relVel.y / 200.0,
            relVel.z / 200.0,
            // Rates (3)
            ps.AngularRate.x,
            ps.AngularRate.y,
            ps.AngularRate.z,
            // Engagement (3)
            Math.Min(range / 5000.0, 2.0),
            aspect / Math.PI,
            closingSpeed / 200.0
        });
    }

    private void StepEvader()
    {
        _evaderManeuverTimer--;
        if (_evaderManeuverTimer <= 0)
        {
            // New random maneuver
            _evaderRoll = (_rng.NextDouble() - 0.5) * 1.6;  // [-0.8, 0.8]
            _evaderPitch = (_rng.NextDouble() - 0.5) * 0.6;  // [-0.3, 0.3]
            _evaderManeuverTimer = 20 + _rng.Next(60);        // hold for 1-4 seconds
        }

        // Altitude safety: pull up if too low
        double evPitch = _evaderPitch;
        if (_evader.State.Altitude < 500) evPitch = 0.5;

        _evader.SetInput(new ControlInput(0.8, evPitch, _evaderRoll, 0));
        _evader.Step(_dt);
    }

    private double GetRange()
    {
        var relPos = _evader.State.Position - _pursuer.State.Position;
        return relPos.GetMagnitude();
    }

    private double GetAspectAngle()
    {
        var relPos = _evader.State.Position - _pursuer.State.Position;
        double range = relPos.GetMagnitude();
        if (range < 0.1) return 0;

        // Aspect angle: angle between pursuer's velocity vector and line-of-sight to evader
        var pVel = _pursuer.State.Velocity;
        double pSpeed = pVel.GetMagnitude();
        if (pSpeed < 0.1) return Math.PI;

        double dot = (pVel.x * relPos.x + pVel.y * relPos.y + pVel.z * relPos.z) / (pSpeed * range);
        return Math.Acos(Math.Clamp(dot, -1, 1));
    }
}
