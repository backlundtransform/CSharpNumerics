using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.RL;

/// <summary>
/// RL environment wrapping the 6DOF flight dynamics engine.
///
/// <b>Observation (12D):</b>
///   [altitude, airspeed, alpha, beta, roll, pitch, yaw, p, q, r, dist_to_waypoint, heading_error]
///
/// <b>Action (4D continuous):</b>
///   [throttle (0–1), pitch (-1 to 1), roll (-1 to 1), yaw (-1 to 1)]
///
/// <b>Reward:</b>
///   - Positive for getting closer to the current waypoint
///   - Bonus for reaching a waypoint (distance &lt; threshold)
///   - Penalty for excessive bank angle, stall, or ground impact
///   - Small fuel-efficiency bonus for lower throttle
///
/// <b>Done:</b>
///   - All waypoints reached, or aircraft crashed (altitude &lt; 0), or max steps exceeded
/// </summary>
public class FlightEnv : IEnvironment
{
    public int ObservationSize => 12;
    public int ActionSize => 4;
    public bool IsDiscrete => false;

    private readonly FlightDynamicsEngine _engine;
    private readonly AircraftConfig _config;
    private readonly double _dt;
    private readonly int _maxSteps;
    private readonly double _waypointRadius;

    private List<Vector> _waypoints;
    private int _currentWaypoint;
    private int _step;
    private double _prevDistance;
    private Random _rng;

    // Initial conditions
    private double _initAltitude;
    private double _initAirspeed;

    /// <summary>
    /// Creates a flight environment.
    /// </summary>
    /// <param name="config">Aircraft configuration. Defaults to GenericLightAircraft.</param>
    /// <param name="waypoints">Waypoints to visit. If null, generates a simple circuit.</param>
    /// <param name="dt">Simulation timestep in seconds.</param>
    /// <param name="maxSteps">Maximum steps per episode.</param>
    /// <param name="waypointRadius">Distance threshold to count waypoint as reached (m).</param>
    /// <param name="initAltitude">Initial altitude in metres.</param>
    /// <param name="initAirspeed">Initial airspeed in m/s.</param>
    public FlightEnv(
        AircraftConfig config = null,
        List<Vector> waypoints = null,
        double dt = 0.05,
        int maxSteps = 2000,
        double waypointRadius = 200.0,
        double initAltitude = 1000.0,
        double initAirspeed = 60.0)
    {
        _config = config ?? AircraftConfig.GenericLightAircraft();
        _engine = new FlightDynamicsEngine(_config);
        _dt = dt;
        _maxSteps = maxSteps;
        _waypointRadius = waypointRadius;
        _initAltitude = initAltitude;
        _initAirspeed = initAirspeed;
        _waypoints = waypoints ?? DefaultWaypoints();
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        _step = 0;
        _currentWaypoint = 0;

        _engine.Reset();
        _engine.Init();

        // Set initial state: level flight at specified altitude and airspeed
        var initState = new AircraftState(
            position: new Vector(0, 0, -_initAltitude),  // NED: z-up is negative
            velocity: new Vector(_initAirspeed, 0, 0),   // flying north
            attitude: Quaternion.Identity,
            angularRate: new Vector(0, 0, 0));

        _engine.SetState(initState);
        _engine.SetInput(new ControlInput(0.5, 0, 0, 0)); // moderate throttle

        _prevDistance = DistanceToCurrentWaypoint();

        return (GetObservation(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        // Map continuous actions to control inputs
        double throttle = Math.Clamp(action[0], 0, 1);
        double pitch = Math.Clamp(action[1], -1, 1);
        double roll = Math.Clamp(action[2], -1, 1);
        double yaw = Math.Clamp(action[3], -1, 1);

        _engine.SetInput(new ControlInput(throttle, pitch, roll, yaw));
        _engine.Step(_dt);
        _step++;

        // Compute reward
        double reward = ComputeReward(throttle);

        // Check termination
        var state = _engine.State;
        bool crashed = state.Altitude < 0;
        bool allWaypointsReached = _currentWaypoint >= _waypoints.Count;
        bool timeout = _step >= _maxSteps;
        bool done = crashed || allWaypointsReached || timeout;

        if (crashed) reward -= 100.0;

        var info = new Dictionary<string, object>
        {
            { "altitude", state.Altitude },
            { "airspeed", state.Airspeed },
            { "waypoint", _currentWaypoint },
            { "crashed", crashed }
        };

        return (GetObservation(), reward, done, info);
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        // Map discrete to continuous: not primary mode but required by interface
        double throttle = (action & 1) != 0 ? 0.8 : 0.4;
        double pitch = ((action >> 1) & 3) switch { 0 => 0, 1 => 0.3, 2 => -0.3, _ => 0 };
        double roll = ((action >> 3) & 3) switch { 0 => 0, 1 => 0.5, 2 => -0.5, _ => 0 };
        return Step(new VectorN(new[] { throttle, pitch, roll, 0.0 }));
    }

    private VectorN GetObservation()
    {
        var s = _engine.State;
        var euler = s.EulerAngles;

        double dist = DistanceToCurrentWaypoint();
        double headingError = ComputeHeadingError();

        return new VectorN(new[]
        {
            s.Altitude / 1000.0,           // normalized altitude (km)
            s.Airspeed / 100.0,            // normalized airspeed
            s.AngleOfAttack,               // radians
            s.SideslipAngle,               // radians
            euler.roll,                     // radians
            euler.pitch,                    // radians
            euler.yaw,                      // radians
            s.AngularRate.x,               // p (rad/s)
            s.AngularRate.y,               // q (rad/s)
            s.AngularRate.z,               // r (rad/s)
            Math.Min(dist / 1000.0, 10.0), // normalized distance to waypoint
            headingError / Math.PI          // normalized heading error [-1, 1]
        });
    }

    private double ComputeReward(double throttle)
    {
        var s = _engine.State;
        double reward = 0;

        // Distance improvement reward
        double dist = DistanceToCurrentWaypoint();
        double distImprovement = _prevDistance - dist;
        reward += distImprovement * 0.01; // scale down

        // Waypoint reached bonus
        if (dist < _waypointRadius && _currentWaypoint < _waypoints.Count)
        {
            reward += 10.0;
            _currentWaypoint++;
        }

        _prevDistance = _currentWaypoint < _waypoints.Count ? DistanceToCurrentWaypoint() : 0;

        // Fuel efficiency: reward lower throttle slightly
        reward += (1.0 - throttle) * 0.001;

        // Stability penalties
        var euler = s.EulerAngles;
        double bankPenalty = Math.Max(0, Math.Abs(euler.roll) - 1.0); // penalty beyond ~57 degrees
        reward -= bankPenalty * 0.1;

        // Stall penalty (low airspeed)
        if (s.Airspeed < 30) reward -= 0.5;

        // Altitude penalty (too low or too high)
        if (s.Altitude < 100) reward -= 1.0;
        if (s.Altitude > 5000) reward -= 0.1;

        return reward;
    }

    private double DistanceToCurrentWaypoint()
    {
        if (_currentWaypoint >= _waypoints.Count)
            return 0;

        var pos = _engine.State.Position;
        var wp = _waypoints[_currentWaypoint];
        double dx = pos.x - wp.x;
        double dy = pos.y - wp.y;
        double dz = pos.z - wp.z;
        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
    }

    private double ComputeHeadingError()
    {
        if (_currentWaypoint >= _waypoints.Count)
            return 0;

        var pos = _engine.State.Position;
        var wp = _waypoints[_currentWaypoint];

        double desiredHeading = Math.Atan2(wp.y - pos.y, wp.x - pos.x);
        double currentHeading = Math.Atan2(_engine.State.Velocity.y, _engine.State.Velocity.x);

        double error = desiredHeading - currentHeading;
        // Normalize to [-π, π]
        while (error > Math.PI) error -= 2 * Math.PI;
        while (error < -Math.PI) error += 2 * Math.PI;
        return error;
    }

    private static List<Vector> DefaultWaypoints()
    {
        // Simple rectangular circuit at 1000m altitude (NED: z = -1000)
        return new List<Vector>
        {
            new Vector(2000, 0, -1000),
            new Vector(2000, 2000, -1000),
            new Vector(0, 2000, -1000),
            new Vector(0, 0, -1000)
        };
    }
}
