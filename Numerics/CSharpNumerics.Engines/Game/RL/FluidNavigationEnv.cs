using CSharpNumerics.Engines.Game.Fluids;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.RL;

/// <summary>
/// RL environment where an agent navigates through a 2D wind/current field.
///
/// The agent controls a point-mass moving through a fluid simulation.
/// Wind/currents affect the agent's movement, and the goal is to reach
/// a target position while minimizing energy expenditure.
///
/// <b>Observation (8D):</b>
///   [agent_x, agent_y, agent_vx, agent_vy, wind_u, wind_v, target_dx, target_dy]
///   All normalized to grid scale.
///
/// <b>Action (2D continuous):</b>
///   [thrust_x, thrust_y] — engine thrust applied to the agent, each in [-1, 1]
///
/// <b>Reward:</b>
///   - Positive for getting closer to target
///   - Bonus for reaching target
///   - Small penalty per step (encourages efficiency)
///   - Penalty proportional to thrust magnitude (fuel cost)
///
/// <b>Done:</b>
///   Target reached, agent exits grid, or max steps exceeded
/// </summary>
public class FluidNavigationEnv : IEnvironment
{
    public int ObservationSize => 8;
    public int ActionSize => 2;
    public bool IsDiscrete => false;

    private readonly GameFluidSolver2D _solver;
    private readonly FluidConfig _fluidConfig;
    private readonly double _dt;
    private readonly int _maxSteps;
    private readonly double _targetRadius;
    private readonly double _agentMass;
    private readonly double _maxThrust;

    // Agent state
    private double _agentX, _agentY;
    private double _agentVx, _agentVy;
    private double _targetX, _targetY;
    private int _step;
    private double _prevDist;
    private Random _rng;

    /// <summary>
    /// Creates a fluid navigation environment.
    /// </summary>
    /// <param name="gridSize">Grid dimension (square grid).</param>
    /// <param name="dt">Timestep.</param>
    /// <param name="maxSteps">Maximum steps per episode.</param>
    /// <param name="targetRadius">Reach radius for the target.</param>
    /// <param name="agentMass">Agent mass (kg).</param>
    /// <param name="maxThrust">Maximum thrust per axis.</param>
    public FluidNavigationEnv(
        int gridSize = 32,
        double dt = 0.1,
        int maxSteps = 500,
        double targetRadius = 2.0,
        double agentMass = 1.0,
        double maxThrust = 2.0)
    {
        _fluidConfig = new FluidConfig
        {
            GridX = gridSize,
            GridY = gridSize,
            Viscosity = 0.001,
            Quality = FluidQuality.Low
        };
        _solver = new GameFluidSolver2D(_fluidConfig);
        _dt = dt;
        _maxSteps = maxSteps;
        _targetRadius = targetRadius;
        _agentMass = agentMass;
        _maxThrust = maxThrust;
    }

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        _step = 0;

        int gs = _fluidConfig.GridX;

        // Place agent near bottom-left, target near top-right
        _agentX = 3 + _rng.NextDouble() * 3;
        _agentY = 3 + _rng.NextDouble() * 3;
        _agentVx = 0;
        _agentVy = 0;
        _targetX = gs - 5 + _rng.NextDouble() * 2;
        _targetY = gs - 5 + _rng.NextDouble() * 2;

        // Reset fluid solver by creating a fresh one — add a cross-wind emitter
        // to make navigation non-trivial
        _solver.ClearEmitters();
        _solver.ClearObstacles();

        var crossWind = new FluidEmitter(
            new Vector(gs / 2.0, 2, 0),
            densityRate: 0,
            radius: 3.0)
        {
            Velocity = new Vector(3.0, 2.0, 0)
        };
        _solver.AddEmitter(crossWind);

        // Pre-simulate fluid to establish flow
        for (int i = 0; i < 20; i++)
            _solver.Step(0.05);

        _prevDist = DistToTarget();

        return (GetObservation(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        double thrustX = Math.Clamp(action[0], -1, 1) * _maxThrust;
        double thrustY = Math.Clamp(action[1], -1, 1) * _maxThrust;

        // Step fluid simulation
        _solver.Step(_dt);

        // Sample wind at agent position
        var (windU, windV) = _solver.SampleVelocity(_agentX, _agentY);

        // Agent dynamics: F = ma, with wind as external velocity field
        // Effective velocity = agent velocity + wind
        double ax = thrustX / _agentMass;
        double ay = thrustY / _agentMass;

        _agentVx += ax * _dt;
        _agentVy += ay * _dt;

        // Add wind influence (advection by fluid)
        _agentX += (_agentVx + windU) * _dt;
        _agentY += (_agentVy + windV) * _dt;

        // Damping (drag from fluid)
        _agentVx *= 0.98;
        _agentVy *= 0.98;

        _step++;

        // Compute reward
        double dist = DistToTarget();
        double distImprovement = _prevDist - dist;
        double reward = distImprovement * 0.5;

        // Fuel cost
        double thrustMag = Math.Sqrt(thrustX * thrustX + thrustY * thrustY);
        reward -= thrustMag * 0.01;

        // Step penalty
        reward -= 0.01;

        // Reached target
        bool reached = dist < _targetRadius;
        if (reached) reward += 20.0;

        // Out of bounds
        int gs = _fluidConfig.GridX;
        bool oob = _agentX < 1 || _agentX >= gs - 1 || _agentY < 1 || _agentY >= gs - 1;
        if (oob) reward -= 5.0;

        bool done = reached || oob || _step >= _maxSteps;
        _prevDist = dist;

        var info = new Dictionary<string, object>
        {
            { "distance", dist },
            { "reached", reached },
            { "wind_u", windU },
            { "wind_v", windV }
        };

        return (GetObservation(), reward, done, info);
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        // Map discrete to 2D continuous: 4 cardinal + diagonal + none
        double tx = action switch { 0 => 1, 1 => -1, 2 => 0, 3 => 0, _ => 0 };
        double ty = action switch { 0 => 0, 1 => 0, 2 => 1, 3 => -1, _ => 0 };
        return Step(new VectorN(new[] { tx, ty }));
    }

    private VectorN GetObservation()
    {
        int gs = _fluidConfig.GridX;
        var (windU, windV) = _solver.SampleVelocity(_agentX, _agentY);

        return new VectorN(new[]
        {
            _agentX / gs,
            _agentY / gs,
            _agentVx / 5.0,
            _agentVy / 5.0,
            windU / 5.0,
            windV / 5.0,
            (_targetX - _agentX) / gs,
            (_targetY - _agentY) / gs
        });
    }

    private double DistToTarget()
    {
        double dx = _targetX - _agentX;
        double dy = _targetY - _agentY;
        return Math.Sqrt(dx * dx + dy * dy);
    }
}
