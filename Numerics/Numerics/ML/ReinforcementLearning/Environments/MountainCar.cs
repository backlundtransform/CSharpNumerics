using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Environments;

/// <summary>
/// Classic MountainCar environment.
/// A car is stuck in a valley and must build momentum to reach the goal on the right.
///
/// State: [position, velocity]  (2 continuous values)
///   position ∈ [-1.2, 0.6],  velocity ∈ [-0.07, 0.07]
/// Actions: 0 = push left, 1 = no push, 2 = push right
/// Reward: -1 per step (sparse reward — only terminates when goal reached or max steps)
/// Done: position >= 0.5  or  steps >= 200
///
/// Based on Moore (1990).
/// </summary>
public class MountainCar : IEnvironment
{
    private const double MinPosition = -1.2;
    private const double MaxPosition = 0.6;
    private const double MaxSpeed = 0.07;
    private const double GoalPosition = 0.5;
    private const double Force = 0.001;
    private const double GravityFactor = 0.0025;

    private double _position;
    private double _velocity;
    private int _step;
    private Random _rng;

    public int ObservationSize => 2;
    public int ActionSize => 3;
    public bool IsDiscrete => true;
    public int MaxSteps { get; set; } = 200;

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        _position = -0.6 + _rng.NextDouble() * 0.2; // [-0.6, -0.4]
        _velocity = 0.0;
        _step = 0;
        return (GetState(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        _velocity += (action - 1) * Force - Math.Cos(3 * _position) * GravityFactor;
        _velocity = Math.Clamp(_velocity, -MaxSpeed, MaxSpeed);

        _position += _velocity;
        _position = Math.Clamp(_position, MinPosition, MaxPosition);

        // Bounce off left wall
        if (_position <= MinPosition && _velocity < 0)
            _velocity = 0;

        _step++;
        bool done = _position >= GoalPosition || _step >= MaxSteps;
        double reward = _position >= GoalPosition ? 0.0 : -1.0;

        return (GetState(), reward, done, new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
        => Step((int)action[0]);

    private VectorN GetState() => new VectorN(new[] { _position, _velocity });
}
