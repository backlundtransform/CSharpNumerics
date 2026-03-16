using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Environments;

/// <summary>
/// Classic Pendulum (inverted pendulum swing-up) environment with continuous control.
///
/// The goal is to swing up a pendulum and balance it upright by applying torque.
///
/// State: [cos(θ), sin(θ), θ_dot]  (3 continuous values)
/// Action: torque ∈ [-MaxTorque, MaxTorque]  (1 continuous value)
/// Reward: -(θ² + 0.1 θ_dot² + 0.001 torque²)  (dense, maximized when upright and still)
/// Done: never (fixed-length episodes, typically 200 steps)
///
/// Based on OpenAI Gym Pendulum-v1.
/// </summary>
public class Pendulum : IEnvironment
{
    private const double G = 10.0;
    private const double Mass = 1.0;
    private const double Length = 1.0;
    private const double Dt = 0.05;

    public double MaxTorque { get; set; } = 2.0;
    public double MaxSpeed { get; set; } = 8.0;
    public int MaxSteps { get; set; } = 200;

    private double _theta;
    private double _thetaDot;
    private int _step;
    private Random _rng;

    public int ObservationSize => 3;
    public int ActionSize => 1;
    public bool IsDiscrete => false;

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        // Random initial angle and angular velocity
        _theta = Uniform(-Math.PI, Math.PI);
        _thetaDot = Uniform(-1.0, 1.0);
        _step = 0;
        return (GetState(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
    {
        double torque = Math.Clamp(action[0], -MaxTorque, MaxTorque);

        double thetaAcc = (-3.0 * G / (2.0 * Length)) * Math.Sin(_theta + Math.PI)
                        + (3.0 / (Mass * Length * Length)) * torque;

        _thetaDot += thetaAcc * Dt;
        _thetaDot = Math.Clamp(_thetaDot, -MaxSpeed, MaxSpeed);
        _theta += _thetaDot * Dt;

        // Normalize angle to [-π, π]
        _theta = NormalizeAngle(_theta);
        _step++;

        // Reward: penalize angle from upright, angular velocity, and torque magnitude
        double reward = -(AngleCost(_theta) + 0.1 * _thetaDot * _thetaDot + 0.001 * torque * torque);

        bool done = _step >= MaxSteps;
        return (GetState(), reward, done, new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        // Map discrete action to continuous: 0 = -MaxTorque, 1 = 0, 2 = +MaxTorque
        double torque = (action - 1) * MaxTorque;
        return Step(new VectorN(new[] { torque }));
    }

    private VectorN GetState() => new VectorN(new[]
    {
        Math.Cos(_theta),
        Math.Sin(_theta),
        _thetaDot
    });

    private static double NormalizeAngle(double angle)
    {
        // Wrap to [-π, π]
        angle %= (2 * Math.PI);
        if (angle > Math.PI) angle -= 2 * Math.PI;
        if (angle < -Math.PI) angle += 2 * Math.PI;
        return angle;
    }

    private static double AngleCost(double theta)
    {
        // θ² where θ is the angle from upright (0 = upright)
        return theta * theta;
    }

    private double Uniform(double lo, double hi) => lo + _rng.NextDouble() * (hi - lo);
}
