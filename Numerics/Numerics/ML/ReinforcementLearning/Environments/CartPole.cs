using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Environments;

/// <summary>
/// Classic CartPole (inverted pendulum) environment.
/// A pole is attached to a cart on a frictionless track.
/// The agent applies a force of +/- 10 N to the cart.
///
/// State: [x, x_dot, theta, theta_dot]  (4 continuous values)
/// Actions: 0 = push left, 1 = push right
/// Reward: +1 for every step the pole stays upright
/// Done: |x| > 2.4  or  |theta| > 12°  or  steps > 500
///
/// Physics based on Barto, Sutton & Anderson (1983).
/// </summary>
public class CartPole : IEnvironment
{
    private const double Gravity = 9.8;
    private const double CartMass = 1.0;
    private const double PoleMass = 0.1;
    private const double TotalMass = CartMass + PoleMass;
    private const double HalfPoleLength = 0.5;
    private const double ForceMagnitude = 10.0;
    private const double Tau = 0.02; // time step
    private const double ThetaThreshold = 12.0 * Math.PI / 180.0; // 12 degrees
    private const double XThreshold = 2.4;

    private double _x, _xDot, _theta, _thetaDot;
    private int _step;
    private Random _rng;

    public int ObservationSize => 4;
    public int ActionSize => 2;
    public bool IsDiscrete => true;

    public (VectorN state, Dictionary<string, object> info) Reset(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
        // Small random initial state
        _x = Uniform(-0.05, 0.05);
        _xDot = Uniform(-0.05, 0.05);
        _theta = Uniform(-0.05, 0.05);
        _thetaDot = Uniform(-0.05, 0.05);
        _step = 0;
        return (GetState(), new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action)
    {
        double force = action == 1 ? ForceMagnitude : -ForceMagnitude;

        double cosTheta = Math.Cos(_theta);
        double sinTheta = Math.Sin(_theta);

        double temp = (force + PoleMass * HalfPoleLength * _thetaDot * _thetaDot * sinTheta) / TotalMass;
        double thetaAcc = (Gravity * sinTheta - cosTheta * temp)
            / (HalfPoleLength * (4.0 / 3.0 - PoleMass * cosTheta * cosTheta / TotalMass));
        double xAcc = temp - PoleMass * HalfPoleLength * thetaAcc * cosTheta / TotalMass;

        // Euler integration
        _x += Tau * _xDot;
        _xDot += Tau * xAcc;
        _theta += Tau * _thetaDot;
        _thetaDot += Tau * thetaAcc;
        _step++;

        bool done = Math.Abs(_x) > XThreshold
                 || Math.Abs(_theta) > ThetaThreshold
                 || _step >= 500;

        double reward = done && _step < 500 ? 0.0 : 1.0;

        return (GetState(), reward, done, new Dictionary<string, object>());
    }

    public (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action)
        => Step((int)action[0]);

    private VectorN GetState() => new VectorN(new[] { _x, _xDot, _theta, _thetaDot });

    private double Uniform(double lo, double hi) => lo + _rng.NextDouble() * (hi - lo);
}
