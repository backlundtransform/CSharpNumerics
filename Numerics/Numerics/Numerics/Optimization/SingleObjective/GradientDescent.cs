using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.Numerics.Optimization.SingleObjective;

/// <summary>
/// Gradient descent with optional momentum and L2 regularisation.
/// Supports vanilla GD, momentum GD, and Nesterov accelerated gradient.
/// </summary>
public class GradientDescent : IOptimizer
{
    public double LearningRate { get; set; }
    public double Momentum { get; set; }
    public double L2 { get; set; }
    public bool Nesterov { get; set; }

    private VectorN _velocity;
    private bool _initialised;

    public GradientDescent(double learningRate = 0.01, double momentum = 0.0,
        double l2 = 0.0, bool nesterov = false)
    {
        LearningRate = learningRate;
        Momentum = momentum;
        L2 = l2;
        Nesterov = nesterov;
    }

    public VectorN Step(VectorN parameters, VectorN gradient)
    {
        if (!_initialised || _velocity.Length != parameters.Length)
        {
            _velocity = new VectorN(parameters.Length);
            _initialised = true;
        }

        // L2 regularisation: add λ·w to gradient
        VectorN g = L2 > 0
            ? gradient + L2 * parameters
            : gradient;

        if (Momentum > 0)
        {
            _velocity = Momentum * _velocity + LearningRate * g;

            if (Nesterov)
                return parameters - (Momentum * _velocity + LearningRate * g);

            return parameters - _velocity;
        }

        return parameters - LearningRate * g;
    }

    public void Reset()
    {
        _initialised = false;
        _velocity = default;
    }
}
