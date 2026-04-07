using CSharpNumerics.Numerics.Optimization.Interfaces;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
using System;

namespace CSharpNumerics.ML.NeuralNetwork;

internal static class OptimizerFactory
{
    public static IOptimizer Clone(IOptimizer optimizer)
    {
        return optimizer switch
        {
            Adam adam => new Adam(
                adam.LearningRate,
                adam.Beta1,
                adam.Beta2,
                adam.Epsilon,
                adam.WeightDecay,
                adam.DecoupledWeightDecay),
            GradientDescent gradientDescent => new GradientDescent(
                gradientDescent.LearningRate,
                gradientDescent.Momentum,
                gradientDescent.L2,
                gradientDescent.Nesterov),
            _ => throw new NotSupportedException($"Optimizer type {optimizer.GetType().Name} is not supported for cloned layer state.")
        };
    }
}