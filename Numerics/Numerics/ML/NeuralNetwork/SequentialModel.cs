using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;
using System.Linq;

namespace CSharpNumerics.ML.NeuralNetwork;

public class SequentialModel
{
    private readonly ILayer[] _layers;

    public SequentialModel(params ILayer[] layers)
    {
        if (layers is null || layers.Length == 0)
        {
            throw new ArgumentException("SequentialModel requires at least one layer.", nameof(layers));
        }

        if (layers.Any(layer => layer is null))
        {
            throw new ArgumentException("SequentialModel cannot contain null layers.", nameof(layers));
        }

        _layers = (ILayer[])layers.Clone();
    }

    public int ParameterCount => _layers.Sum(layer => layer.ParameterCount);

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        var activations = input;

        foreach (var layer in _layers)
        {
            activations = layer.Forward(activations, training);
        }

        return activations;
    }

    public VectorN ForwardSingle(VectorN[] input, bool training = true)
    {
        var output = Forward(input, training);
        if (output.Length != 1)
        {
            throw new InvalidOperationException("ForwardSingle requires the final layer to emit exactly one timestep.");
        }

        return output[0];
    }

    public VectorN[] Backward(VectorN[] lossGradient)
    {
        var gradients = lossGradient;

        for (int i = _layers.Length - 1; i >= 0; i--)
        {
            gradients = _layers[i].Backward(gradients);
        }

        return gradients;
    }

    public VectorN BackwardSingle(VectorN lossGradient)
    {
        var gradients = Backward(new[] { lossGradient });
        if (gradients.Length != 1)
        {
            throw new InvalidOperationException("BackwardSingle requires the first layer to consume exactly one timestep.");
        }

        return gradients[0];
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        foreach (var layer in _layers)
        {
            layer.ApplyGradients(weightOptimizer, biasOptimizer, batchSize);
        }
    }
}