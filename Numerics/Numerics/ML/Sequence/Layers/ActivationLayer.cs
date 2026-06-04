using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// Applies a pointwise activation function to every element of a sequence.
/// A parameter-free <see cref="ILayer"/> that lets activations be composed
/// between other layers (e.g. Conv → BatchNorm → ReLU → Dropout in a TCN block).
/// </summary>
public class ActivationLayer : ILayer
{
    private readonly ActivationType _activation;
    private VectorN[] _lastActivatedOutput = Array.Empty<VectorN>();

    public ActivationLayer(ActivationType activation)
    {
        _activation = activation;
    }

    public int ParameterCount => 0;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        var output = new VectorN[input.Length];
        for (int t = 0; t < input.Length; t++)
            output[t] = Activations.Apply(input[t], _activation);

        _lastActivatedOutput = CloneVectors(output);
        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_lastActivatedOutput.Length == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != _lastActivatedOutput.Length)
            throw new ArgumentException("Gradient sequence must match the output sequence length.", nameof(gradOutput));

        var gradInput = new VectorN[gradOutput.Length];
        for (int t = 0; t < gradOutput.Length; t++)
            gradInput[t] = gradOutput[t].Hadamard(Activations.Derivative(_lastActivatedOutput[t], _activation));

        return gradInput;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
    }

    private static VectorN[] CloneVectors(VectorN[] vectors)
    {
        var clone = new VectorN[vectors.Length];
        for (int i = 0; i < vectors.Length; i++)
            clone[i] = new VectorN(vectors[i].Values);
        return clone;
    }
}
