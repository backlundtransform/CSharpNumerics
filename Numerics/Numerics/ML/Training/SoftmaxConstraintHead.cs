using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Training;

/// <summary>
/// Output head that maps raw scores to a constrained ratio via softmax: every output lies in
/// [0, 1] and the outputs sum to 1. Useful for predicting the fraction of total flow assigned to
/// each component, guaranteeing the partition is valid by construction. Parameter-free.
/// </summary>
public class SoftmaxConstraintHead : ILayer
{
    private VectorN[] _lastOutput = Array.Empty<VectorN>();

    public int ParameterCount => 0;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        var output = new VectorN[input.Length];
        for (int t = 0; t < input.Length; t++)
            output[t] = Activations.Softmax(input[t]);

        _lastOutput = CloneVectors(output);
        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_lastOutput.Length == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != _lastOutput.Length)
            throw new ArgumentException("Gradient sequence must match the output sequence length.", nameof(gradOutput));

        var gradInput = new VectorN[gradOutput.Length];
        for (int t = 0; t < gradOutput.Length; t++)
        {
            var s = _lastOutput[t];
            var g = gradOutput[t];

            // Softmax Jacobian-vector product: dz_i = s_i (g_i - Σ_j g_j s_j).
            double dot = 0.0;
            for (int j = 0; j < s.Length; j++)
                dot += g[j] * s[j];

            var values = new double[s.Length];
            for (int i = 0; i < s.Length; i++)
                values[i] = s[i] * (g[i] - dot);

            gradInput[t] = new VectorN(values);
        }

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
