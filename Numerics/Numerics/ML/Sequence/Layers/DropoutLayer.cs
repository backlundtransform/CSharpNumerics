using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// Inverted dropout regulariser. During training each element is independently zeroed
/// with probability <c>rate</c> and the survivors are scaled by <c>1/(1-rate)</c> so the
/// expected activation is unchanged. At inference the layer is a pass-through, so no
/// rescaling is needed downstream. Parameter-free.
/// </summary>
public class DropoutLayer : ILayer
{
    private readonly double _rate;
    private readonly Random _random;
    private bool[][] _mask = Array.Empty<bool[]>();

    public DropoutLayer(double rate, int seed = 123)
    {
        if (rate < 0.0 || rate >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(rate), "Dropout rate must be in [0, 1).");

        _rate = rate;
        _random = new Random(seed);
    }

    public int ParameterCount => 0;

    /// <summary>Dropout probability.</summary>
    public double Rate => _rate;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        if (!training || _rate == 0.0)
        {
            _mask = Array.Empty<bool[]>();
            return CloneVectors(input);
        }

        double scale = 1.0 / (1.0 - _rate);
        var output = new VectorN[input.Length];
        _mask = new bool[input.Length][];

        for (int t = 0; t < input.Length; t++)
        {
            int width = input[t].Length;
            var keep = new bool[width];
            var values = new double[width];

            for (int c = 0; c < width; c++)
            {
                bool kept = _random.NextDouble() >= _rate;
                keep[c] = kept;
                values[c] = kept ? input[t][c] * scale : 0.0;
            }

            _mask[t] = keep;
            output[t] = new VectorN(values);
        }

        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (gradOutput is null || gradOutput.Length == 0)
            throw new ArgumentException("Gradient sequence must contain at least one timestep.", nameof(gradOutput));

        // Inference / rate 0: identity gradient.
        if (_mask.Length == 0)
            return CloneVectors(gradOutput);

        if (gradOutput.Length != _mask.Length)
            throw new ArgumentException("Gradient sequence must match the cached input length.", nameof(gradOutput));

        double scale = 1.0 / (1.0 - _rate);
        var gradInput = new VectorN[gradOutput.Length];

        for (int t = 0; t < gradOutput.Length; t++)
        {
            int width = gradOutput[t].Length;
            var values = new double[width];
            for (int c = 0; c < width; c++)
                values[c] = _mask[t][c] ? gradOutput[t][c] * scale : 0.0;
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
