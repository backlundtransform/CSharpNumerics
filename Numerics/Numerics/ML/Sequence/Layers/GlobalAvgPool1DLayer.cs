using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

public class GlobalAvgPool1DLayer : ILayer
{
    private int _inputLength;
    private int _channels;

    public int ParameterCount => 0;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        _inputLength = input.Length;
        _channels = input[0].Length;
        var values = new double[_channels];

        for (int timestep = 0; timestep < input.Length; timestep++)
        {
            if (input[timestep].Length != _channels)
                throw new ArgumentException("All timesteps must have the same channel count.", nameof(input));

            for (int channel = 0; channel < _channels; channel++)
                values[channel] += input[timestep][channel];
        }

        for (int channel = 0; channel < _channels; channel++)
            values[channel] /= _inputLength;

        return new[] { new VectorN(values) };
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_inputLength == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != 1)
            throw new ArgumentException("Global average pooling expects a single pooled gradient vector.", nameof(gradOutput));
        if (gradOutput[0].Length != _channels)
            throw new ArgumentException("Gradient width must match the pooled channel count.", nameof(gradOutput));

        var result = new VectorN[_inputLength];
        for (int timestep = 0; timestep < _inputLength; timestep++)
        {
            var values = new double[_channels];
            for (int channel = 0; channel < _channels; channel++)
                values[channel] = gradOutput[0][channel] / _inputLength;
            result[timestep] = new VectorN(values);
        }

        return result;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
    }
}