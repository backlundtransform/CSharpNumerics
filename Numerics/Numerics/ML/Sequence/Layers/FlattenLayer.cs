using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

public class FlattenLayer : ILayer
{
    private int _timesteps;
    private int _channels;

    public int ParameterCount => 0;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        _timesteps = input.Length;
        _channels = input[0].Length;
        var values = new double[_timesteps * _channels];
        int index = 0;

        for (int timestep = 0; timestep < _timesteps; timestep++)
        {
            if (input[timestep].Length != _channels)
                throw new ArgumentException("All timesteps must have the same channel count.", nameof(input));

            for (int channel = 0; channel < _channels; channel++)
                values[index++] = input[timestep][channel];
        }

        return new[] { new VectorN(values) };
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_timesteps == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != 1)
            throw new ArgumentException("Flatten expects a single gradient vector.", nameof(gradOutput));
        if (gradOutput[0].Length != _timesteps * _channels)
            throw new ArgumentException("Gradient width must match the flattened output width.", nameof(gradOutput));

        var result = new VectorN[_timesteps];
        int index = 0;
        for (int timestep = 0; timestep < _timesteps; timestep++)
        {
            var values = new double[_channels];
            for (int channel = 0; channel < _channels; channel++)
                values[channel] = gradOutput[0][index++];
            result[timestep] = new VectorN(values);
        }

        return result;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
    }
}