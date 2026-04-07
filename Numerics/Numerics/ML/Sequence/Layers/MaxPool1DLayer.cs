using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

public class MaxPool1DLayer : ILayer
{
    private readonly int _poolSize;
    private readonly int _stride;
    private int[,] _argmaxIndices = new int[0, 0];
    private int _inputLength;
    private int _channels;

    public MaxPool1DLayer(int poolSize = 2, int stride = 2)
    {
        if (poolSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(poolSize), "Pool size must be positive.");
        if (stride <= 0)
            throw new ArgumentOutOfRangeException(nameof(stride), "Stride must be positive.");

        _poolSize = poolSize;
        _stride = stride;
    }

    public int ParameterCount => 0;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        _inputLength = input.Length;
        _channels = input[0].Length;
        if (_inputLength < _poolSize)
            throw new InvalidOperationException("Pool size cannot exceed the sequence length.");

        int outputLength = ((_inputLength - _poolSize) / _stride) + 1;
        var output = new VectorN[outputLength];
        _argmaxIndices = new int[outputLength, _channels];

        for (int outputIndex = 0; outputIndex < outputLength; outputIndex++)
        {
            int start = outputIndex * _stride;
            var pooled = new double[_channels];

            for (int channel = 0; channel < _channels; channel++)
            {
                double maxValue = double.NegativeInfinity;
                int maxIndex = start;

                for (int offset = 0; offset < _poolSize; offset++)
                {
                    int index = start + offset;
                    if (input[index][channel] > maxValue)
                    {
                        maxValue = input[index][channel];
                        maxIndex = index;
                    }
                }

                pooled[channel] = maxValue;
                _argmaxIndices[outputIndex, channel] = maxIndex;
            }

            output[outputIndex] = new VectorN(pooled);
        }

        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_inputLength == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != _argmaxIndices.GetLength(0))
            throw new ArgumentException("Gradient sequence must match the pooled output length.", nameof(gradOutput));

        var gradients = new double[_inputLength, _channels];

        for (int outputIndex = 0; outputIndex < gradOutput.Length; outputIndex++)
        {
            if (gradOutput[outputIndex].Length != _channels)
                throw new ArgumentException("Gradient width must match the pooled channel count.", nameof(gradOutput));

            for (int channel = 0; channel < _channels; channel++)
            {
                gradients[_argmaxIndices[outputIndex, channel], channel] += gradOutput[outputIndex][channel];
            }
        }

        var result = new VectorN[_inputLength];
        for (int timestep = 0; timestep < _inputLength; timestep++)
        {
            var values = new double[_channels];
            for (int channel = 0; channel < _channels; channel++)
                values[channel] = gradients[timestep, channel];
            result[timestep] = new VectorN(values);
        }

        return result;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
    }
}