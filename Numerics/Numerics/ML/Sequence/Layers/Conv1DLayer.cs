using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.ML.Sequence.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

public class Conv1DLayer : ILayer
{
    private readonly int _inputChannels;
    private readonly int _filters;
    private readonly int _kernelSize;
    private readonly int _stride;
    private readonly ConvolutionPaddingMode _padding;
    private readonly ActivationType _activation;

    private Matrix _weights;
    private VectorN _biases;
    private Matrix _weightGradients;
    private VectorN _biasGradients;
    private VectorN[] _lastPaddedInput = Array.Empty<VectorN>();
    private VectorN[] _lastActivatedOutput = Array.Empty<VectorN>();
    private VectorN[] _lastWindows = Array.Empty<VectorN>();
    private int _lastInputLength;
    private int _padLeft;
    private IOptimizer _weightOptimizerState;
    private IOptimizer _biasOptimizerState;

    public Conv1DLayer(
        int inputChannels,
        int filters,
        int kernelSize,
        int stride = 1,
        ConvolutionPaddingMode padding = ConvolutionPaddingMode.Same,
        ActivationType activation = ActivationType.ReLU,
        int seed = 123)
    {
        if (inputChannels <= 0)
            throw new ArgumentOutOfRangeException(nameof(inputChannels), "Input channels must be positive.");
        if (filters <= 0)
            throw new ArgumentOutOfRangeException(nameof(filters), "Filter count must be positive.");
        if (kernelSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(kernelSize), "Kernel size must be positive.");
        if (stride <= 0)
            throw new ArgumentOutOfRangeException(nameof(stride), "Stride must be positive.");

        _inputChannels = inputChannels;
        _filters = filters;
        _kernelSize = kernelSize;
        _stride = stride;
        _padding = padding;
        _activation = activation;

        _weights = CreateXavierMatrix(kernelSize * inputChannels, filters, seed);
        _biases = new VectorN(filters);
        _weightGradients = new Matrix(_weights.rowLength, _weights.columnLength);
        _biasGradients = new VectorN(filters);
    }

    public Conv1DLayer(
        Matrix weights,
        VectorN biases,
        int inputChannels,
        int kernelSize,
        int stride = 1,
        ConvolutionPaddingMode padding = ConvolutionPaddingMode.Same,
        ActivationType activation = ActivationType.ReLU)
    {
        if (inputChannels <= 0)
            throw new ArgumentOutOfRangeException(nameof(inputChannels), "Input channels must be positive.");
        if (kernelSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(kernelSize), "Kernel size must be positive.");
        if (stride <= 0)
            throw new ArgumentOutOfRangeException(nameof(stride), "Stride must be positive.");
        if (weights.rowLength != kernelSize * inputChannels)
            throw new ArgumentException("Weight rows must equal kernelSize * inputChannels.", nameof(weights));
        if (weights.columnLength != biases.Length)
            throw new ArgumentException("Bias length must match the filter count.", nameof(biases));

        _inputChannels = inputChannels;
        _filters = biases.Length;
        _kernelSize = kernelSize;
        _stride = stride;
        _padding = padding;
        _activation = activation;
        _weights = new Matrix((double[,])weights.values.Clone());
        _biases = new VectorN(biases.Values);
        _weightGradients = new Matrix(_weights.rowLength, _weights.columnLength);
        _biasGradients = new VectorN(_filters);
    }

    public int ParameterCount => (_weights.rowLength * _weights.columnLength) + _biases.Length;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        ValidateInput(input);

        _lastInputLength = input.Length;
        var (padLeft, padRight, outputLength) = ComputePadding(input.Length);
        _padLeft = padLeft;

        _lastPaddedInput = PadInput(input, padLeft, padRight);
        _lastWindows = new VectorN[outputLength];
        _lastActivatedOutput = new VectorN[outputLength];

        for (int outputIndex = 0; outputIndex < outputLength; outputIndex++)
        {
            int start = outputIndex * _stride;
            var window = FlattenWindow(_lastPaddedInput, start);
            var z = (_weights.Transpose() * window) + _biases;
            _lastWindows[outputIndex] = window;
            _lastActivatedOutput[outputIndex] = Activations.Apply(z, _activation);
        }

        return CloneVectors(_lastActivatedOutput);
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_lastActivatedOutput.Length == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != _lastActivatedOutput.Length)
            throw new ArgumentException("Gradient sequence must match the output sequence length.", nameof(gradOutput));

        var gradPadded = new double[_lastPaddedInput.Length, _inputChannels];

        for (int outputIndex = 0; outputIndex < gradOutput.Length; outputIndex++)
        {
            if (gradOutput[outputIndex].Length != _filters)
                throw new ArgumentException("Gradient width must match the filter count.", nameof(gradOutput));

            var localGradient = gradOutput[outputIndex].Hadamard(
                Activations.Derivative(_lastActivatedOutput[outputIndex], _activation));

            _weightGradients += _lastWindows[outputIndex].Outer(localGradient);
            _biasGradients += localGradient;

            var gradWindow = _weights * localGradient;
            int start = outputIndex * _stride;
            int index = 0;
            for (int kernelOffset = 0; kernelOffset < _kernelSize; kernelOffset++)
            {
                for (int channel = 0; channel < _inputChannels; channel++)
                {
                    gradPadded[start + kernelOffset, channel] += gradWindow[index++];
                }
            }
        }

        var gradInput = new VectorN[_lastInputLength];
        for (int timestep = 0; timestep < _lastInputLength; timestep++)
        {
            var values = new double[_inputChannels];
            for (int channel = 0; channel < _inputChannels; channel++)
            {
                values[channel] = gradPadded[timestep + _padLeft, channel];
            }

            gradInput[timestep] = new VectorN(values);
        }

        return gradInput;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        if (batchSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(batchSize), "Batch size must be positive.");

        _weightOptimizerState ??= OptimizerFactory.Clone(weightOptimizer);
        _biasOptimizerState ??= OptimizerFactory.Clone(biasOptimizer);

        int rows = _weights.rowLength;
        int cols = _weights.columnLength;
        var flatWeights = new double[rows * cols];
        var flatGradients = new double[rows * cols];

        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                int index = (row * cols) + col;
                flatWeights[index] = _weights.values[row, col];
                flatGradients[index] = _weightGradients.values[row, col] / batchSize;
            }
        }

        var updatedWeights = _weightOptimizerState.Step(new VectorN(flatWeights), new VectorN(flatGradients));
        var values = new double[rows, cols];
        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                values[row, col] = updatedWeights[(row * cols) + col];
            }
        }

        _weights = new Matrix(values);
        _biases = _biasOptimizerState.Step(_biases, _biasGradients / batchSize);
        _weightGradients = new Matrix(rows, cols);
        _biasGradients = new VectorN(_filters);
    }

    private (int padLeft, int padRight, int outputLength) ComputePadding(int inputLength)
    {
        if (_padding == ConvolutionPaddingMode.Valid)
        {
            int numerator = inputLength - _kernelSize;
            if (numerator < 0)
                throw new InvalidOperationException("Kernel size cannot exceed sequence length when using valid padding.");

            return (0, 0, (numerator / _stride) + 1);
        }

        int outputLength = (inputLength + _stride - 1) / _stride;
        int totalPadding = Math.Max(((outputLength - 1) * _stride) + _kernelSize - inputLength, 0);
        int padLeft = totalPadding / 2;
        return (padLeft, totalPadding - padLeft, outputLength);
    }

    private VectorN[] PadInput(VectorN[] input, int padLeft, int padRight)
    {
        var padded = new VectorN[input.Length + padLeft + padRight];

        for (int i = 0; i < padLeft; i++)
            padded[i] = new VectorN(_inputChannels);

        for (int i = 0; i < input.Length; i++)
            padded[i + padLeft] = new VectorN(input[i].Values);

        for (int i = 0; i < padRight; i++)
            padded[input.Length + padLeft + i] = new VectorN(_inputChannels);

        return padded;
    }

    private VectorN FlattenWindow(VectorN[] input, int start)
    {
        var values = new double[_kernelSize * _inputChannels];
        int index = 0;

        for (int kernelOffset = 0; kernelOffset < _kernelSize; kernelOffset++)
        {
            var timestep = input[start + kernelOffset];
            for (int channel = 0; channel < _inputChannels; channel++)
            {
                values[index++] = timestep[channel];
            }
        }

        return new VectorN(values);
    }

    private void ValidateInput(VectorN[] input)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        for (int i = 0; i < input.Length; i++)
        {
            if (input[i].Length != _inputChannels)
                throw new ArgumentException("Each timestep must match the configured input channel count.", nameof(input));
        }
    }

    private static Matrix CreateXavierMatrix(int inputSize, int outputSize, int seed)
    {
        var random = new Random(seed);
        var values = new double[inputSize, outputSize];
        double limit = Math.Sqrt(6.0 / (inputSize + outputSize));

        for (int row = 0; row < inputSize; row++)
        {
            for (int col = 0; col < outputSize; col++)
            {
                values[row, col] = (random.NextDouble() * 2.0 - 1.0) * limit;
            }
        }

        return new Matrix(values);
    }

    private static VectorN[] CloneVectors(VectorN[] vectors)
    {
        var clone = new VectorN[vectors.Length];
        for (int i = 0; i < vectors.Length; i++)
            clone[i] = new VectorN(vectors[i].Values);
        return clone;
    }
}