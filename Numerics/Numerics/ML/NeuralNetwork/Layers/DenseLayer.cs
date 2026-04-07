using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.NeuralNetwork.Layers;

public class DenseLayer : ILayer
{
    private Matrix _weights;
    private VectorN _biases;
    private Matrix _weightGradients;
    private VectorN _biasGradients;
    private VectorN[] _lastInput = Array.Empty<VectorN>();
    private VectorN[] _lastActivatedOutput = Array.Empty<VectorN>();
    private IOptimizer _weightOptimizerState;
    private IOptimizer _biasOptimizerState;

    public DenseLayer(int inputSize, int outputSize, ActivationType activation = ActivationType.Linear, int seed = 123)
        : this(CreateXavierMatrix(inputSize, outputSize, seed), new VectorN(outputSize), activation)
    {
    }

    public DenseLayer(Matrix weights, VectorN biases, ActivationType activation = ActivationType.Linear)
    {
        if (weights.rowLength <= 0 || weights.columnLength <= 0)
        {
            throw new ArgumentException("Weights must have a positive shape.", nameof(weights));
        }

        if (weights.columnLength != biases.Length)
        {
            throw new ArgumentException("Bias length must match the weight output dimension.", nameof(biases));
        }

        _weights = new Matrix((double[,])weights.values.Clone());
        _biases = new VectorN(biases.Values);
        _weightGradients = new Matrix(weights.rowLength, weights.columnLength);
        _biasGradients = new VectorN(biases.Values.Length);
        Activation = activation;
    }

    public ActivationType Activation { get; }

    public int ParameterCount => (_weights.rowLength * _weights.columnLength) + _biases.Length;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
        {
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));
        }

        var output = new VectorN[input.Length];

        for (int i = 0; i < input.Length; i++)
        {
            ValidateInputWidth(input[i]);
            var z = (_weights.Transpose() * input[i]) + _biases;
            output[i] = Activations.Apply(z, Activation);
        }

        _lastInput = CloneVectors(input);
        _lastActivatedOutput = CloneVectors(output);
        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_lastInput.Length == 0)
        {
            throw new InvalidOperationException("Call Forward before Backward.");
        }

        if (gradOutput is null || gradOutput.Length != _lastInput.Length)
        {
            throw new ArgumentException("Gradient sequence must match the cached input sequence length.", nameof(gradOutput));
        }

        var gradInput = new VectorN[gradOutput.Length];

        for (int i = 0; i < gradOutput.Length; i++)
        {
            if (gradOutput[i].Length != _biases.Length)
            {
                throw new ArgumentException("Gradient width must match the layer output dimension.", nameof(gradOutput));
            }

            var localGradient = gradOutput[i].Hadamard(Activations.Derivative(_lastActivatedOutput[i], Activation));
            _weightGradients += _lastInput[i].Outer(localGradient);
            _biasGradients += localGradient;
            gradInput[i] = _weights * localGradient;
        }

        return gradInput;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        if (batchSize <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(batchSize), "Batch size must be positive.");
        }

        _weightOptimizerState ??= OptimizerFactory.Clone(weightOptimizer);
        _biasOptimizerState ??= OptimizerFactory.Clone(biasOptimizer);

        int rows = _weights.rowLength;
        int cols = _weights.columnLength;
        var flattenedWeights = new double[rows * cols];
        var flattenedGradients = new double[rows * cols];

        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                int index = (row * cols) + col;
                flattenedWeights[index] = _weights.values[row, col];
                flattenedGradients[index] = _weightGradients.values[row, col] / batchSize;
            }
        }

        var updatedWeights = _weightOptimizerState.Step(new VectorN(flattenedWeights), new VectorN(flattenedGradients));
        var updatedMatrix = new double[rows, cols];
        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                updatedMatrix[row, col] = updatedWeights[(row * cols) + col];
            }
        }

        _weights = new Matrix(updatedMatrix);
        _biases = _biasOptimizerState.Step(_biases, _biasGradients / batchSize);

        _weightGradients = new Matrix(rows, cols);
        _biasGradients = new VectorN(_biases.Length);
    }

    private static Matrix CreateXavierMatrix(int inputSize, int outputSize, int seed)
    {
        if (inputSize <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(inputSize), "Input size must be positive.");
        }

        if (outputSize <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(outputSize), "Output size must be positive.");
        }

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
        {
            clone[i] = new VectorN(vectors[i].Values);
        }

        return clone;
    }

    private void ValidateInputWidth(VectorN input)
    {
        if (input.Length != _weights.rowLength)
        {
            throw new ArgumentException("Input width must match the layer input dimension.", nameof(input));
        }
    }
}