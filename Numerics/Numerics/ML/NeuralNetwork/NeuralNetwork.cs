using CSharpNumerics.ML.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.NeuralNetwork;

/// <summary>
/// A feedforward neural network with configurable hidden layers, activation functions,
/// and output behaviour. Shared foundation for MLPClassifier, MLPRegressor, and RL agents.
/// </summary>
public class NeuralNetwork
{
    private List<Matrix> _weights;
    private List<VectorN> _biases;
    private List<VectorN> _lastActivations;

    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public OutputMode Output { get; set; } = OutputMode.Linear;
    public int InputSize { get; private set; }
    public int OutputSize { get; private set; }
    public int[] HiddenLayers { get; private set; }
    public bool IsInitialized => _weights != null;

    public enum OutputMode { Linear, Softmax }

    public NeuralNetwork(int inputSize, int[] hiddenLayers, int outputSize,
        ActivationType activation = ActivationType.ReLU,
        OutputMode output = OutputMode.Linear,
        int seed = 123)
    {
        InputSize = inputSize;
        OutputSize = outputSize;
        HiddenLayers = (int[])hiddenLayers.Clone();
        Activation = activation;
        Output = output;

        Initialize(seed);
    }

    private NeuralNetwork() { }

    private void Initialize(int seed)
    {
        _weights = new List<Matrix>();
        _biases = new List<VectorN>();

        var layers = new List<int> { InputSize };
        layers.AddRange(HiddenLayers);
        layers.Add(OutputSize);

        var rnd = new Random(seed);
        for (int i = 0; i < layers.Count - 1; i++)
        {
            _weights.Add(XavierMatrix(layers[i], layers[i + 1], rnd));
            _biases.Add(new VectorN(layers[i + 1]));
        }
    }

    // ── Forward pass ────────────────────────────────────────────

    public VectorN Forward(VectorN input)
    {
        return Forward(input, out _);
    }

    public VectorN Forward(VectorN input, out List<VectorN> activations)
    {
        activations = new List<VectorN> { input };

        for (int i = 0; i < _weights.Count; i++)
        {
            var z = (_weights[i].Transpose() * activations[^1]) + _biases[i];

            VectorN a;
            if (i == _weights.Count - 1)
                a = Output == OutputMode.Softmax ? Softmax(z) : z;
            else
                a = Activate(z);

            activations.Add(a);
        }

        _lastActivations = activations;
        return activations[^1];
    }

    // ── Backward pass ───────────────────────────────────────────

    public void ComputeGradients(VectorN target, List<VectorN> activations,
        out List<Matrix> dW, out List<VectorN> dB)
    {
        var deltas = new List<VectorN>();
        var error = activations[^1] - target;
        deltas.Add(error);

        for (int i = _weights.Count - 2; i >= 0; i--)
        {
            var w = _weights[i + 1];
            var delta = (w * deltas[^1]).Hadamard(ActivationDerivative(activations[i + 1]));
            deltas.Add(delta);
        }
        deltas.Reverse();

        dW = new List<Matrix>();
        dB = new List<VectorN>();
        for (int i = 0; i < _weights.Count; i++)
        {
            dW.Add(activations[i].Outer(deltas[i]));
            dB.Add(deltas[i]);
        }
    }

    /// <summary>
    /// Backward pass from an external loss gradient (for RL policy gradient).
    /// Computes gradients and applies a single weight update.
    /// </summary>
    public void Backward(VectorN lossGradient, double learningRate)
    {
        if (_lastActivations == null)
            throw new InvalidOperationException("Call Forward before Backward.");

        var activations = _lastActivations;
        var deltas = new List<VectorN>();
        deltas.Add(lossGradient);

        for (int i = _weights.Count - 2; i >= 0; i--)
        {
            var w = _weights[i + 1];
            var delta = (w * deltas[^1]).Hadamard(ActivationDerivative(activations[i + 1]));
            deltas.Add(delta);
        }
        deltas.Reverse();

        for (int i = 0; i < _weights.Count; i++)
        {
            var dw = activations[i].Outer(deltas[i]);
            _weights[i] -= learningRate * dw;
            _biases[i] -= learningRate * deltas[i];
        }
    }

    // ── Weight update ───────────────────────────────────────────

    public void ApplyGradients(List<Matrix> dW, List<VectorN> dB,
        double learningRate, int batchSize, double l2 = 0.0)
    {
        for (int l = 0; l < _weights.Count; l++)
        {
            _weights[l] -= (learningRate / batchSize) * (dW[l] + l2 * _weights[l]);
            _biases[l] -= (learningRate / batchSize) * dB[l];
        }
    }

    /// <summary>
    /// Apply accumulated gradients using an <see cref="IOptimizer"/>.
    /// Each layer's weights and biases are flattened to a VectorN,
    /// stepped through the optimizer, then reshaped back.
    /// </summary>
    public void ApplyGradients(List<Matrix> dW, List<VectorN> dB,
        IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        for (int l = 0; l < _weights.Count; l++)
        {
            // Flatten weight matrix → VectorN for optimizer
            int rows = _weights[l].rowLength;
            int cols = _weights[l].columnLength;
            var wFlat = new double[rows * cols];
            var gFlat = new double[rows * cols];
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++)
                {
                    int idx = r * cols + c;
                    wFlat[idx] = _weights[l].values[r, c];
                    gFlat[idx] = dW[l].values[r, c] / batchSize;
                }

            var updatedW = weightOptimizer.Step(new VectorN(wFlat), new VectorN(gFlat));

            // Reshape back to matrix
            var newVals = new double[rows, cols];
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++)
                    newVals[r, c] = updatedW[r * cols + c];
            _weights[l] = new Matrix { values = newVals, rowLength = rows, columnLength = cols };

            // Bias update
            var bGrad = (1.0 / batchSize) * dB[l];
            _biases[l] = biasOptimizer.Step(_biases[l], bGrad);
        }
    }

    /// <summary>
    /// Backward pass from an external loss gradient using an <see cref="IOptimizer"/>.
    /// </summary>
    public void Backward(VectorN lossGradient, IOptimizer weightOptimizer, IOptimizer biasOptimizer)
    {
        if (_lastActivations == null)
            throw new InvalidOperationException("Call Forward before Backward.");

        var activations = _lastActivations;
        var deltas = new List<VectorN>();
        deltas.Add(lossGradient);

        for (int i = _weights.Count - 2; i >= 0; i--)
        {
            var w = _weights[i + 1];
            var delta = (w * deltas[^1]).Hadamard(ActivationDerivative(activations[i + 1]));
            deltas.Add(delta);
        }
        deltas.Reverse();

        var dW = new List<Matrix>();
        var dB = new List<VectorN>();
        for (int i = 0; i < _weights.Count; i++)
        {
            dW.Add(activations[i].Outer(deltas[i]));
            dB.Add(deltas[i]);
        }

        ApplyGradients(dW, dB, weightOptimizer, biasOptimizer, 1);
    }

    // ── Weight transfer (DQN target sync, DDPG Polyak) ─────────

    public void CopyWeightsFrom(NeuralNetwork source)
    {
        for (int i = 0; i < _weights.Count; i++)
        {
            _weights[i] = new Matrix(source._weights[i].values);
            _biases[i] = new VectorN(source._biases[i].Values);
        }
    }

    public void SoftUpdate(NeuralNetwork source, double tau)
    {
        for (int i = 0; i < _weights.Count; i++)
        {
            var sw = source._weights[i];
            var tw = _weights[i];
            var sb = source._biases[i];
            var tb = _biases[i];

            var newW = new double[tw.rowLength, tw.columnLength];
            for (int r = 0; r < tw.rowLength; r++)
                for (int c = 0; c < tw.columnLength; c++)
                    newW[r, c] = tau * sw.values[r, c] + (1 - tau) * tw.values[r, c];

            _weights[i] = new Matrix { values = newW, rowLength = tw.rowLength, columnLength = tw.columnLength };

            var newB = new double[tb.Length];
            for (int j = 0; j < tb.Length; j++)
                newB[j] = tau * sb.Values[j] + (1 - tau) * tb.Values[j];

            _biases[i] = new VectorN(newB);
        }
    }

    // ── Gradient helpers ────────────────────────────────────────

    public List<Matrix> InitWeightGrads() =>
        _weights.Select(m => new Matrix
        {
            values = new double[m.rowLength, m.columnLength],
            rowLength = m.rowLength,
            columnLength = m.columnLength
        }).ToList();

    public List<VectorN> InitBiasGrads() =>
        _biases.Select(v => new VectorN(v.Length)).ToList();

    public void AccumulateGradients(List<Matrix> accumW, List<VectorN> accumB,
        List<Matrix> dW, List<VectorN> dB)
    {
        for (int l = 0; l < _weights.Count; l++)
        {
            accumW[l] += dW[l];
            accumB[l] += dB[l];
        }
    }

    // ── Snapshot / restore ──────────────────────────────────────

    public (List<Matrix> weights, List<VectorN> biases) SnapshotWeights() =>
    (
        _weights.Select(m => new Matrix(m.values)).ToList(),
        _biases.Select(v => new VectorN(v.Values)).ToList()
    );

    public void RestoreWeights(List<Matrix> weights, List<VectorN> biases)
    {
        _weights = weights;
        _biases = biases;
    }

    // ── Clone ───────────────────────────────────────────────────

    public NeuralNetwork Clone()
    {
        var clone = new NeuralNetwork
        {
            InputSize = InputSize,
            OutputSize = OutputSize,
            HiddenLayers = (int[])HiddenLayers.Clone(),
            Activation = Activation,
            Output = Output,
        };
        clone._weights = _weights.Select(m => new Matrix(m.values)).ToList();
        clone._biases = _biases.Select(v => new VectorN(v.Values)).ToList();
        return clone;
    }

    // ── Activation functions ────────────────────────────────────

    private VectorN Activate(VectorN v)
    {
        return Activation switch
        {
            ActivationType.ReLU => ApplyElementwise(v, x => Math.Max(0, x)),
            ActivationType.Sigmoid => ApplyElementwise(v, x => 1.0 / (1.0 + Math.Exp(-x))),
            ActivationType.Tanh => ApplyElementwise(v, Math.Tanh),
            ActivationType.Linear => v,
            _ => throw new ArgumentOutOfRangeException()
        };
    }

    private VectorN ActivationDerivative(VectorN v)
    {
        return Activation switch
        {
            ActivationType.ReLU => ApplyElementwise(v, x => x > 0 ? 1.0 : 0.0),
            ActivationType.Sigmoid => ApplyElementwise(v, x => x * (1.0 - x)),
            ActivationType.Tanh => ApplyElementwise(v, x => 1.0 - x * x),
            ActivationType.Linear => ApplyElementwise(v, _ => 1.0),
            _ => throw new ArgumentOutOfRangeException()
        };
    }

    private static VectorN Softmax(VectorN z)
    {
        var values = z.Values;
        double max = values.Max();
        double[] exp = new double[values.Length];
        double sum = 0;
        for (int i = 0; i < values.Length; i++)
        {
            exp[i] = Math.Exp(values[i] - max);
            sum += exp[i];
        }
        for (int i = 0; i < values.Length; i++) exp[i] /= sum;
        return new VectorN(exp);
    }

    private static VectorN ApplyElementwise(VectorN v, Func<double, double> f)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = f(v.Values[i]);
        return new VectorN(result);
    }

    private static Matrix XavierMatrix(int rows, int cols, Random rnd)
    {
        var values = new double[rows, cols];
        double limit = Math.Sqrt(6.0 / (rows + cols));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i, j] = (rnd.NextDouble() * 2 - 1) * limit;
        return new Matrix { values = values, rowLength = rows, columnLength = cols };
    }
}
