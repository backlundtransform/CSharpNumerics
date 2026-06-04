using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// Channel-wise batch normalisation for sequence data. Each channel is normalised across
/// the timesteps of the sequence to zero mean and unit variance, then rescaled by a learnable
/// scale <c>γ</c> and shift <c>β</c>. During training the running mean/variance are accumulated
/// as an exponential moving average and used (instead of per-sequence statistics) at inference.
/// </summary>
public class BatchNorm1DLayer : ILayer
{
    private readonly int _channels;
    private readonly double _epsilon;
    private readonly double _momentum;

    private VectorN _gamma;
    private VectorN _beta;
    private VectorN _gammaGradients;
    private VectorN _betaGradients;

    private double[] _runningMean;
    private double[] _runningVariance;

    private VectorN[] _lastInput = Array.Empty<VectorN>();
    private double[] _lastMean;
    private double[] _lastStd;
    private double[][] _lastNormalized;

    private IOptimizer _gammaOptimizerState;
    private IOptimizer _betaOptimizerState;

    public BatchNorm1DLayer(int channels, double momentum = 0.9, double epsilon = 1e-5)
    {
        if (channels <= 0)
            throw new ArgumentOutOfRangeException(nameof(channels), "Channel count must be positive.");
        if (momentum < 0.0 || momentum > 1.0)
            throw new ArgumentOutOfRangeException(nameof(momentum), "Momentum must be in [0, 1].");

        _channels = channels;
        _momentum = momentum;
        _epsilon = epsilon;

        _gamma = OnesVector(channels);
        _beta = new VectorN(channels);
        _gammaGradients = new VectorN(channels);
        _betaGradients = new VectorN(channels);

        _runningMean = new double[channels];
        _runningVariance = new double[channels];
        for (int c = 0; c < channels; c++)
            _runningVariance[c] = 1.0;
    }

    public int ParameterCount => 2 * _channels;

    /// <summary>Running per-channel mean used at inference.</summary>
    public double[] RunningMean => (double[])_runningMean.Clone();

    /// <summary>Running per-channel variance used at inference.</summary>
    public double[] RunningVariance => (double[])_runningVariance.Clone();

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        ValidateInput(input);
        int n = input.Length;
        var output = new VectorN[n];

        if (!training)
        {
            // Inference: normalise with the running statistics.
            for (int t = 0; t < n; t++)
            {
                var values = new double[_channels];
                for (int c = 0; c < _channels; c++)
                {
                    double std = Math.Sqrt(_runningVariance[c] + _epsilon);
                    double normalized = (input[t][c] - _runningMean[c]) / std;
                    values[c] = (_gamma[c] * normalized) + _beta[c];
                }
                output[t] = new VectorN(values);
            }
            return output;
        }

        // Training: per-channel statistics over the time dimension.
        var mean = new double[_channels];
        var variance = new double[_channels];

        for (int t = 0; t < n; t++)
            for (int c = 0; c < _channels; c++)
                mean[c] += input[t][c];
        for (int c = 0; c < _channels; c++)
            mean[c] /= n;

        for (int t = 0; t < n; t++)
            for (int c = 0; c < _channels; c++)
            {
                double d = input[t][c] - mean[c];
                variance[c] += d * d;
            }
        for (int c = 0; c < _channels; c++)
            variance[c] /= n;

        var std2 = new double[_channels];
        var normalizedCache = new double[n][];
        for (int c = 0; c < _channels; c++)
            std2[c] = Math.Sqrt(variance[c] + _epsilon);

        for (int t = 0; t < n; t++)
        {
            var values = new double[_channels];
            var norm = new double[_channels];
            for (int c = 0; c < _channels; c++)
            {
                double normalized = (input[t][c] - mean[c]) / std2[c];
                norm[c] = normalized;
                values[c] = (_gamma[c] * normalized) + _beta[c];
            }
            normalizedCache[t] = norm;
            output[t] = new VectorN(values);
        }

        // Update running statistics (EMA).
        for (int c = 0; c < _channels; c++)
        {
            _runningMean[c] = (_momentum * _runningMean[c]) + ((1 - _momentum) * mean[c]);
            _runningVariance[c] = (_momentum * _runningVariance[c]) + ((1 - _momentum) * variance[c]);
        }

        _lastInput = CloneVectors(input);
        _lastMean = mean;
        _lastStd = std2;
        _lastNormalized = normalizedCache;

        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_lastInput.Length == 0)
            throw new InvalidOperationException("Call Forward (training) before Backward.");
        if (gradOutput is null || gradOutput.Length != _lastInput.Length)
            throw new ArgumentException("Gradient sequence must match the cached input length.", nameof(gradOutput));

        int n = _lastInput.Length;
        var gradInput = new double[n][];
        for (int t = 0; t < n; t++)
            gradInput[t] = new double[_channels];

        // Backprop is independent per channel over the N timesteps.
        for (int c = 0; c < _channels; c++)
        {
            double std = _lastStd[c];
            double invStd = 1.0 / std;

            double dGamma = 0.0;
            double dBeta = 0.0;
            double dVar = 0.0;
            double sumDxhat = 0.0;

            // First pass: parameter grads and accumulate dxhat terms.
            var dxhat = new double[n];
            for (int t = 0; t < n; t++)
            {
                double xhat = _lastNormalized[t][c];
                double dy = gradOutput[t][c];
                dGamma += dy * xhat;
                dBeta += dy;

                double dxh = dy * _gamma[c];
                dxhat[t] = dxh;
                sumDxhat += dxh;
                dVar += dxh * (_lastInput[t][c] - _lastMean[c]);
            }

            dVar *= -0.5 * invStd * invStd * invStd;            // dVar = Σ dxhat·(x-μ)·(-0.5)(σ²+ε)^-1.5
            double dMean = -invStd * sumDxhat;                  // Σ(x-μ) = 0 cancels the dVar term

            for (int t = 0; t < n; t++)
            {
                double centered = _lastInput[t][c] - _lastMean[c];
                gradInput[t][c] = (dxhat[t] * invStd) + (dVar * 2.0 * centered / n) + (dMean / n);
            }

            _gammaGradients[c] += dGamma;
            _betaGradients[c] += dBeta;
        }

        var result = new VectorN[n];
        for (int t = 0; t < n; t++)
            result[t] = new VectorN(gradInput[t]);
        return result;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        if (batchSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(batchSize), "Batch size must be positive.");

        _gammaOptimizerState ??= OptimizerFactory.Clone(weightOptimizer);
        _betaOptimizerState ??= OptimizerFactory.Clone(biasOptimizer);

        _gamma = _gammaOptimizerState.Step(_gamma, _gammaGradients / batchSize);
        _beta = _betaOptimizerState.Step(_beta, _betaGradients / batchSize);

        _gammaGradients = new VectorN(_channels);
        _betaGradients = new VectorN(_channels);
    }

    private void ValidateInput(VectorN[] input)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));
        for (int t = 0; t < input.Length; t++)
            if (input[t].Length != _channels)
                throw new ArgumentException("Each timestep must match the configured channel count.", nameof(input));
    }

    private static VectorN OnesVector(int n)
    {
        var values = new double[n];
        for (int i = 0; i < n; i++) values[i] = 1.0;
        return new VectorN(values);
    }

    private static VectorN[] CloneVectors(VectorN[] vectors)
    {
        var clone = new VectorN[vectors.Length];
        for (int i = 0; i < vectors.Length; i++)
            clone[i] = new VectorN(vectors[i].Values);
        return clone;
    }
}
