using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.ML.Sequence.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// Temporal residual block (Bai et al., 2018): two dilated causal convolutions, each followed
/// by batch normalisation, a ReLU activation, and dropout, wrapped in a residual (skip)
/// connection. When the input and output channel counts differ, the skip path uses a 1×1
/// convolution to match dimensions. The skip connection lets gradients bypass the convolutions,
/// keeping deep dilation stacks trainable.
/// <code>
///   x → Conv(causal,dilated) → BN → ReLU → Dropout
///     → Conv(causal,dilated) → BN → Dropout → (+ skip(x)) → ReLU
/// </code>
/// </summary>
public class ResidualBlock : ILayer
{
    private readonly List<ILayer> _mainPath;
    private readonly Conv1DLayer _downsample;     // 1×1 conv on the skip path, or null for identity
    private readonly int _outputChannels;

    private VectorN[] _lastActivatedOutput = Array.Empty<VectorN>();

    public ResidualBlock(
        int inputChannels,
        int outputChannels,
        int kernelSize,
        int dilation,
        double dropoutRate = 0.0,
        ActivationType activation = ActivationType.ReLU,
        int seed = 123)
    {
        if (inputChannels <= 0) throw new ArgumentOutOfRangeException(nameof(inputChannels));
        if (outputChannels <= 0) throw new ArgumentOutOfRangeException(nameof(outputChannels));
        if (kernelSize <= 0) throw new ArgumentOutOfRangeException(nameof(kernelSize));
        if (dilation <= 0) throw new ArgumentOutOfRangeException(nameof(dilation));

        _outputChannels = outputChannels;

        _mainPath = new List<ILayer>
        {
            new Conv1DLayer(inputChannels, outputChannels, kernelSize, 1, ConvolutionPaddingMode.Causal, ActivationType.Linear, seed, dilation),
            new BatchNorm1DLayer(outputChannels),
            new ActivationLayer(activation),
            new DropoutLayer(dropoutRate, seed + 1),
            new Conv1DLayer(outputChannels, outputChannels, kernelSize, 1, ConvolutionPaddingMode.Causal, ActivationType.Linear, seed + 2, dilation),
            new BatchNorm1DLayer(outputChannels),
            new DropoutLayer(dropoutRate, seed + 3),
        };

        // 1×1 convolution to align channels on the skip path when they change.
        _downsample = inputChannels == outputChannels
            ? null
            : new Conv1DLayer(inputChannels, outputChannels, 1, 1, ConvolutionPaddingMode.Same, ActivationType.Linear, seed + 4);

        FinalActivation = activation;
    }

    public ActivationType FinalActivation { get; }

    /// <summary>Number of output channels emitted by the block.</summary>
    public int OutputChannels => _outputChannels;

    public int ParameterCount
    {
        get
        {
            int sum = 0;
            foreach (var layer in _mainPath) sum += layer.ParameterCount;
            if (_downsample != null) sum += _downsample.ParameterCount;
            return sum;
        }
    }

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must contain at least one timestep.", nameof(input));

        var main = input;
        foreach (var layer in _mainPath)
            main = layer.Forward(main, training);

        var skip = _downsample != null ? _downsample.Forward(input, training) : input;

        if (main.Length != skip.Length)
            throw new InvalidOperationException("Residual and main paths produced different sequence lengths.");

        var output = new VectorN[main.Length];
        for (int t = 0; t < main.Length; t++)
        {
            var sum = main[t] + skip[t];
            output[t] = Activations.Apply(sum, FinalActivation);
        }

        _lastActivatedOutput = CloneVectors(output);
        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_lastActivatedOutput.Length == 0)
            throw new InvalidOperationException("Call Forward before Backward.");
        if (gradOutput is null || gradOutput.Length != _lastActivatedOutput.Length)
            throw new ArgumentException("Gradient sequence must match the output sequence length.", nameof(gradOutput));

        // Gradient through the final activation, shared by both branches.
        var gradSum = new VectorN[gradOutput.Length];
        for (int t = 0; t < gradOutput.Length; t++)
            gradSum[t] = gradOutput[t].Hadamard(Activations.Derivative(_lastActivatedOutput[t], FinalActivation));

        // Main path (reverse order).
        var gradMain = CloneVectors(gradSum);
        for (int i = _mainPath.Count - 1; i >= 0; i--)
            gradMain = _mainPath[i].Backward(gradMain);

        // Skip path.
        VectorN[] gradSkip = _downsample != null
            ? _downsample.Backward(CloneVectors(gradSum))
            : gradSum;

        var gradInput = new VectorN[gradMain.Length];
        for (int t = 0; t < gradMain.Length; t++)
            gradInput[t] = gradMain[t] + gradSkip[t];

        return gradInput;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        foreach (var layer in _mainPath)
            layer.ApplyGradients(weightOptimizer, biasOptimizer, batchSize);
        _downsample?.ApplyGradients(weightOptimizer, biasOptimizer, batchSize);
    }

    private static VectorN[] CloneVectors(VectorN[] vectors)
    {
        var clone = new VectorN[vectors.Length];
        for (int i = 0; i < vectors.Length; i++)
            clone[i] = new VectorN(vectors[i].Values);
        return clone;
    }
}
