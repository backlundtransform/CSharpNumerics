using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// Temporal Convolutional Network block: a stack of <see cref="ResidualBlock"/>s whose dilation
/// grows exponentially (1, 2, 4, 8, …). Each level doubles the dilation, so the receptive field
/// grows exponentially with depth while the sequence length is preserved (causal, stride 1).
/// This is the primary feature extractor for long-range temporal dependencies.
/// </summary>
public class TCNBlock : ILayer
{
    private readonly List<ResidualBlock> _blocks;
    private readonly int _kernelSize;
    private readonly int[] _dilations;

    /// <summary>
    /// Builds a TCN with <paramref name="levels"/> residual blocks of equal channel width.
    /// Dilation at level i is 2^i.
    /// </summary>
    /// <param name="inputChannels">Number of input channels (features).</param>
    /// <param name="channels">Channel width of every residual block.</param>
    /// <param name="kernelSize">Convolution kernel size (≥ 2 to grow the receptive field).</param>
    /// <param name="levels">Number of residual blocks (dilation doubles each level).</param>
    /// <param name="dropoutRate">Dropout rate inside each residual block.</param>
    /// <param name="activation">Activation used inside the residual blocks.</param>
    /// <param name="seed">Base random seed.</param>
    public TCNBlock(
        int inputChannels,
        int channels,
        int kernelSize,
        int levels,
        double dropoutRate = 0.0,
        ActivationType activation = ActivationType.ReLU,
        int seed = 123)
    {
        if (inputChannels <= 0) throw new ArgumentOutOfRangeException(nameof(inputChannels));
        if (channels <= 0) throw new ArgumentOutOfRangeException(nameof(channels));
        if (kernelSize <= 0) throw new ArgumentOutOfRangeException(nameof(kernelSize));
        if (levels <= 0) throw new ArgumentOutOfRangeException(nameof(levels));

        _kernelSize = kernelSize;
        _blocks = new List<ResidualBlock>(levels);
        _dilations = new int[levels];

        int inChannels = inputChannels;
        for (int level = 0; level < levels; level++)
        {
            int dilation = 1 << level;   // 2^level
            _dilations[level] = dilation;
            _blocks.Add(new ResidualBlock(inChannels, channels, kernelSize, dilation, dropoutRate, activation, seed + (level * 10)));
            inChannels = channels;
        }
    }

    /// <summary>Number of output channels.</summary>
    public int OutputChannels => _blocks[_blocks.Count - 1].OutputChannels;

    /// <summary>Number of residual blocks (levels).</summary>
    public int Levels => _blocks.Count;

    /// <summary>The dilation factor used at each level.</summary>
    public int[] Dilations => (int[])_dilations.Clone();

    /// <summary>
    /// Total receptive field in timesteps. Each residual block contains two dilated convolutions,
    /// so it contributes <c>2 × (kernelSize - 1) × dilation</c>; summed over all levels and
    /// offset by 1 for the centre tap:
    /// <c>RF = 1 + 2(kernelSize - 1) Σ dilations</c>.
    /// </summary>
    public int ReceptiveField
    {
        get
        {
            int dilationSum = 0;
            foreach (int d in _dilations) dilationSum += d;
            return 1 + (2 * (_kernelSize - 1) * dilationSum);
        }
    }

    public int ParameterCount
    {
        get
        {
            int sum = 0;
            foreach (var block in _blocks) sum += block.ParameterCount;
            return sum;
        }
    }

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        var activations = input;
        foreach (var block in _blocks)
            activations = block.Forward(activations, training);
        return activations;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        var gradients = gradOutput;
        for (int i = _blocks.Count - 1; i >= 0; i--)
            gradients = _blocks[i].Backward(gradients);
        return gradients;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        foreach (var block in _blocks)
            block.ApplyGradients(weightOptimizer, biasOptimizer, batchSize);
    }
}
