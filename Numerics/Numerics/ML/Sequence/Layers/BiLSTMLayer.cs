using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// Bidirectional LSTM layer. Runs two internal LSTMLayer instances — one forward
/// and one on the time-reversed input — then concatenates their hidden states per
/// timestep so that output_t = [h_fwd_t ∥ h_bwd_t] with dimension 2 × hiddenSize.
///
/// When returnSequences=false only the final concatenated vector is returned:
///   [h_fwd_T ∥ h_bwd_1]
/// </summary>
public class BiLSTMLayer : ILayer
{
    private readonly int _hiddenSize;
    private readonly bool _returnSequences;
    private readonly LSTMLayer _forward;
    private readonly LSTMLayer _backward;
    private int _timeSteps;

    public BiLSTMLayer(int inputSize, int hiddenSize, bool returnSequences = true, double clipNorm = 5.0, int seed = 123)
    {
        if (inputSize <= 0) throw new ArgumentOutOfRangeException(nameof(inputSize));
        if (hiddenSize <= 0) throw new ArgumentOutOfRangeException(nameof(hiddenSize));

        _hiddenSize = hiddenSize;
        _returnSequences = returnSequences;

        // Two LSTMLayer instances — both return sequences so we can concatenate per timestep
        _forward = new LSTMLayer(inputSize, hiddenSize, returnSequences: true, clipNorm: clipNorm, seed: seed);
        _backward = new LSTMLayer(inputSize, hiddenSize, returnSequences: true, clipNorm: clipNorm, seed: seed + 1);
    }

    public int ParameterCount => _forward.ParameterCount + _backward.ParameterCount;

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must not be empty.", nameof(input));

        _timeSteps = input.Length;
        int T = _timeSteps;

        // Forward LSTM processes input in order
        var fwdOutput = _forward.Forward(input, training);

        // Backward LSTM processes time-reversed input
        var reversedInput = ReverseSequence(input);
        var bwdOutputReversed = _backward.Forward(reversedInput, training);
        // Reverse the backward output to align with the original time axis
        var bwdOutput = ReverseSequence(bwdOutputReversed);

        // Concatenate per timestep: [h_fwd_t ∥ h_bwd_t]
        if (_returnSequences)
        {
            var output = new VectorN[T];
            for (int t = 0; t < T; t++)
                output[t] = fwdOutput[t].Concat(bwdOutput[t]);
            return output;
        }
        else
        {
            // Return single vector: [h_fwd_T ∥ h_bwd_1 (= first element of bwdOutput after reversal)]
            return new[] { fwdOutput[T - 1].Concat(bwdOutput[0]) };
        }
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        int T = _timeSteps;

        // Expand gradOutput to full sequence if returnSequences=false
        VectorN[] fullGrad;
        if (_returnSequences)
        {
            if (gradOutput.Length != T)
                throw new ArgumentException("Gradient length must match sequence length when returnSequences=true.");
            fullGrad = gradOutput;
        }
        else
        {
            if (gradOutput.Length != 1)
                throw new ArgumentException("Gradient must have length 1 when returnSequences=false.");

            fullGrad = new VectorN[T];
            for (int t = 0; t < T; t++)
                fullGrad[t] = new VectorN(2 * _hiddenSize);

            // Last forward timestep gets fwd gradient, first backward (index 0 after reversal) gets bwd gradient
            var fwdGrad = SliceVector(gradOutput[0], 0, _hiddenSize);
            var bwdGrad = SliceVector(gradOutput[0], _hiddenSize, _hiddenSize);

            fullGrad[T - 1] = fwdGrad.Concat(new VectorN(_hiddenSize));
            fullGrad[0] = new VectorN(_hiddenSize).Concat(bwdGrad);
        }

        // Split concatenated gradients into forward and backward parts
        var gradFwd = new VectorN[T];
        var gradBwd = new VectorN[T];
        for (int t = 0; t < T; t++)
        {
            gradFwd[t] = SliceVector(fullGrad[t], 0, _hiddenSize);
            gradBwd[t] = SliceVector(fullGrad[t], _hiddenSize, _hiddenSize);
        }

        // Forward LSTM backward pass
        var gradInputFwd = _forward.Backward(gradFwd);

        // Backward LSTM: its output was reversed, so its gradient must also be reversed
        var gradBwdReversed = ReverseSequence(gradBwd);
        var gradInputBwdReversed = _backward.Backward(gradBwdReversed);
        var gradInputBwd = ReverseSequence(gradInputBwdReversed);

        // Sum input gradients from both directions
        var gradInput = new VectorN[T];
        for (int t = 0; t < T; t++)
            gradInput[t] = gradInputFwd[t] + gradInputBwd[t];

        return gradInput;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        _forward.ApplyGradients(weightOptimizer, biasOptimizer, batchSize);
        _backward.ApplyGradients(weightOptimizer, biasOptimizer, batchSize);
    }

    private static VectorN[] ReverseSequence(VectorN[] sequence)
    {
        int length = sequence.Length;
        var reversed = new VectorN[length];
        for (int i = 0; i < length; i++)
            reversed[i] = sequence[length - 1 - i];
        return reversed;
    }

    private static VectorN SliceVector(VectorN v, int start, int length)
    {
        var result = new double[length];
        for (int i = 0; i < length; i++)
            result[i] = v[start + i];
        return new VectorN(result);
    }
}
