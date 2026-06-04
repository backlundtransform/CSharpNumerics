using System;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Sequence.Enums;
using CSharpNumerics.ML.Sequence.Layers;
using CSharpNumerics.ML.Sequence.Models.Regression;
using CSharpNumerics.Numerics.Objects;

namespace NumericsTests;

[TestClass]
public class TCNTests
{
    private static VectorN[] RandomSequence(int length, int channels, Random rng)
    {
        var seq = new VectorN[length];
        for (int t = 0; t < length; t++)
        {
            var values = new double[channels];
            for (int c = 0; c < channels; c++)
                values[c] = rng.NextDouble() * 2.0 - 1.0;
            seq[t] = new VectorN(values);
        }
        return seq;
    }

    #region Causal convolution

    [TestMethod]
    public void CausalConv1D_OutputDoesNotDependOnFutureTimesteps()
    {
        // Linear activation so differences are not masked by ReLU.
        var conv = new Conv1DLayer(
            inputChannels: 1, filters: 1, kernelSize: 3, stride: 1,
            padding: ConvolutionPaddingMode.Causal, activation: ActivationType.Linear, seed: 7);

        var rng = new Random(1);
        var signal = RandomSequence(10, 1, rng);

        var baseline = conv.Forward(signal, training: false);

        // Perturb a future timestep (index 7) and re-run.
        var perturbed = (VectorN[])signal.Clone();
        perturbed[7] = new VectorN(new[] { signal[7][0] + 5.0 });
        var modified = conv.Forward(perturbed, training: false);

        // Outputs at indices < 7 must be unchanged (strictly causal).
        for (int t = 0; t <= 6; t++)
            Assert.AreEqual(baseline[t][0], modified[t][0], 1e-12,
                $"Causal output at t={t} should not depend on a future input at t=7.");

        // The perturbed timestep itself must change.
        Assert.IsTrue(Math.Abs(baseline[7][0] - modified[7][0]) > 1e-9,
            "Output at the perturbed timestep should change.");
    }

    [TestMethod]
    public void DilatedConv1D_HasExpectedReceptiveField()
    {
        var conv = new Conv1DLayer(
            inputChannels: 1, filters: 1, kernelSize: 3, stride: 1,
            padding: ConvolutionPaddingMode.Causal, activation: ActivationType.Linear, seed: 1, dilation: 4);

        // (kernelSize - 1) * dilation + 1 = 2*4 + 1 = 9
        Assert.AreEqual(9, conv.ReceptiveField);
    }

    #endregion

    #region TCN block

    [TestMethod]
    public void TCNBlock_EightLevels_ReceptiveFieldAtLeast512()
    {
        var tcn = new TCNBlock(inputChannels: 1, channels: 4, kernelSize: 3, levels: 8);

        Assert.AreEqual(128, tcn.Dilations[^1], "Eighth level should have dilation 128.");
        Assert.IsTrue(tcn.ReceptiveField >= 512,
            $"Receptive field ({tcn.ReceptiveField}) should cover at least 512 timesteps.");
    }

    [TestMethod]
    public void ResidualBlock_GradientFlowsThroughSkipConnection()
    {
        // Identity skip (in == out channels), no dropout.
        var block = new ResidualBlock(inputChannels: 4, outputChannels: 4, kernelSize: 2, dilation: 1, dropoutRate: 0.0, seed: 3);

        var rng = new Random(5);
        var input = RandomSequence(8, 4, rng);

        var output = block.Forward(input, training: true);
        Assert.AreEqual(input.Length, output.Length);
        Assert.AreEqual(4, output[0].Length);

        // Upstream gradient of all ones.
        var gradOut = new VectorN[output.Length];
        for (int t = 0; t < output.Length; t++)
            gradOut[t] = new VectorN(new[] { 1.0, 1.0, 1.0, 1.0 });

        var gradIn = block.Backward(gradOut);

        double maxAbs = 0.0;
        foreach (var g in gradIn)
            foreach (var v in g.Values)
            {
                Assert.IsFalse(double.IsNaN(v) || double.IsInfinity(v), "Gradient must be finite.");
                maxAbs = Math.Max(maxAbs, Math.Abs(v));
            }

        // The skip connection passes gradient ~1:1, so the input gradient cannot vanish.
        Assert.IsTrue(maxAbs > 0.5,
            $"Input gradient magnitude ({maxAbs:F4}) indicates the skip connection is not propagating gradient.");
    }

    #endregion

    #region Dropout & BatchNorm

    [TestMethod]
    public void Dropout_ZerosApproximatelyRateFraction_DuringTraining()
    {
        const double rate = 0.5;
        var drop = new DropoutLayer(rate, seed: 42);

        int length = 200, channels = 10;
        var input = new VectorN[length];
        for (int t = 0; t < length; t++)
        {
            var values = new double[channels];
            for (int c = 0; c < channels; c++) values[c] = 1.0;
            input[t] = new VectorN(values);
        }

        var trained = drop.Forward(input, training: true);

        int zeros = 0, total = length * channels;
        foreach (var v in trained)
            foreach (var x in v.Values)
                if (x == 0.0) zeros++;

        double fraction = (double)zeros / total;
        Assert.IsTrue(Math.Abs(fraction - rate) < 0.05,
            $"Zeroed fraction ({fraction:P1}) should be close to the dropout rate ({rate:P0}).");

        // Inference is a pass-through (no zeros, unchanged values).
        var inference = drop.Forward(input, training: false);
        foreach (var v in inference)
            foreach (var x in v.Values)
                Assert.AreEqual(1.0, x, 1e-12, "Dropout should be a no-op at inference.");
    }

    [TestMethod]
    public void BatchNorm_NormalizesPerChannel_AndRunningStatsConverge()
    {
        var bn = new BatchNorm1DLayer(channels: 2, momentum: 0.9);

        int length = 50;
        var rng = new Random(11);
        var input = new VectorN[length];
        for (int t = 0; t < length; t++)
        {
            double ch0 = 5.0 + 2.0 * (rng.NextDouble() * 2.0 - 1.0);   // mean ≈ 5
            double ch1 = -3.0 + 0.5 * (rng.NextDouble() * 2.0 - 1.0);  // mean ≈ -3
            input[t] = new VectorN(new[] { ch0, ch1 });
        }

        var output = bn.Forward(input, training: true);

        // Output per channel should be ~zero mean, ~unit variance (γ=1, β=0).
        for (int c = 0; c < 2; c++)
        {
            double mean = 0.0;
            for (int t = 0; t < length; t++) mean += output[t][c];
            mean /= length;

            double var = 0.0;
            for (int t = 0; t < length; t++) var += (output[t][c] - mean) * (output[t][c] - mean);
            var /= length;

            Assert.AreEqual(0.0, mean, 1e-9, $"Channel {c} output mean should be ~0.");
            Assert.AreEqual(1.0, var, 1e-3, $"Channel {c} output variance should be ~1.");
        }

        // Running statistics converge to the data statistics over repeated forwards.
        for (int i = 0; i < 300; i++) bn.Forward(input, training: true);

        Assert.AreEqual(5.0, bn.RunningMean[0], 0.3, "Running mean of channel 0 should approach the data mean.");
        Assert.AreEqual(-3.0, bn.RunningMean[1], 0.3, "Running mean of channel 1 should approach the data mean.");
    }

    #endregion

    #region TCN regressor

    [TestMethod]
    public void TCNRegressor_TrainsAndConvergesOnSineWave()
    {
        // Next-step prediction of a clean sine. The model should beat a mean predictor.
        const int window = 16;
        int total = 180;
        var signal = new double[total];
        for (int i = 0; i < total; i++) signal[i] = Math.Sin(2.0 * Math.PI * i / 20.0);

        int samples = total - window;
        var x = new double[samples, window];
        var yValues = new double[samples];
        for (int i = 0; i < samples; i++)
        {
            for (int j = 0; j < window; j++) x[i, j] = signal[i + j];
            yValues[i] = signal[i + window];
        }

        var X = new Matrix(x);
        var y = new VectorN(yValues);

        var model = new TCNRegressor
        {
            TimeSteps = window,
            Features = 1,
            Channels = 8,
            KernelSize = 2,
            Levels = 3,
            HiddenUnits = 8,
            DropoutRate = 0.0,
            LearningRate = 0.01,
            Epochs = 150,
            BatchSize = 16,
            ValidationSplit = 0.0,
            Seed = 123
        };

        model.Fit(X, y);
        var predictions = model.Predict(X);

        double meanY = 0.0;
        for (int i = 0; i < samples; i++) meanY += yValues[i];
        meanY /= samples;

        double sseModel = 0.0, sseBaseline = 0.0;
        for (int i = 0; i < samples; i++)
        {
            sseModel += Math.Pow(predictions[i] - yValues[i], 2);
            sseBaseline += Math.Pow(meanY - yValues[i], 2);
        }
        double rmseModel = Math.Sqrt(sseModel / samples);
        double rmseBaseline = Math.Sqrt(sseBaseline / samples);

        // Global average pooling discards temporal position, which caps absolute accuracy on
        // next-step phase prediction — but the TCN should still clearly beat a mean predictor.
        Assert.IsTrue(rmseModel < 0.9 * rmseBaseline,
            $"TCN RMSE ({rmseModel:F4}) should clearly beat the mean-predictor baseline ({rmseBaseline:F4}).");
    }

    [TestMethod]
    public void TCNRegressor_FullForwardBackwardWithDilationStack_DoesNotCrash()
    {
        const int window = 32;
        int samples = 40;
        var rng = new Random(2);
        var x = new double[samples, window];
        var yValues = new double[samples];
        for (int i = 0; i < samples; i++)
        {
            for (int j = 0; j < window; j++) x[i, j] = rng.NextDouble();
            yValues[i] = rng.NextDouble();
        }

        var model = new TCNRegressor
        {
            TimeSteps = window,
            Features = 1,
            Channels = 8,
            KernelSize = 3,
            Levels = 4,            // dilations 1, 2, 4, 8
            HiddenUnits = 8,
            DropoutRate = 0.1,
            Epochs = 2,
            BatchSize = 8,
            ValidationSplit = 0.0,
            Seed = 9
        };

        model.Fit(new Matrix(x), new VectorN(yValues));
        var predictions = model.Predict(new Matrix(x));

        Assert.AreEqual(samples, predictions.Length);
        foreach (var p in predictions.Values)
            Assert.IsFalse(double.IsNaN(p) || double.IsInfinity(p), "Predictions must be finite.");
    }

    #endregion
}
