using CSharpNumerics.ML;
using CSharpNumerics.ML.CrossValidators;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Experiment;
using CSharpNumerics.ML.Sequence.Layers;
using CSharpNumerics.ML.Sequence.Models.Classification;
using CSharpNumerics.ML.Sequence.Models.Regression;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.SingleObjective;

namespace NumericTest;

[TestClass]
public class BiLSTMTests
{
    // ──────────────────────────────────────────
    // Layer-level tests
    // ──────────────────────────────────────────

    [TestMethod]
    public void BiLSTMLayer_Forward_ReturnSequences_OutputDimIs2xHidden()
    {
        int inputSize = 3;
        int hiddenSize = 4;
        var layer = new BiLSTMLayer(inputSize, hiddenSize, returnSequences: true, seed: 42);

        var input = new[]
        {
            new VectorN(new[] { 1.0, 0.5, -0.5 }),
            new VectorN(new[] { 0.2, 0.8, -0.1 }),
            new VectorN(new[] { -0.3, 0.0, 1.0 })
        };

        var output = layer.Forward(input);

        Assert.AreEqual(3, output.Length, "Expected one output per timestep.");
        foreach (var h in output)
            Assert.AreEqual(2 * hiddenSize, h.Length, "Each output should have 2 × hiddenSize dimensions.");
    }

    [TestMethod]
    public void BiLSTMLayer_Forward_ReturnSequencesFalse_SingleConcatenatedOutput()
    {
        int inputSize = 2;
        int hiddenSize = 3;
        var layer = new BiLSTMLayer(inputSize, hiddenSize, returnSequences: false, seed: 7);

        var input = new[]
        {
            new VectorN(new[] { 1.0, -1.0 }),
            new VectorN(new[] { 0.5, 0.5 }),
            new VectorN(new[] { -0.2, 0.8 })
        };

        var output = layer.Forward(input);

        Assert.AreEqual(1, output.Length, "returnSequences=false should yield one output.");
        Assert.AreEqual(2 * hiddenSize, output[0].Length, "Output should be [h_fwd_T || h_bwd_1].");
    }

    [TestMethod]
    public void BiLSTMLayer_ForwardOutputDiffersFromUnidirectional()
    {
        int inputSize = 2;
        int hiddenSize = 4;
        int seed = 99;

        var uniLayer = new LSTMLayer(inputSize, hiddenSize, returnSequences: true, seed: seed);
        var biLayer = new BiLSTMLayer(inputSize, hiddenSize, returnSequences: true, seed: seed);

        var input = new[]
        {
            new VectorN(new[] { 1.0, 0.0 }),
            new VectorN(new[] { 0.0, 1.0 }),
            new VectorN(new[] { -1.0, 0.5 })
        };

        var uniOut = uniLayer.Forward(input);
        var biOut = biLayer.Forward(input);

        // BiLSTM output is wider (2×hidden vs hidden)
        Assert.AreEqual(hiddenSize, uniOut[0].Length);
        Assert.AreEqual(2 * hiddenSize, biOut[0].Length);

        // First half of BiLSTM output should match unidirectional (same seed for forward LSTM)
        for (int t = 0; t < input.Length; t++)
        {
            for (int i = 0; i < hiddenSize; i++)
                Assert.AreEqual(uniOut[t][i], biOut[t][i], 1e-12,
                    $"Forward half of BiLSTM should match unidirectional at t={t}, i={i}.");
        }
    }

    [TestMethod]
    public void BiLSTMLayer_Backward_ProducesCorrectShapeGradients()
    {
        int inputSize = 2;
        int hiddenSize = 3;
        var layer = new BiLSTMLayer(inputSize, hiddenSize, returnSequences: true, seed: 55);

        var input = new[]
        {
            new VectorN(new[] { 1.0, 0.5 }),
            new VectorN(new[] { -0.5, 1.0 }),
            new VectorN(new[] { 0.3, -0.3 })
        };

        layer.Forward(input);

        var gradOutput = new[]
        {
            new VectorN(new double[2 * hiddenSize]),
            new VectorN(new double[2 * hiddenSize]),
            new VectorN(new[] { 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 })
        };

        var gradInput = layer.Backward(gradOutput);

        Assert.AreEqual(3, gradInput.Length);
        foreach (var grad in gradInput)
            Assert.AreEqual(inputSize, grad.Length, "Gradient width should match input size.");
    }

    [TestMethod]
    public void BiLSTMLayer_Backward_ReturnSequencesFalse_ProducesGradients()
    {
        int inputSize = 2;
        int hiddenSize = 3;
        var layer = new BiLSTMLayer(inputSize, hiddenSize, returnSequences: false, seed: 33);

        var input = new[]
        {
            new VectorN(new[] { 1.0, -0.5 }),
            new VectorN(new[] { 0.5, 0.5 })
        };

        layer.Forward(input);

        var gradOutput = new[] { new VectorN(new[] { 1.0, -1.0, 0.5, 0.2, -0.3, 0.8 }) };
        var gradInput = layer.Backward(gradOutput);

        Assert.AreEqual(2, gradInput.Length);
        foreach (var grad in gradInput)
            Assert.AreEqual(inputSize, grad.Length);

        // Gradients should be non-zero (both directions contribute)
        bool anyNonZero = false;
        foreach (var grad in gradInput)
            for (int i = 0; i < grad.Length; i++)
                if (Math.Abs(grad[i]) > 1e-15)
                { anyNonZero = true; break; }

        Assert.IsTrue(anyNonZero, "BPTT from both directions should produce non-zero gradients.");
    }

    [TestMethod]
    public void BiLSTMLayer_ParameterCount_IsDoubleUnidirectional()
    {
        int inputSize = 3;
        int hiddenSize = 5;
        var uniLayer = new LSTMLayer(inputSize, hiddenSize, seed: 1);
        var biLayer = new BiLSTMLayer(inputSize, hiddenSize, seed: 1);

        Assert.AreEqual(2 * uniLayer.ParameterCount, biLayer.ParameterCount,
            "BiLSTM should have exactly twice the parameters of a unidirectional LSTM.");
    }

    // ──────────────────────────────────────────
    // End-to-end model tests
    // ──────────────────────────────────────────

    [TestMethod]
    public void BiLSTMClassifier_SupervisedExperiment_ShouldClassifySequencePattern()
    {
        var (X, y) = CreateClassificationDataset(samples: 80, timeSteps: 10);

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<BiLSTMClassifier>(g => g
                    .Add("TimeSteps", 10)
                    .Add("Features", 1)
                    .Add("HiddenSize", 8)
                    .Add("HiddenUnits", 8)
                    .Add("LearningRate", 0.01)
                    .Add("Epochs", 200)
                    .Add("BatchSize", 8)
                    .Add("ClipNorm", 5.0)
                    .Add("ValidationSplit", 0.0)
                    .Add("Activation", ActivationType.ReLU)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 3))
            .Run();

        Assert.IsNotNull(result.BestPipeline);
        Assert.AreEqual(nameof(BiLSTMClassifier), result.BestModelName);
        Assert.IsTrue(result.BestScore > 0.70,
            $"Expected BiLSTMClassifier accuracy above 0.70, got {result.BestScore}.");
    }

    [TestMethod]
    public void BiLSTMRegressor_KFold_ShouldFitSpikeAmplitude()
    {
        var (X, y) = CreateRegressionDataset(samples: 84, timeSteps: 10);

        var grid = new PipelineGrid()
            .AddModel<BiLSTMRegressor>(g => g
                .Add("TimeSteps", 10)
                .Add("Features", 1)
                .Add("HiddenSize", 8)
                .Add("HiddenUnits", 8)
                .Add("LearningRate", 0.01)
                .Add("Epochs", 200)
                .Add("BatchSize", 8)
                .Add("ClipNorm", 5.0)
                .Add("ValidationSplit", 0.0)
                .Add("Activation", ActivationType.ReLU));

        var result = new KFoldCrossValidator(grid, folds: 3).Run(X, y);

        Assert.IsTrue(result.CoefficientOfDetermination > 0.50,
            $"Expected BiLSTMRegressor R² above 0.50, got {result.CoefficientOfDetermination}.");
    }

    // ──────────────────────────────────────────
    // Test data generators
    // ──────────────────────────────────────────

    private static (Matrix X, VectorN y) CreateClassificationDataset(int samples, int timeSteps)
    {
        var X = new Matrix(samples, timeSteps);
        var y = new VectorN(samples);
        var random = new Random(7);

        for (int sample = 0; sample < samples; sample++)
        {
            bool hasSpike = sample >= samples / 2;
            y[sample] = hasSpike ? 1.0 : 0.0;

            for (int timestep = 0; timestep < timeSteps; timestep++)
                X.values[sample, timestep] = (random.NextDouble() - 0.5) * 0.02;

            if (hasSpike)
            {
                int start = random.Next(1, timeSteps - 5);
                X.values[sample, start] += 0.5;
                X.values[sample, start + 1] += 1.5;
                X.values[sample, start + 2] += 3.0;
                X.values[sample, start + 3] += 1.5;
                X.values[sample, start + 4] += 0.5;
            }
        }

        return (X, y);
    }

    private static (Matrix X, VectorN y) CreateRegressionDataset(int samples, int timeSteps)
    {
        var X = new Matrix(samples, timeSteps);
        var y = new VectorN(samples);
        var random = new Random(11);

        for (int sample = 0; sample < samples; sample++)
        {
            double amplitude = 0.25 + (1.75 * random.NextDouble());
            y[sample] = amplitude;

            for (int timestep = 0; timestep < timeSteps; timestep++)
                X.values[sample, timestep] = (random.NextDouble() - 0.5) * 0.01;

            int start = random.Next(1, timeSteps - 5);
            X.values[sample, start] += 0.25 * amplitude;
            X.values[sample, start + 1] += amplitude;
            X.values[sample, start + 2] += 2.0 * amplitude;
            X.values[sample, start + 3] += amplitude;
            X.values[sample, start + 4] += 0.25 * amplitude;
        }

        return (X, y);
    }
}
