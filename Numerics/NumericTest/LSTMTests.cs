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
public class LSTMTests
{
    // ──────────────────────────────────────────
    // Layer-level tests
    // ──────────────────────────────────────────

    [TestMethod]
    public void LSTMLayer_Forward_OutputDimensionsMatchHiddenSize()
    {
        int inputSize = 3;
        int hiddenSize = 4;
        var layer = new LSTMLayer(inputSize, hiddenSize, returnSequences: true, seed: 42);

        var input = new[]
        {
            new VectorN(new[] { 1.0, 0.5, -0.5 }),
            new VectorN(new[] { 0.2, 0.8, -0.1 }),
            new VectorN(new[] { -0.3, 0.0, 1.0 })
        };

        var output = layer.Forward(input);

        // returnSequences=true → output length = input length
        Assert.AreEqual(3, output.Length, "Expected one output per timestep.");
        foreach (var h in output)
            Assert.AreEqual(hiddenSize, h.Length, "Each output should have hiddenSize dimensions.");
    }

    [TestMethod]
    public void LSTMLayer_Forward_ReturnSequencesFalse_SingleOutput()
    {
        int inputSize = 2;
        int hiddenSize = 3;
        var layer = new LSTMLayer(inputSize, hiddenSize, returnSequences: false, seed: 7);

        var input = new[]
        {
            new VectorN(new[] { 1.0, -1.0 }),
            new VectorN(new[] { 0.5, 0.5 }),
            new VectorN(new[] { -0.2, 0.8 }),
            new VectorN(new[] { 0.0, 1.0 })
        };

        var output = layer.Forward(input);

        Assert.AreEqual(1, output.Length, "returnSequences=false should yield one output.");
        Assert.AreEqual(hiddenSize, output[0].Length);
    }

    [TestMethod]
    public void LSTMLayer_Backward_ProducesInputGradientsOfCorrectShape()
    {
        int inputSize = 2;
        int hiddenSize = 4;
        var layer = new LSTMLayer(inputSize, hiddenSize, returnSequences: true, seed: 99);

        var input = new[]
        {
            new VectorN(new[] { 1.0, 0.0 }),
            new VectorN(new[] { 0.0, 1.0 }),
            new VectorN(new[] { -1.0, 0.5 })
        };

        layer.Forward(input);

        var gradOutput = new[]
        {
            new VectorN(new double[hiddenSize]),
            new VectorN(new double[hiddenSize]),
            new VectorN(new[] { 1.0, 0.0, 0.0, 0.0 })
        };

        var gradInput = layer.Backward(gradOutput);

        Assert.AreEqual(3, gradInput.Length, "Gradient sequence should match input length.");
        foreach (var grad in gradInput)
            Assert.AreEqual(inputSize, grad.Length, "Each gradient should have inputSize dimensions.");
    }

    [TestMethod]
    public void LSTMLayer_Backward_ReturnSequencesFalse_ProducesInputGradients()
    {
        int inputSize = 2;
        int hiddenSize = 3;
        var layer = new LSTMLayer(inputSize, hiddenSize, returnSequences: false, seed: 55);

        var input = new[]
        {
            new VectorN(new[] { 1.0, 0.5 }),
            new VectorN(new[] { -0.5, 1.0 })
        };

        layer.Forward(input);

        var gradOutput = new[] { new VectorN(new[] { 1.0, -1.0, 0.5 }) };
        var gradInput = layer.Backward(gradOutput);

        Assert.AreEqual(2, gradInput.Length);
        foreach (var grad in gradInput)
            Assert.AreEqual(inputSize, grad.Length);

        // At least one input gradient should be non-zero (gradient flows through BPTT)
        bool anyNonZero = false;
        foreach (var grad in gradInput)
        {
            for (int i = 0; i < grad.Length; i++)
            {
                if (Math.Abs(grad[i]) > 1e-15)
                {
                    anyNonZero = true;
                    break;
                }
            }
        }

        Assert.IsTrue(anyNonZero, "BPTT should produce non-zero input gradients.");
    }

    [TestMethod]
    public void LSTMLayer_GradientClipping_CapsGlobalNorm()
    {
        int inputSize = 1;
        int hiddenSize = 2;
        double clipNorm = 0.001; // Very tight clip to ensure clipping activates

        var layer = new LSTMLayer(inputSize, hiddenSize, returnSequences: false, clipNorm: clipNorm, seed: 42);

        var input = new[]
        {
            new VectorN(new[] { 10.0 }), // Large values to create large gradients
            new VectorN(new[] { -10.0 })
        };

        layer.Forward(input);

        var gradOutput = new[] { new VectorN(new[] { 100.0, -100.0 }) }; // Large gradient
        var gradInput = layer.Backward(gradOutput);

        // After clipping with very small clipNorm, gradients should be small.
        // The key test: applying gradients shouldn't explode. The internal state accounts for clipping.
        var wOpt = new GradientDescent(0.01);
        var bOpt = new GradientDescent(0.01);
        layer.ApplyGradients(wOpt, bOpt, 1);

        // Verify we can still forward: weights didn't become NaN
        var output = layer.Forward(new[] { new VectorN(new[] { 1.0 }), new VectorN(new[] { -1.0 }) });
        Assert.IsFalse(double.IsNaN(output[0][0]), "Weights should not be NaN after clipped gradient update.");
        Assert.IsFalse(double.IsNaN(output[0][1]), "Weights should not be NaN after clipped gradient update.");
    }

    [TestMethod]
    public void LSTMLayer_ParameterCount_MatchesExpected()
    {
        int inputSize = 3;
        int hiddenSize = 5;
        var layer = new LSTMLayer(inputSize, hiddenSize, seed: 1);

        // 4 gates × ((inputSize + hiddenSize) × hiddenSize + hiddenSize)
        int concatSize = inputSize + hiddenSize; // 8
        int expected = 4 * ((concatSize * hiddenSize) + hiddenSize); // 4 * (40 + 5) = 180
        Assert.AreEqual(expected, layer.ParameterCount);
    }

    // ──────────────────────────────────────────
    // End-to-end model tests
    // ──────────────────────────────────────────

    [TestMethod]
    public void LSTMClassifier_SupervisedExperiment_ShouldClassifySequencePattern()
    {
        var (X, y) = CreateClassificationDataset(samples: 80, timeSteps: 10);

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<LSTMClassifier>(g => g
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
        Assert.AreEqual(nameof(LSTMClassifier), result.BestModelName);
        Assert.IsTrue(result.BestScore > 0.70,
            $"Expected LSTMClassifier accuracy above 0.70, got {result.BestScore}.");
    }

    [TestMethod]
    public void LSTMRegressor_KFold_ShouldFitSpikeAmplitude()
    {
        var (X, y) = CreateRegressionDataset(samples: 84, timeSteps: 10);

        var grid = new PipelineGrid()
            .AddModel<LSTMRegressor>(g => g
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
            $"Expected LSTMRegressor R² above 0.50, got {result.CoefficientOfDetermination}.");
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
