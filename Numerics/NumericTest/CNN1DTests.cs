using CSharpNumerics.ML;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Experiment;
using CSharpNumerics.ML.Sequence.Enums;
using CSharpNumerics.ML.Sequence.Layers;
using CSharpNumerics.ML.Sequence.Models.Classification;
using CSharpNumerics.ML.Sequence.Models.Regression;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class CNN1DTests
{
    [TestMethod]
    public void Conv1DLayer_Forward_ShouldMatchKnownValidConvolution()
    {
        var layer = new Conv1DLayer(
            new Matrix(new double[,]
            {
                { 1.0 },
                { 0.0 },
                { -1.0 }
            }),
            new VectorN(new[] { 0.0 }),
            inputChannels: 1,
            kernelSize: 3,
            stride: 1,
            padding: ConvolutionPaddingMode.Valid,
            activation: ActivationType.Linear);

        var output = layer.Forward(new[]
        {
            new VectorN(new[] { 1.0 }),
            new VectorN(new[] { 2.0 }),
            new VectorN(new[] { 4.0 }),
            new VectorN(new[] { 8.0 })
        });

        Assert.AreEqual(2, output.Length);
        Assert.AreEqual(-3.0, output[0][0], 1e-12);
        Assert.AreEqual(-6.0, output[1][0], 1e-12);
    }

    [TestMethod]
    public void MaxPool1DLayer_Backward_ShouldRouteGradientToWindowMaxima()
    {
        var layer = new MaxPool1DLayer(poolSize: 2, stride: 2);

        layer.Forward(new[]
        {
            new VectorN(new[] { 1.0 }),
            new VectorN(new[] { 3.0 }),
            new VectorN(new[] { 2.0 }),
            new VectorN(new[] { 4.0 })
        });

        var inputGradient = layer.Backward(new[]
        {
            new VectorN(new[] { 10.0 }),
            new VectorN(new[] { 20.0 })
        });

        CollectionAssert.AreEqual(new[] { 0.0, 10.0, 0.0, 20.0 }, inputGradient.Select(v => v[0]).ToArray());
    }

    [TestMethod]
    public void GlobalAvgPool1DLayer_ForwardBackward_ShouldAverageAndBroadcast()
    {
        var layer = new GlobalAvgPool1DLayer();

        var pooled = layer.Forward(new[]
        {
            new VectorN(new[] { 1.0, 2.0 }),
            new VectorN(new[] { 3.0, 4.0 })
        });

        Assert.AreEqual(1, pooled.Length);
        CollectionAssert.AreEqual(new[] { 2.0, 3.0 }, pooled[0].Values);

        var gradInput = layer.Backward(new[] { new VectorN(new[] { 6.0, 8.0 }) });
        CollectionAssert.AreEqual(new[] { 3.0, 4.0 }, gradInput[0].Values);
        CollectionAssert.AreEqual(new[] { 3.0, 4.0 }, gradInput[1].Values);
    }

    [TestMethod]
    public void FlattenLayer_ForwardBackward_ShouldPreserveTimestepOrder()
    {
        var layer = new FlattenLayer();

        var flattened = layer.Forward(new[]
        {
            new VectorN(new[] { 1.0, 2.0 }),
            new VectorN(new[] { 3.0, 4.0 })
        });

        Assert.AreEqual(1, flattened.Length);
        CollectionAssert.AreEqual(new[] { 1.0, 2.0, 3.0, 4.0 }, flattened[0].Values);

        var gradInput = layer.Backward(new[] { new VectorN(new[] { 10.0, 20.0, 30.0, 40.0 }) });
        CollectionAssert.AreEqual(new[] { 10.0, 20.0 }, gradInput[0].Values);
        CollectionAssert.AreEqual(new[] { 30.0, 40.0 }, gradInput[1].Values);
    }

    [TestMethod]
    public void CNN1DClassifier_SupervisedExperimentGrid_ShouldClassifyLocalSpikePatterns()
    {
        var (X, y) = CreateClassificationDataset(samples: 72, timeSteps: 12);

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<CNN1DClassifier>(g => g
                    .Add("TimeSteps", 12)
                    .Add("Features", 1)
                    .Add("Filters", 4, 6)
                    .Add("KernelSize", 5)
                    .Add("HiddenUnits", 8)
                    .Add("LearningRate", 0.03)
                    .Add("Epochs", 140)
                    .Add("BatchSize", 8)
                    .Add("ValidationSplit", 0.0)
                    .Add("UseGlobalAveragePooling", true)
                    .Add("Activation", ActivationType.ReLU)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 3))
            .Run();

        Assert.IsNotNull(result.BestPipeline);
        Assert.AreEqual(nameof(CNN1DClassifier), result.BestModelName);
        Assert.IsTrue(result.BestScore > 0.85, $"Expected CNN1DClassifier accuracy above 0.85, got {result.BestScore}.");
    }

    [TestMethod]
    public void CNN1DRegressor_KFold_ShouldFitSpikeAmplitude()
    {
        var (X, y) = CreateRegressionDataset(samples: 84, timeSteps: 12);

        var grid = new PipelineGrid()
            .AddModel<CNN1DRegressor>(g => g
                .Add("TimeSteps", 12)
                .Add("Features", 1)
                .Add("Filters", 6)
                .Add("KernelSize", 5)
                .Add("HiddenUnits", 8)
                .Add("LearningRate", 0.02)
                .Add("Epochs", 160)
                .Add("BatchSize", 8)
                .Add("ValidationSplit", 0.0)
                .Add("UseGlobalAveragePooling", true)
                .Add("Activation", ActivationType.ReLU));

        var result = new CSharpNumerics.ML.CrossValidators.KFoldCrossValidator(grid, folds: 3).Run(X, y);

        Assert.IsTrue(result.CoefficientOfDetermination > 0.75,
            $"Expected CNN1DRegressor R2 above 0.75, got {result.CoefficientOfDetermination}.");
    }

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