using CSharpNumerics.ML;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Experiment;
using CSharpNumerics.ML.Sequence;
using CSharpNumerics.ML.Sequence.Models.Classification;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Linq;

namespace NumericTest;

/// <summary>
/// End-to-end example: exoplanet-transit detection on a synthetic light curve.
///
/// A light curve records stellar flux over time. When an exoplanet transits its
/// host star, flux dips briefly (~1–4 hours out of an orbital period of days).
/// The task: classify fixed-length windows of the light curve as transit (1) or
/// no-transit (0).
///
/// Pipeline:
///   1. Synthesise a Kepler-like light curve using TimeSeries
///   2. Window the data using SequenceDataHelper.CreateWindows
///   3. Train a CNN1DClassifier via SupervisedExperiment grid search
///   4. Assert accuracy above threshold
/// </summary>
[TestClass]
public class ExoplanetTransitExampleTests
{
    [TestMethod]
    public void ExoplanetTransit_CNN1D_ShouldDetectTransitDips()
    {
        // ── 1. Synthesise a light curve ──────────────────────────────────
        // 2000 timesteps at ~30 min cadence (≈ 41 days of Kepler-like data).
        // Baseline flux ≈ 1.0 with Gaussian noise.
        // Inject periodic transit dips: box-shaped drop of ~0.01 (1%) lasting 6 points,
        // every 80 points (≈ 40 h orbital period).

        int totalTimesteps = 2000;
        int transitPeriod = 80;
        int transitDuration = 6;
        double transitDepth = 0.01;
        double noiseStd = 0.001;
        var rng = new Random(42);

        var times = new DateTime[totalTimesteps];
        var flux = new double[totalTimesteps];
        var labels = new double[totalTimesteps]; // 0 = no transit, 1 = transit

        DateTime start = new DateTime(2025, 1, 1);
        for (int i = 0; i < totalTimesteps; i++)
        {
            times[i] = start.AddMinutes(30 * i);

            // Baseline + Gaussian noise (Box-Muller)
            double u1 = 1.0 - rng.NextDouble();
            double u2 = rng.NextDouble();
            double noise = noiseStd * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            flux[i] = 1.0 + noise;

            // Transit dip
            int phase = i % transitPeriod;
            if (phase >= 0 && phase < transitDuration)
            {
                flux[i] -= transitDepth;
                labels[i] = 1.0;
            }
        }

        // ── 2. Create a TimeSeries and window it ────────────────────────
        var timeSeries = new TimeSeries(
            times,
            new[] { flux, labels },
            new[] { "Flux", "Label" });

        int windowSize = 20;
        int labelColumnIndex = 1; // "Label" column
        int stride = 5;

        var (X, y) = SequenceDataHelper.CreateWindows(timeSeries, windowSize, labelColumnIndex, stride);

        // Verify windowing produced the expected dimensions
        int expectedWindows = ((totalTimesteps - windowSize) / stride) + 1;
        Assert.AreEqual(expectedWindows, X.rowLength,
            $"Expected {expectedWindows} windows, got {X.rowLength}.");
        Assert.AreEqual(windowSize * 1, X.columnLength,
            "Each row should be windowSize × 1 feature (flux only).");

        // ── 3. Train a CNN1DClassifier ──────────────────────────────────
        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<CNN1DClassifier>(g => g
                    .Add("TimeSteps", windowSize)
                    .Add("Features", 1)
                    .Add("Filters", 8)
                    .Add("KernelSize", 5)
                    .Add("HiddenUnits", 8)
                    .Add("LearningRate", 0.02)
                    .Add("Epochs", 150)
                    .Add("BatchSize", 16)
                    .Add("ValidationSplit", 0.0)
                    .Add("UseGlobalAveragePooling", true)
                    .Add("Activation", ActivationType.ReLU)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 3))
            .Run();

        // ── 4. Assert ───────────────────────────────────────────────────
        Assert.IsNotNull(result.BestPipeline);
        Assert.AreEqual(nameof(CNN1DClassifier), result.BestModelName);
        Assert.IsTrue(result.BestScore > 0.80,
            $"Expected transit detection accuracy above 0.80, got {result.BestScore}.");
    }

    [TestMethod]
    public void SequenceDataHelper_CreateWindows_FromTimeSeries_CorrectDimensions()
    {
        // Simple verification that windowing math is correct
        int rows = 100;
        var times = Enumerable.Range(0, rows)
            .Select(i => new DateTime(2025, 1, 1).AddHours(i))
            .ToArray();
        var col0 = Enumerable.Range(0, rows).Select(i => (double)i).ToArray();
        var col1 = Enumerable.Range(0, rows).Select(i => i % 2 == 0 ? 0.0 : 1.0).ToArray();

        var ts = new TimeSeries(times, new[] { col0, col1 }, new[] { "Feature", "Label" });

        int windowSize = 10;
        int stride = 3;
        var (X, y) = SequenceDataHelper.CreateWindows(ts, windowSize, labelColumnIndex: 1, stride: stride);

        int expectedWindows = ((rows - windowSize) / stride) + 1; // (100-10)/3 + 1 = 31
        Assert.AreEqual(expectedWindows, X.rowLength);
        Assert.AreEqual(windowSize * 1, X.columnLength); // 1 feature column (Feature), label excluded

        // First window: timesteps 0..9, feature col0 = [0,1,2,...,9]
        for (int t = 0; t < windowSize; t++)
            Assert.AreEqual((double)t, X.values[0, t], 1e-12);

        // First window label = col1[9] = 1.0 (9 is odd)
        Assert.AreEqual(1.0, y[0], 1e-12);

        // Second window starts at stride=3, label = col1[3+9]=col1[12] = 0.0 (12 is even)
        Assert.AreEqual(0.0, y[1], 1e-12);
    }

    [TestMethod]
    public void SequenceDataHelper_CreateWindows_FromRawArrays_CorrectDimensions()
    {
        int rows = 50;
        var flux = Enumerable.Range(0, rows).Select(i => 1.0 + (0.01 * i)).ToArray();
        var labels = Enumerable.Range(0, rows).Select(i => i >= 25 ? 1.0 : 0.0).ToArray();

        int windowSize = 8;
        var (X, y) = SequenceDataHelper.CreateWindows(new[] { flux }, labels, windowSize, stride: 1);

        int expectedWindows = rows - windowSize + 1; // 43
        Assert.AreEqual(expectedWindows, X.rowLength);
        Assert.AreEqual(windowSize, X.columnLength);

        // Window at index 0: flux[0..7], label = labels[7] = 0.0
        Assert.AreEqual(0.0, y[0], 1e-12);

        // Window near the middle: starts at 18, label = labels[25] = 1.0
        Assert.AreEqual(1.0, y[18], 1e-12);
    }
}
