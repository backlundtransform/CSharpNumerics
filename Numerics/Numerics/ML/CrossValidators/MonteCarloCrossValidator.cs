using CSharpNumerics.ML.CrossValidators.Interfaces;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics;
using CSharpNumerics.Statistics.MonteCarlo;
using CSharpNumerics.Statistics.Random;
using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.CrossValidators;

/// <summary>
/// Monte Carlo Cross-Validator: evaluates pipelines by running many random
/// train/test splits and collecting the resulting score distribution.
/// 
/// Unlike <see cref="ShuffleSplitCrossValidator"/> which returns a single
/// aggregate score, this validator leverages <see cref="MonteCarloSimulator"/>
/// to produce a full <see cref="MonteCarloResult"/> per pipeline — with
/// confidence intervals, percentiles, histograms, and standard error.
/// 
/// Implements <see cref="ICrossValidator"/> for drop-in compatibility with
/// existing code, plus an extended <see cref="RunDetailed"/> method that
/// returns the rich <see cref="MonteCarloValidationResult"/>.
/// 
/// <example>
/// <code>
/// var cv = new MonteCarloCrossValidator(pipelines, iterations: 200, testSize: 0.2, seed: 42);
/// 
/// // Standard interface — single best score
/// CrossValidationResult result = cv.Run(X, y);
/// 
/// // Extended — full score distributions + confidence intervals
/// MonteCarloValidationResult detailed = cv.RunDetailed(X, y);
/// var ci = detailed.BestConfidenceInterval; // (lower, upper)
/// </code>
/// </example>
/// </summary>
public class MonteCarloCrossValidator : ICrossValidator
{
    /// <summary>Pipelines to evaluate.</summary>
    public List<Pipeline> Pipelines { get; }

    /// <summary>Number of Monte Carlo iterations (random splits). Default 100.</summary>
    public int Iterations { get; }

    /// <summary>Fraction of data used for testing in each split. Default 0.2.</summary>
    public double TestSize { get; }

    /// <summary>Fraction of data used for training. Default: 1 - TestSize.</summary>
    public double TrainSize { get; }

    /// <summary>Optional seed for reproducible results.</summary>
    public int? Seed { get; }

    /// <summary>Confidence level used for the reported interval (default 0.95).</summary>
    public double ConfidenceLevel { get; }

    /// <summary>
    /// After <see cref="Run"/> or <see cref="RunDetailed"/> has been called,
    /// contains the detailed Monte Carlo result. Null before execution.
    /// </summary>
    public MonteCarloValidationResult DetailedResult { get; private set; }

    // ── Constructors ─────────────────────────────────────────────

    /// <summary>
    /// Creates a Monte Carlo cross-validator from an explicit list of pipelines.
    /// </summary>
    /// <param name="pipelines">The pipelines to evaluate.</param>
    /// <param name="iterations">Number of random splits (default 100).</param>
    /// <param name="testSize">Fraction of data held out for testing (default 0.2).</param>
    /// <param name="trainSize">Fraction of data used for training (default = 1 - testSize).</param>
    /// <param name="seed">Optional random seed for reproducibility.</param>
    /// <param name="confidenceLevel">CI level, e.g. 0.95 for 95 % (default 0.95).</param>
    public MonteCarloCrossValidator(
        List<Pipeline> pipelines,
        int iterations = 100,
        double testSize = 0.2,
        double? trainSize = null,
        int? seed = null,
        double confidenceLevel = 0.95)
    {
        Pipelines = pipelines ?? throw new ArgumentNullException(nameof(pipelines));
        Iterations = iterations > 0 ? iterations : throw new ArgumentException("iterations must be positive.");
        TestSize = testSize;
        TrainSize = trainSize ?? (1.0 - testSize);
        Seed = seed;
        ConfidenceLevel = confidenceLevel;
    }

    /// <summary>
    /// Creates a Monte Carlo cross-validator from a <see cref="PipelineGrid"/>.
    /// </summary>
    public MonteCarloCrossValidator(
        PipelineGrid pipelineGrid,
        int iterations = 100,
        double testSize = 0.2,
        double? trainSize = null,
        int? seed = null,
        double confidenceLevel = 0.95)
        : this(pipelineGrid.Expand().ToList(), iterations, testSize, trainSize, seed, confidenceLevel)
    {
    }

    // ── ICrossValidator ──────────────────────────────────────────

    /// <summary>
    /// Runs Monte Carlo cross-validation and returns a standard
    /// <see cref="CrossValidationResult"/>.
    /// The detailed result is stored in <see cref="DetailedResult"/>.
    /// </summary>
    public CrossValidationResult Run(Matrix X, VectorN y)
    {
        var detailed = RunDetailed(X, y);

        // Map to the standard CrossValidationResult for compatibility
        var result = new CrossValidationResult();

        foreach (var kvp in detailed.DetailedScores)
            result.Scores[kvp.Key] = kvp.Value.Mean;

        result.BestPipeline = detailed.BestPipeline;
        result.BestScore = detailed.BestMeanScore;

        // Refit best pipeline on full data for final predictions
        var bestClone = detailed.BestPipeline.Clone();
        bestClone.Fit(X, y);
        var predictions = bestClone.Predict(X);

        result.ActualValues = y;
        result.PredictedValues = predictions;

        var data = y.Values
            .Select((v, i) => (Pred: predictions.Values[i], Actual: v));

        result.CoefficientOfDetermination = data.CoefficientOfDetermination(p => (p.Pred, p.Actual));

        if (!(detailed.BestPipeline.Model is IRegressionModel))
        {
            result.ConfusionMatrix = y.BuildConfusionMatrix(predictions);
        }

        return result;
    }

    // ── Extended API ─────────────────────────────────────────────

    /// <summary>
    /// Runs Monte Carlo cross-validation and returns a rich
    /// <see cref="MonteCarloValidationResult"/> with full score
    /// distributions, confidence intervals and convergence curves.
    /// </summary>
    public MonteCarloValidationResult RunDetailed(Matrix X, VectorN y)
    {
        int n = y.Length;
        int testCount = Math.Max(1, (int)(TestSize * n));
        int trainCount = Math.Max(1, (int)(TrainSize * n));

        var detailedScores = new Dictionary<Pipeline, MonteCarloResult>();
        Pipeline bestPipeline = null;
        double bestMean = double.NegativeInfinity;

        // Use a master RNG to derive per-pipeline seeds so that every pipeline
        // sees the exact same sequence of splits (fair comparison).
        var masterRng = Seed.HasValue ? new RandomGenerator(Seed.Value) : new RandomGenerator();

        // Pre-generate split indices shared across all pipelines
        var splitIndices = GenerateSplitIndices(masterRng, n, testCount, trainCount, Iterations);

        foreach (var pipe in Pipelines)
        {
            var scores = new double[Iterations];

            for (int iter = 0; iter < Iterations; iter++)
            {
                var (trainIdx, testIdx) = splitIndices[iter];

                var Xtrain = X.SubMatrix(trainIdx);
                var Ytrain = new VectorN(trainIdx.Select(i => y.Values[i]).ToArray());
                var Xval = X.SubMatrix(testIdx);
                var Yval = new VectorN(testIdx.Select(i => y.Values[i]).ToArray());

                var cloned = pipe.Clone();
                cloned.Fit(Xtrain, Ytrain);
                var pred = cloned.Predict(Xval);

                scores[iter] = ScoreOptimized(pipe.Model is IRegressionModel, pred, Yval);
            }

            var mcResult = new MonteCarloResult(scores);
            detailedScores[pipe] = mcResult;

            if (mcResult.Mean > bestMean)
            {
                bestMean = mcResult.Mean;
                bestPipeline = pipe;
            }
        }

        // Build convergence curve for best pipeline
        var bestScores = detailedScores[bestPipeline].Samples;
        var convergence = new double[Iterations];
        double runningSum = 0;
        for (int i = 0; i < Iterations; i++)
        {
            runningSum += bestScores[i];
            convergence[i] = runningSum / (i + 1);
        }

        var bestMcResult = detailedScores[bestPipeline];
        var ci = bestMcResult.ConfidenceInterval(ConfidenceLevel);

        var detailed = new MonteCarloValidationResult
        {
            DetailedScores = detailedScores,
            BestPipeline = bestPipeline,
            BestMeanScore = bestMcResult.Mean,
            BestScoreStdDev = bestMcResult.StandardDeviation,
            BestConfidenceInterval = ci,
            ConvergenceCurve = convergence,
            Iterations = Iterations
        };

        DetailedResult = detailed;
        return detailed;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Pre-generates all random train/test splits so every pipeline
    /// is evaluated on the exact same data partitions.
    /// </summary>
    private static (int[] trainIdx, int[] testIdx)[] GenerateSplitIndices(
        RandomGenerator rng,
        int n,
        int testCount,
        int trainCount,
        int iterations)
    {
        var splits = new (int[] trainIdx, int[] testIdx)[iterations];

        for (int iter = 0; iter < iterations; iter++)
        {
            // Create a shuffled index array
            var indices = Enumerable.Range(0, n).ToArray();
            rng.Shuffle(indices);

            var testIdx = new int[testCount];
            Array.Copy(indices, 0, testIdx, 0, testCount);

            int actualTrain = Math.Min(trainCount, n - testCount);
            var trainIdx = new int[actualTrain];
            Array.Copy(indices, testCount, trainIdx, 0, actualTrain);

            splits[iter] = (trainIdx, testIdx);
        }

        return splits;
    }

    /// <summary>
    /// Computes the score for a single fold.
    /// Regression: negative MSE (higher is better).
    /// Classification: accuracy (higher is better).
    /// </summary>
    private static double ScoreOptimized(bool isRegression, VectorN pred, VectorN y)
    {
        if (isRegression)
        {
            double sum = 0;
            var predVals = pred.Values;
            var yVals = y.Values;
            for (int i = 0; i < pred.Length; i++)
            {
                double diff = predVals[i] - yVals[i];
                sum += diff * diff;
            }
            return -sum / pred.Length;
        }
        else
        {
            int correct = 0;
            var predVals = pred.Values;
            var yVals = y.Values;
            for (int i = 0; i < pred.Length; i++)
            {
                if (Math.Round(predVals[i]) == yVals[i])
                    correct++;
            }
            return (double)correct / pred.Length;
        }
    }
}
