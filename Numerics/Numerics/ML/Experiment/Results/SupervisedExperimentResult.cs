using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Experiment.Results;

/// <summary>
/// Full supervised experiment output: all pipelines ranked by the primary cross-validator.
/// Analogous to <see cref="CSharpNumerics.ML.Clustering.Results.ClusteringExperimentResult"/>
/// but for supervised models with cross-validation.
/// </summary>
public class SupervisedExperimentResult
{
    /// <summary>All results ordered by primary cross-validator score (descending).</summary>
    public List<SupervisedResult> Rankings { get; set; } = new();

    /// <summary>The primary cross-validator used for ranking (first one added).</summary>
    public string PrimaryCrossValidator { get; set; }

    /// <summary>
    /// Raw <see cref="CrossValidationResult"/> per cross-validator (CV name → result).
    /// Each result contains per-pipeline scores, best pipeline, confusion matrix, etc.
    /// </summary>
    public Dictionary<string, CrossValidationResult> CVResults { get; set; } = new();

    // ── Convenience: best overall ────────────────────────────────

    /// <summary>Best result (highest primary CV score).</summary>
    public SupervisedResult Best => Rankings.FirstOrDefault();

    /// <summary>Best primary cross-validator score.</summary>
    public double BestScore => Best?.Scores.Values.FirstOrDefault() ?? 0;

    /// <summary>The best pipeline configuration.</summary>
    public Pipeline BestPipeline => Best?.Pipeline;

    /// <summary>Best model class name.</summary>
    public string BestModelName => Best?.ModelName;

    /// <summary>Actual values from the best pipeline's primary CV run.</summary>
    public VectorN BestActualValues
        => CVResults.TryGetValue(PrimaryCrossValidator, out var r) ? r.ActualValues : new VectorN(0);

    /// <summary>Predicted values from the best pipeline's primary CV run.</summary>
    public VectorN BestPredictedValues
        => CVResults.TryGetValue(PrimaryCrossValidator, out var r) ? r.PredictedValues : new VectorN(0);

    /// <summary>Confusion matrix from the best pipeline's primary CV run (classification only).</summary>
    public Matrix BestConfusionMatrix
        => CVResults.TryGetValue(PrimaryCrossValidator, out var r) ? r.ConfusionMatrix : new Matrix();

    /// <summary>R² from the best pipeline's primary CV run (regression only).</summary>
    public double BestR2
        => CVResults.TryGetValue(PrimaryCrossValidator, out var r) ? r.CoefficientOfDetermination : 0;

    // ── Best by specific cross-validator ─────────────────────────

    /// <summary>
    /// Find the best result according to a specific cross-validator name.
    /// </summary>
    public SupervisedResult BestBy(string cvName)
    {
        return Rankings
            .Where(r => r.Scores.ContainsKey(cvName))
            .OrderByDescending(r => r.Scores[cvName])
            .FirstOrDefault();
    }

    // ── Total experiment time ────────────────────────────────────

    /// <summary>Total wall-clock time for the entire experiment.</summary>
    public TimeSpan TotalDuration { get; set; }

    // ── Score distribution statistics ────────────────────────────

    /// <summary>
    /// Descriptive statistics (mean, median, std, IQR, skewness, kurtosis, …)
    /// of all pipeline scores for the given cross-validator.
    /// Answers: "how sensitive is the score to pipeline configuration?"
    /// </summary>
    /// <param name="cvName">
    /// Cross-validator name (e.g. "KFold"). Defaults to the primary CV.
    /// </param>
    public ScoreDistributionSummary ScoreSummary(string cvName = null)
    {
        cvName ??= PrimaryCrossValidator;
        var scores = Rankings
            .Where(r => r.Scores.ContainsKey(cvName))
            .Select(r => r.Scores[cvName])
            .ToList();
        return new ScoreDistributionSummary(scores);
    }

    /// <summary>
    /// Returns the percentile rank (0–100) of the given result among
    /// all rankings for the specified cross-validator.
    /// E.g. 95 means this pipeline scored better than 95 % of tested configurations.
    /// </summary>
    public double ScorePercentile(SupervisedResult result, string cvName = null)
    {
        cvName ??= PrimaryCrossValidator;
        if (!result.Scores.ContainsKey(cvName))
            return 0;

        double score = result.Scores[cvName];
        int belowCount = Rankings.Count(r =>
            r.Scores.ContainsKey(cvName) && r.Scores[cvName] < score);
        return 100.0 * belowCount / Rankings.Count;
    }

    /// <summary>
    /// Spearman rank correlation between the pipeline rankings produced
    /// by two different cross-validators.
    /// High ρ (close to 1) means the CVs agree on which pipelines are best.
    /// Low ρ means the ranking is sensitive to the CV strategy chosen.
    /// </summary>
    /// <returns>(ρ, p-value) from Spearman's rank correlation test.</returns>
    public (double Rho, double PValue) RankCorrelation(string cv1, string cv2)
    {
        var paired = Rankings
            .Where(r => r.Scores.ContainsKey(cv1) && r.Scores.ContainsKey(cv2))
            .Select(r => (x: r.Scores[cv1], y: r.Scores[cv2]))
            .ToList();

        if (paired.Count < 3)
            return (double.NaN, double.NaN);

        return paired.SpearmanCorrelation(p => (p.x, p.y));
    }

    /// <summary>
    /// Standard deviation of a single pipeline's scores across all
    /// cross-validators. Low = the pipeline performs consistently
    /// regardless of CV strategy. High = performance is CV-dependent.
    /// </summary>
    public double ScoreConsistency(SupervisedResult result)
    {
        var scores = result.Scores.Values.ToList();
        if (scores.Count < 2) return 0;
        double mean = scores.Average();
        return Math.Sqrt(scores.Sum(s => (s - mean) * (s - mean)) / (scores.Count - 1));
    }
}
