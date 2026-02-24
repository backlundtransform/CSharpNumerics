using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Clustering.Results;

/// <summary>
/// Full experiment output: all runs ranked by primary evaluator.
/// </summary>
public class ClusteringExperimentResult
{
    /// <summary>All results ordered by primary evaluator score (descending).</summary>
    public List<ClusteringResult> Rankings { get; set; } = new();

    /// <summary>The primary evaluator used for ranking.</summary>
    public string PrimaryEvaluator { get; set; }

    // ── Convenience: best overall ────────────────────────────────

    /// <summary>Best result (highest primary score).</summary>
    public ClusteringResult Best => Rankings.FirstOrDefault();

    /// <summary>Best cluster count found.</summary>
    public int BestClusterCount => Best?.ClusterCount ?? 0;

    /// <summary>Best primary evaluator score.</summary>
    public double BestScore => Best?.Scores.Values.FirstOrDefault() ?? 0;

    /// <summary>Cluster labels from the best run.</summary>
    public VectorN BestLabels => Best?.Labels ?? new VectorN(0);

    /// <summary>The best fitted model.</summary>
    public IClusteringModel BestModel => Best?.Model;

    // ── Best by specific evaluator ───────────────────────────────

    /// <summary>
    /// Find the best result according to a specific evaluator type.
    /// </summary>
    public ClusteringResult BestBy<TEvaluator>() where TEvaluator : IClusteringEvaluator, new()
    {
        string name = new TEvaluator().Name;
        return Rankings
            .Where(r => r.Scores.ContainsKey(name))
            .OrderByDescending(r => r.Scores[name])
            .FirstOrDefault();
    }

    /// <summary>
    /// Find the best result according to a named evaluator.
    /// </summary>
    public ClusteringResult BestBy(string evaluatorName)
    {
        return Rankings
            .Where(r => r.Scores.ContainsKey(evaluatorName))
            .OrderByDescending(r => r.Scores[evaluatorName])
            .FirstOrDefault();
    }

    // ── Elbow curve ──────────────────────────────────────────────

    /// <summary>
    /// K → raw inertia pairs for elbow method plotting.
    /// Only populated if <see cref="InertiaEvaluator"/> was included.
    /// </summary>
    public List<(int K, double Inertia)> ElbowCurve { get; set; } = new();

    // ── Total experiment time ────────────────────────────────────

    /// <summary>Total wall-clock time for the entire experiment.</summary>
    public TimeSpan TotalDuration { get; set; }

    // ── Monte Carlo uncertainty (optional) ───────────────────────

    /// <summary>
    /// Monte Carlo uncertainty analysis result. Only populated when
    /// <see cref="ClusteringExperimentBuilder.WithMonteCarloUncertainty"/>
    /// was called before <c>Run()</c>.
    /// Contains score distributions with confidence intervals,
    /// optimal-K distribution, and (for bootstrap) a consensus matrix.
    /// </summary>
    public MonteCarloClusteringResult MonteCarloResult { get; set; }
}
