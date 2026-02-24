using CSharpNumerics.Statistics.MonteCarlo;
using System.Collections.Generic;

namespace CSharpNumerics.ML.CrossValidators.Result;

/// <summary>
/// Extended result from <see cref="MonteCarloCrossValidator"/> that
/// exposes the full <see cref="MonteCarloResult"/> per pipeline.
/// This gives access to score distributions, confidence intervals,
/// histograms and convergence information — not just a single mean score.
/// </summary>
public class MonteCarloValidationResult
{
    /// <summary>
    /// Per-pipeline Monte Carlo result containing the full distribution
    /// of scores across all iterations.
    /// </summary>
    public Dictionary<Pipeline, MonteCarloResult> DetailedScores { get; set; } = new();

    /// <summary>Best pipeline based on mean score.</summary>
    public Pipeline BestPipeline { get; set; }

    /// <summary>Mean score of the best pipeline.</summary>
    public double BestMeanScore { get; set; }

    /// <summary>Standard deviation of the best pipeline's score distribution.</summary>
    public double BestScoreStdDev { get; set; }

    /// <summary>
    /// Confidence interval for the best pipeline's score.
    /// Default 95 % confidence level.
    /// </summary>
    public (double lower, double upper) BestConfidenceInterval { get; set; }

    /// <summary>
    /// Convergence curve (running mean) for the best pipeline.
    /// Index i = mean of first i+1 iterations.
    /// Useful for verifying that enough MC iterations were run.
    /// </summary>
    public double[] ConvergenceCurve { get; set; }

    /// <summary>Number of Monte Carlo iterations that were executed.</summary>
    public int Iterations { get; set; }
}
