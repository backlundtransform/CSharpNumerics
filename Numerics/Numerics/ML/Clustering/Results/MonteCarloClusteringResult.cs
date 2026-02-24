using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.MonteCarlo;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering.Results;

/// <summary>
/// Result of a Monte Carlo clustering analysis.
/// Contains score distributions with confidence intervals, a consensus matrix
/// for cluster stability analysis, and an optimal-K distribution.
/// 
/// <example>
/// <code>
/// var mc = new MonteCarloClustering { Iterations = 200, Seed = 42 };
/// var result = mc.RunBootstrap(data, new KMeans { K = 3 }, new SilhouetteEvaluator());
/// 
/// // Score uncertainty
/// var ci = result.ScoreConfidenceInterval();   // e.g. (0.68, 0.74)
/// double se = result.ScoreDistribution.StandardError;
/// 
/// // Point stability
/// double[] stability = result.PointStability;  // per-point cluster membership stability
/// 
/// // Consensus matrix
/// Matrix consensus = result.ConsensusMatrix;   // N×N co-assignment frequency matrix
/// </code>
/// </example>
/// </summary>
public class MonteCarloClusteringResult
{
    /// <summary>
    /// Distribution of evaluator scores across all MC iterations.
    /// Provides Mean, StdDev, Percentile, ConfidenceInterval, Histogram, StandardError.
    /// </summary>
    public MonteCarloResult ScoreDistribution { get; set; }

    /// <summary>
    /// Consensus matrix (N × N): element (i, j) = fraction of iterations
    /// where data point i and data point j were assigned to the same cluster.
    /// Values in [0, 1]. High diagonal blocks indicate stable clusters.
    /// Only populated by <see cref="MonteCarloClustering.RunBootstrap"/>.
    /// </summary>
    public Matrix ConsensusMatrix { get; set; }

    /// <summary>
    /// Per-point stability score in [0, 1].
    /// For each point, this is the average consensus value with all other points
    /// that share its final (majority) cluster assignment.
    /// A value close to 1.0 means the point consistently belongs to the same cluster.
    /// Only populated by <see cref="MonteCarloClustering.RunBootstrap"/>.
    /// </summary>
    public double[] PointStability { get; set; }

    /// <summary>
    /// Distribution of optimal K across MC iterations.
    /// Key = K value, Value = number of times that K was selected as best.
    /// Only populated by <see cref="MonteCarloClustering.RunExperiment"/>.
    /// </summary>
    public Dictionary<int, int> OptimalKDistribution { get; set; }

    /// <summary>
    /// Convergence curve: running mean of the primary score.
    /// Index i = mean of first i+1 iterations.
    /// Useful for verifying that enough iterations were run.
    /// </summary>
    public double[] ConvergenceCurve { get; set; }

    /// <summary>Number of Monte Carlo iterations executed.</summary>
    public int Iterations { get; set; }

    /// <summary>Name of the evaluator used for scoring.</summary>
    public string EvaluatorName { get; set; }

    /// <summary>
    /// Returns a confidence interval for the primary score.
    /// </summary>
    /// <param name="level">Confidence level, e.g. 0.95 for 95 %. Default 0.95.</param>
    public (double lower, double upper) ScoreConfidenceInterval(double level = 0.95)
        => ScoreDistribution.ConfidenceInterval(level);
}
