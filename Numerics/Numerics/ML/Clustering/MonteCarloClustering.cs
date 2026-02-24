using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Clustering.Results;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.MonteCarlo;
using CSharpNumerics.Statistics.Random;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CSharpNumerics.ML.Clustering;

/// <summary>
/// Monte Carlo uncertainty estimation for clustering.
/// 
/// Provides two analysis modes:
/// <list type="number">
///   <item>
///     <see cref="RunBootstrap"/> — Runs a single algorithm N times on bootstrap
///     samples, producing a <b>consensus matrix</b>, per-point <b>stability scores</b>,
///     and a <b>score distribution</b> with confidence intervals.
///   </item>
///   <item>
///     <see cref="RunExperiment"/> — Runs a full K-range experiment N times on
///     bootstrap samples, producing an <b>optimal-K distribution</b> that shows
///     how often each K "wins" across random resamples.
///   </item>
/// </list>
/// 
/// Both modes leverage the library's <see cref="MonteCarloSimulator"/> and
/// <see cref="MonteCarloResult"/> for statistical summaries.
/// 
/// <example>
/// <code>
/// // Bootstrap: consensus matrix + score distribution
/// var mc = new MonteCarloClustering { Iterations = 200, Seed = 42 };
/// var result = mc.RunBootstrap(
///     data,
///     new KMeans { K = 3 },
///     new SilhouetteEvaluator(),
///     new StandardScaler());
/// 
/// var ci = result.ScoreConfidenceInterval();       // (0.68, 0.74)
/// Matrix consensus = result.ConsensusMatrix;        // N×N
/// double[] stability = result.PointStability;       // per-point [0,1]
/// 
/// // Experiment: optimal-K distribution
/// var kResult = mc.RunExperiment(
///     data,
///     new KMeans(),
///     new SilhouetteEvaluator(),
///     minK: 2, maxK: 8);
/// 
/// foreach (var (k, count) in kResult.OptimalKDistribution.OrderByDescending(x => x.Value))
///     Console.WriteLine($"K={k}: chosen {count}/{mc.Iterations} times");
/// </code>
/// </example>
/// </summary>
public class MonteCarloClustering
{
    /// <summary>Number of Monte Carlo iterations. Default 100.</summary>
    public int Iterations { get; set; } = 100;

    /// <summary>Optional random seed for reproducibility.</summary>
    public int? Seed { get; set; }

    /// <summary>Confidence level for reported intervals. Default 0.95.</summary>
    public double ConfidenceLevel { get; set; } = 0.95;

    // ══════════════════════════════════════════════════════════════
    //  RunBootstrap — consensus matrix + score distribution
    // ══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs the clustering algorithm <paramref name="iterations"/> times
    /// on bootstrap-resampled data. Produces:
    /// <list type="bullet">
    ///   <item>A <b>score distribution</b> (e.g. Silhouette) as <see cref="MonteCarloResult"/>.</item>
    ///   <item>A <b>consensus matrix</b> (N×N) where cell (i,j) = fraction of times i and j co-clustered.</item>
    ///   <item>Per-point <b>stability</b> — how consistently each point lands in the same cluster.</item>
    ///   <item>A <b>convergence curve</b> for the score.</item>
    /// </list>
    /// </summary>
    /// <param name="data">The data matrix (n rows × d columns).</param>
    /// <param name="algorithm">Clustering algorithm to use (will be cloned per iteration).</param>
    /// <param name="evaluator">Evaluator to score each iteration.</param>
    /// <param name="scaler">Optional scaler applied before clustering.</param>
    public MonteCarloClusteringResult RunBootstrap(
        Matrix data,
        IClusteringModel algorithm,
        IClusteringEvaluator evaluator,
        IScaler scaler = null)
    {
        if (evaluator == null) throw new ArgumentNullException(nameof(evaluator));
        if (algorithm == null) throw new ArgumentNullException(nameof(algorithm));

        int n = data.rowLength;
        var rng = Seed.HasValue ? new RandomGenerator(Seed.Value) : new RandomGenerator();

        var scores = new double[Iterations];

        // Consensus: count how many times each pair co-occurs in identical cluster
        var coAssignment = new double[n, n];
        var coOccurrence = new double[n, n]; // how many times both i and j appear in bootstrap

        for (int iter = 0; iter < Iterations; iter++)
        {
            // 1. Bootstrap sample (with replacement)
            int[] indices = SampleWithReplacement(rng, n, n);

            Matrix xBoot = data.SubMatrix(indices);

            // 2. Scale if needed
            if (scaler != null)
            {
                var scalerClone = scaler.Clone();
                xBoot = scalerClone.FitTransform(xBoot);
            }

            // 3. Fit & predict
            var modelClone = algorithm.Clone();
            VectorN labels = modelClone.FitPredict(xBoot);

            // 4. Score
            scores[iter] = evaluator.Score(xBoot, labels);

            // 5. Update consensus tracking
            // Map bootstrap position → original index and label
            for (int a = 0; a < indices.Length; a++)
            {
                for (int b = a + 1; b < indices.Length; b++)
                {
                    int origA = indices[a];
                    int origB = indices[b];

                    coOccurrence[origA, origB] += 1;
                    coOccurrence[origB, origA] += 1;

                    if ((int)labels[a] == (int)labels[b] && (int)labels[a] >= 0)
                    {
                        coAssignment[origA, origB] += 1;
                        coAssignment[origB, origA] += 1;
                    }
                }
            }
        }

        // Build consensus matrix (normalize)
        var consensusValues = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            consensusValues[i, i] = 1.0; // point always co-clusters with itself
            for (int j = i + 1; j < n; j++)
            {
                double val = coOccurrence[i, j] > 0
                    ? coAssignment[i, j] / coOccurrence[i, j]
                    : 0.0;
                consensusValues[i, j] = val;
                consensusValues[j, i] = val;
            }
        }
        var consensusMatrix = new Matrix(consensusValues);

        // Compute final labels from consensus (run algorithm once on full data)
        Matrix dataForLabels = data;
        if (scaler != null)
        {
            var scalerFinal = scaler.Clone();
            dataForLabels = scalerFinal.FitTransform(data);
        }
        var finalModel = algorithm.Clone();
        VectorN finalLabels = finalModel.FitPredict(dataForLabels);

        // Compute per-point stability
        double[] pointStability = ComputePointStability(consensusValues, finalLabels, n);

        // Build convergence curve
        var convergence = new double[Iterations];
        double runningSum = 0;
        for (int i = 0; i < Iterations; i++)
        {
            runningSum += scores[i];
            convergence[i] = runningSum / (i + 1);
        }

        var scoreResult = new MonteCarloResult(scores);

        return new MonteCarloClusteringResult
        {
            ScoreDistribution = scoreResult,
            ConsensusMatrix = consensusMatrix,
            PointStability = pointStability,
            ConvergenceCurve = convergence,
            Iterations = Iterations,
            EvaluatorName = evaluator.Name
        };
    }

    // ══════════════════════════════════════════════════════════════
    //  RunExperiment — optimal-K distribution
    // ══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs a full K-range clustering experiment N times on bootstrap samples.
    /// For each iteration, determines which K produces the best evaluator score.
    /// Returns the distribution of "winning" K values across all iterations,
    /// plus the score distribution for the overall best K.
    /// </summary>
    /// <param name="data">The data matrix (n rows × d columns).</param>
    /// <param name="algorithm">Clustering algorithm template (must support a settable K property).</param>
    /// <param name="evaluator">Evaluator used to rank K values.</param>
    /// <param name="minK">Minimum number of clusters to try.</param>
    /// <param name="maxK">Maximum number of clusters to try.</param>
    /// <param name="scaler">Optional scaler applied before clustering.</param>
    public MonteCarloClusteringResult RunExperiment(
        Matrix data,
        IClusteringModel algorithm,
        IClusteringEvaluator evaluator,
        int minK = 2,
        int maxK = 8,
        IScaler scaler = null)
    {
        if (evaluator == null) throw new ArgumentNullException(nameof(evaluator));
        if (algorithm == null) throw new ArgumentNullException(nameof(algorithm));
        if (minK < 2) minK = 2;
        if (maxK < minK) maxK = minK;

        int n = data.rowLength;
        var rng = Seed.HasValue ? new RandomGenerator(Seed.Value) : new RandomGenerator();

        var kDistribution = new Dictionary<int, int>();
        for (int k = minK; k <= maxK; k++)
            kDistribution[k] = 0;

        // Track scores for each K across iterations
        var kScores = new Dictionary<int, List<double>>();
        for (int k = minK; k <= maxK; k++)
            kScores[k] = new List<double>(Iterations);

        for (int iter = 0; iter < Iterations; iter++)
        {
            // Bootstrap sample
            int[] indices = SampleWithReplacement(rng, n, n);
            Matrix xBoot = data.SubMatrix(indices);

            if (scaler != null)
            {
                var scalerClone = scaler.Clone();
                xBoot = scalerClone.FitTransform(xBoot);
            }

            double bestScore = double.NegativeInfinity;
            int bestK = minK;

            for (int k = minK; k <= maxK; k++)
            {
                var modelClone = algorithm.Clone();
                SetK(modelClone, k);

                VectorN labels = modelClone.FitPredict(xBoot);
                double score = evaluator.Score(xBoot, labels);

                kScores[k].Add(score);

                if (score > bestScore)
                {
                    bestScore = score;
                    bestK = k;
                }
            }

            kDistribution[bestK]++;
        }

        // Find the K that won most often
        int overallBestK = kDistribution.OrderByDescending(x => x.Value).First().Key;
        var bestKScores = kScores[overallBestK].ToArray();
        var scoreResult = new MonteCarloResult(bestKScores);

        // Convergence curve for the most-chosen K
        var convergence = new double[bestKScores.Length];
        double runningSum = 0;
        for (int i = 0; i < bestKScores.Length; i++)
        {
            runningSum += bestKScores[i];
            convergence[i] = runningSum / (i + 1);
        }

        return new MonteCarloClusteringResult
        {
            ScoreDistribution = scoreResult,
            OptimalKDistribution = kDistribution,
            ConvergenceCurve = convergence,
            Iterations = Iterations,
            EvaluatorName = evaluator.Name
        };
    }

    // ══════════════════════════════════════════════════════════════
    //  Private helpers
    // ══════════════════════════════════════════════════════════════

    /// <summary>
    /// Bootstrap sampling: draw n indices from [0, n) with replacement.
    /// </summary>
    private static int[] SampleWithReplacement(RandomGenerator rng, int populationSize, int sampleSize)
    {
        var result = new int[sampleSize];
        for (int i = 0; i < sampleSize; i++)
            result[i] = rng.NextInt(populationSize);
        return result;
    }

    /// <summary>
    /// Computes per-point stability from the consensus matrix.
    /// For each point, averages the consensus value with all other points
    /// sharing the same final cluster label.
    /// </summary>
    private static double[] ComputePointStability(double[,] consensus, VectorN labels, int n)
    {
        var stability = new double[n];

        for (int i = 0; i < n; i++)
        {
            int ci = (int)labels[i];
            if (ci < 0) { stability[i] = 0; continue; } // noise

            double sum = 0;
            int count = 0;
            for (int j = 0; j < n; j++)
            {
                if (j == i) continue;
                if ((int)labels[j] == ci)
                {
                    sum += consensus[i, j];
                    count++;
                }
            }

            stability[i] = count > 0 ? sum / count : 1.0;
        }

        return stability;
    }

    /// <summary>
    /// Sets the K property on a clustering model via reflection.
    /// </summary>
    private static void SetK(IClusteringModel model, int k)
    {
        var prop = model.GetType().GetProperty("K");
        if (prop != null && prop.CanWrite)
            prop.SetValue(model, k);
    }
}
