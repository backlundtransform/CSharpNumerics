using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.ML.Clustering.Results;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.MonteCarlo;

namespace NumericTest;

[TestClass]
public class MonteCarloClusteringTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Test data helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Three well-separated 2D clusters:
    ///   Cluster 0: around (0, 0)
    ///   Cluster 1: around (10, 10)
    ///   Cluster 2: around (20, 0)
    /// </summary>
    private static Matrix ThreeClusterData()
    {
        var data = new double[,]
        {
            // Cluster 0
            { 0.0, 0.0 }, { 0.5, 0.3 }, { -0.3, 0.5 }, { 0.2, -0.4 }, { -0.1, 0.1 },
            { 0.4, -0.2 }, { -0.5, -0.3 }, { 0.1, 0.6 }, { -0.2, -0.1 }, { 0.3, 0.2 },
            // Cluster 1
            { 10.0, 10.0 }, { 10.5, 10.3 }, { 9.7, 10.5 }, { 10.2, 9.6 }, { 9.9, 10.1 },
            { 10.4, 9.8 }, { 9.5, 9.7 }, { 10.1, 10.6 }, { 9.8, 9.9 }, { 10.3, 10.2 },
            // Cluster 2
            { 20.0, 0.0 }, { 20.5, 0.3 }, { 19.7, 0.5 }, { 20.2, -0.4 }, { 19.9, 0.1 },
            { 20.4, -0.2 }, { 19.5, -0.3 }, { 20.1, 0.6 }, { 19.8, -0.1 }, { 20.3, 0.2 },
        };
        return new Matrix(data);
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunBootstrap — Score Distribution
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_ShouldReturnScoreDistribution()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        Assert.IsNotNull(result.ScoreDistribution);
        Assert.AreEqual(50, result.ScoreDistribution.Count);
        Assert.AreEqual(50, result.Iterations);
        Assert.AreEqual("Silhouette", result.EvaluatorName);
    }

    [TestMethod]
    public void RunBootstrap_ScoreDistribution_ShouldHavePositiveStdDev()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        // There IS variance across bootstrap samples
        Assert.IsTrue(result.ScoreDistribution.StandardDeviation > 0,
            "Score distribution should have non-zero std dev");
    }

    [TestMethod]
    public void RunBootstrap_ShouldProvideConfidenceInterval()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 100, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        var ci = result.ScoreConfidenceInterval(0.95);
        double mean = result.ScoreDistribution.Mean;

        Assert.IsTrue(ci.lower < mean, $"CI lower ({ci.lower}) should be < mean ({mean})");
        Assert.IsTrue(ci.upper > mean, $"CI upper ({ci.upper}) should be > mean ({mean})");
        Assert.IsTrue(ci.upper - ci.lower > 0, "CI width must be positive");
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunBootstrap — Consensus Matrix
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_ShouldReturnConsensusMatrix()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        Assert.IsNotNull(result.ConsensusMatrix);
        Assert.AreEqual(X.rowLength, result.ConsensusMatrix.rowLength);
        Assert.AreEqual(X.rowLength, result.ConsensusMatrix.columnLength);
    }

    [TestMethod]
    public void RunBootstrap_ConsensusMatrix_DiagonalShouldBeOne()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        // A point always co-clusters with itself
        for (int i = 0; i < X.rowLength; i++)
            Assert.AreEqual(1.0, result.ConsensusMatrix.values[i, i], 1e-10,
                $"Diagonal element [{i},{i}] should be 1.0");
    }

    [TestMethod]
    public void RunBootstrap_ConsensusMatrix_ShouldBeSymmetric()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        int n = X.rowLength;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                Assert.AreEqual(
                    result.ConsensusMatrix.values[i, j],
                    result.ConsensusMatrix.values[j, i],
                    1e-10,
                    $"Consensus[{i},{j}] should equal Consensus[{j},{i}]");
    }

    [TestMethod]
    public void RunBootstrap_ConsensusMatrix_SameClusterShouldBeHigh()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 100, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        // Points 0-9 are in cluster 0 — their consensus should be high
        double avgWithin = 0;
        int count = 0;
        for (int i = 0; i < 10; i++)
            for (int j = i + 1; j < 10; j++)
            {
                avgWithin += result.ConsensusMatrix.values[i, j];
                count++;
            }
        avgWithin /= count;

        Assert.IsTrue(avgWithin > 0.5,
            $"Average within-cluster consensus should be > 0.5, got {avgWithin:F3}");
    }

    [TestMethod]
    public void RunBootstrap_ConsensusMatrix_DifferentClustersShouldBeLow()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 100, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        // Points 0-9 (cluster 0) vs 10-19 (cluster 1) — should rarely co-cluster
        double avgBetween = 0;
        int count = 0;
        for (int i = 0; i < 10; i++)
            for (int j = 10; j < 20; j++)
            {
                avgBetween += result.ConsensusMatrix.values[i, j];
                count++;
            }
        avgBetween /= count;

        Assert.IsTrue(avgBetween < 0.2,
            $"Average between-cluster consensus should be < 0.2, got {avgBetween:F3}");
    }

    [TestMethod]
    public void RunBootstrap_ConsensusMatrix_ValuesShouldBeInZeroOne()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        int n = X.rowLength;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double val = result.ConsensusMatrix.values[i, j];
                Assert.IsTrue(val >= 0.0 && val <= 1.0,
                    $"Consensus[{i},{j}] = {val} is out of [0,1]");
            }
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunBootstrap — Point Stability
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_ShouldReturnPointStability()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        Assert.IsNotNull(result.PointStability);
        Assert.AreEqual(X.rowLength, result.PointStability.Length);

        // Well-separated clusters → high stability for all points
        foreach (double s in result.PointStability)
            Assert.IsTrue(s > 0.3, $"Expected moderate-to-high stability, got {s:F3}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunBootstrap — Convergence Curve
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_ShouldReturnConvergenceCurve()
    {
        int iterations = 60;
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = iterations, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        Assert.AreEqual(iterations, result.ConvergenceCurve.Length);

        // Last value should equal the mean
        double lastVal = result.ConvergenceCurve[iterations - 1];
        Assert.AreEqual(result.ScoreDistribution.Mean, lastVal, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunBootstrap — With Scaler
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_WithScaler_ShouldWork()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 30, Seed = 42 };

        var result = mc.RunBootstrap(
            X,
            new KMeans { K = 3 },
            new SilhouetteEvaluator(),
            new StandardScaler());

        Assert.IsNotNull(result.ScoreDistribution);
        Assert.AreEqual(30, result.ScoreDistribution.Count);
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunBootstrap — Reproducibility
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_SameSeed_ShouldBeReproducible()
    {
        var X = ThreeClusterData();

        var mc1 = new MonteCarloClustering { Iterations = 30, Seed = 99 };
        var mc2 = new MonteCarloClustering { Iterations = 30, Seed = 99 };

        var r1 = mc1.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());
        var r2 = mc2.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        Assert.AreEqual(r1.ScoreDistribution.Mean, r2.ScoreDistribution.Mean, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunExperiment — Optimal-K Distribution
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunExperiment_ShouldReturnKDistribution()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 30, Seed = 42 };

        var result = mc.RunExperiment(
            X,
            new KMeans(),
            new SilhouetteEvaluator(),
            minK: 2,
            maxK: 5);

        Assert.IsNotNull(result.OptimalKDistribution);
        Assert.IsTrue(result.OptimalKDistribution.Count > 0);

        // Sum of distribution should equal iterations
        int totalVotes = 0;
        foreach (var kvp in result.OptimalKDistribution)
            totalVotes += kvp.Value;
        Assert.AreEqual(30, totalVotes, "K distribution should sum to iteration count");
    }

    [TestMethod]
    public void RunExperiment_WellSeparated_K3ShouldDominate()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunExperiment(
            X,
            new KMeans(),
            new SilhouetteEvaluator(),
            minK: 2,
            maxK: 6);

        // K=3 should win most of the time for 3 well-separated clusters
        int k3Wins = result.OptimalKDistribution.ContainsKey(3) ? result.OptimalKDistribution[3] : 0;
        Assert.IsTrue(k3Wins > 20,
            $"K=3 should dominate for 3-cluster data, but only won {k3Wins}/{50} times");
    }

    [TestMethod]
    public void RunExperiment_ShouldReturnScoreDistribution()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 30, Seed = 42 };

        var result = mc.RunExperiment(
            X,
            new KMeans(),
            new SilhouetteEvaluator(),
            minK: 2,
            maxK: 5);

        Assert.IsNotNull(result.ScoreDistribution);
        Assert.AreEqual(30, result.ScoreDistribution.Count);
        Assert.IsTrue(result.ScoreDistribution.Mean > 0);
    }

    [TestMethod]
    public void RunExperiment_ShouldReturnConvergenceCurve()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 30, Seed = 42 };

        var result = mc.RunExperiment(
            X,
            new KMeans(),
            new SilhouetteEvaluator(),
            minK: 2,
            maxK: 4);

        Assert.IsNotNull(result.ConvergenceCurve);
        Assert.AreEqual(30, result.ConvergenceCurve.Length);
    }

    [TestMethod]
    public void RunExperiment_WithScaler_ShouldWork()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 20, Seed = 42 };

        var result = mc.RunExperiment(
            X,
            new KMeans(),
            new SilhouetteEvaluator(),
            minK: 2,
            maxK: 4,
            scaler: new StandardScaler());

        Assert.IsNotNull(result.ScoreDistribution);
        Assert.IsTrue(result.OptimalKDistribution.Values.Sum() == 20);
    }

    [TestMethod]
    public void RunExperiment_SameSeed_ShouldBeReproducible()
    {
        var X = ThreeClusterData();

        var mc1 = new MonteCarloClustering { Iterations = 20, Seed = 77 };
        var mc2 = new MonteCarloClustering { Iterations = 20, Seed = 77 };

        var r1 = mc1.RunExperiment(X, new KMeans(), new SilhouetteEvaluator(), 2, 5);
        var r2 = mc2.RunExperiment(X, new KMeans(), new SilhouetteEvaluator(), 2, 5);

        Assert.AreEqual(r1.ScoreDistribution.Mean, r2.ScoreDistribution.Mean, 1e-10);

        foreach (var k in r1.OptimalKDistribution.Keys)
            Assert.AreEqual(r1.OptimalKDistribution[k], r2.OptimalKDistribution[k]);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API Integration
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluentAPI_WithMonteCarloUncertainty_ShouldPopulateMCResult()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 5)
            .WithEvaluator(new SilhouetteEvaluator())
            .WithMonteCarloUncertainty(iterations: 30, seed: 42)
            .Run();

        // Standard result should still work
        Assert.IsNotNull(result.Best);
        Assert.IsTrue(result.Rankings.Count > 0);

        // MC result should be populated
        Assert.IsNotNull(result.MonteCarloResult, "MonteCarloResult should be populated");
        Assert.IsNotNull(result.MonteCarloResult.ScoreDistribution);
        Assert.IsNotNull(result.MonteCarloResult.OptimalKDistribution);
        Assert.AreEqual(30, result.MonteCarloResult.Iterations);
    }

    [TestMethod]
    public void FluentAPI_WithoutMonteCarlo_MCResultShouldBeNull()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 4)
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        Assert.IsNull(result.MonteCarloResult, "MonteCarloResult should be null without WithMonteCarloUncertainty");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Histogram & StandardError via MonteCarloResult
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RunBootstrap_MonteCarloResult_ShouldSupportHistogram()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 50, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        var histogram = result.ScoreDistribution.Histogram(10);
        Assert.IsTrue(histogram.Length > 0);

        int totalCount = histogram.Sum(h => h.count);
        Assert.AreEqual(50, totalCount);
    }

    [TestMethod]
    public void RunBootstrap_StandardError_ShouldEqualStdDevOverSqrtN()
    {
        var X = ThreeClusterData();
        var mc = new MonteCarloClustering { Iterations = 100, Seed = 42 };

        var result = mc.RunBootstrap(X, new KMeans { K = 3 }, new SilhouetteEvaluator());

        double se = result.ScoreDistribution.StandardError;
        double expected = result.ScoreDistribution.StandardDeviation / Math.Sqrt(100);

        Assert.IsTrue(se > 0, "Standard error should be positive");
        Assert.AreEqual(expected, se, 1e-10,
            "SE should equal StdDev / sqrt(N)");
    }
}
