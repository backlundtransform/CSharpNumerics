using CSharpNumerics.ML;
using CSharpNumerics.ML.CrossValidators;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.MonteCarlo;

namespace NumericTest;

[TestClass]
public class MonteCarloCrossValidatorTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>y = 3*x1 + 2*x2 + noise — a clean regression problem.</summary>
    private static (Matrix X, VectorN y) RegressionData(int n = 100, int seed = 123)
    {
        var rnd = new Random(seed);
        var X = new Matrix(n, 2);
        var y = new VectorN(n);

        for (int i = 0; i < n; i++)
        {
            double x1 = rnd.NextDouble() * 10;
            double x2 = rnd.NextDouble() * 10;
            X.values[i, 0] = x1;
            X.values[i, 1] = x2;
            y[i] = 3 * x1 + 2 * x2 + rnd.NextDouble() * 0.1;
        }

        return (X, y);
    }

    /// <summary>
    /// Three-class problem separated by x1+x2 thresholds.
    /// </summary>
    private static (Matrix X, VectorN y) ClassificationData(int n = 120, int seed = 123)
    {
        var rnd = new Random(seed);
        var X = new Matrix(n, 2);
        var y = new VectorN(n);

        for (int i = 0; i < n; i++)
        {
            double x1 = rnd.NextDouble() * 10;
            double x2 = rnd.NextDouble() * 10;
            X.values[i, 0] = x1;
            X.values[i, 1] = x2;
            double s = x1 + x2;
            y[i] = s < 7 ? 0 : s < 13 ? 1 : 2;
        }

        return (X, y);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Basic Run (ICrossValidator interface)
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_Regression_Run_ShouldReturnValidResult()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cv = new MonteCarloCrossValidator(
            pipelines: new List<Pipeline> { pipeline },
            iterations: 30,
            testSize: 0.2,
            seed: 42);

        CrossValidationResult result = cv.Run(X, y);

        Assert.IsNotNull(result.BestPipeline);
        Assert.IsTrue(result.BestScore < 0,
            "Regression score (neg MSE) should be negative");
    }

    [TestMethod]
    public void MCCV_Classification_Run_ShouldReturnValidResult()
    {
        var (X, y) = ClassificationData();

        var pipeline = new Pipeline(
            new DecisionTree(),
            new Dictionary<string, object> { ["MaxDepth"] = 5 });

        var cv = new MonteCarloCrossValidator(
            pipelines: new List<Pipeline> { pipeline },
            iterations: 30,
            testSize: 0.2,
            seed: 42);

        CrossValidationResult result = cv.Run(X, y);

        Assert.IsNotNull(result.BestPipeline);
        Assert.IsTrue(result.BestScore > 0.5,
            $"Expected decent accuracy, got {result.BestScore}");
        Assert.IsNotNull(result.ConfusionMatrix);
    }

    // ═══════════════════════════════════════════════════════════════
    //  RunDetailed — MonteCarloValidationResult
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_RunDetailed_ShouldReturnScoreDistribution()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cv = new MonteCarloCrossValidator(
            pipelines: new List<Pipeline> { pipeline },
            iterations: 50,
            testSize: 0.2,
            seed: 42);

        MonteCarloValidationResult detailed = cv.RunDetailed(X, y);

        // Should have one entry per pipeline
        Assert.AreEqual(1, detailed.DetailedScores.Count);

        // MonteCarloResult should contain exactly 50 samples
        MonteCarloResult mcResult = detailed.DetailedScores[pipeline];
        Assert.AreEqual(50, mcResult.Count);

        // Standard deviation should be positive (there IS variance across splits)
        Assert.IsTrue(mcResult.StandardDeviation > 0,
            "Score distribution should have positive std dev");

        // StandardError should be less than StdDev (SE = σ/√n)
        Assert.IsTrue(mcResult.StandardError < mcResult.StandardDeviation);
    }

    [TestMethod]
    public void MCCV_RunDetailed_ShouldProvideConfidenceInterval()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cv = new MonteCarloCrossValidator(
            pipelines: new List<Pipeline> { pipeline },
            iterations: 100,
            testSize: 0.2,
            seed: 42);

        MonteCarloValidationResult detailed = cv.RunDetailed(X, y);

        // CI should be a proper interval: lower < mean < upper
        var ci = detailed.BestConfidenceInterval;
        Assert.IsTrue(ci.lower < detailed.BestMeanScore,
            $"CI lower ({ci.lower}) should be less than mean ({detailed.BestMeanScore})");
        Assert.IsTrue(ci.upper > detailed.BestMeanScore,
            $"CI upper ({ci.upper}) should be greater than mean ({detailed.BestMeanScore})");

        // CI width should be positive
        Assert.IsTrue(ci.upper - ci.lower > 0, "CI width must be positive");
    }

    [TestMethod]
    public void MCCV_RunDetailed_ShouldProvideConvergenceCurve()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        int iterations = 60;
        var cv = new MonteCarloCrossValidator(
            pipelines: new List<Pipeline> { pipeline },
            iterations: iterations,
            testSize: 0.2,
            seed: 42);

        MonteCarloValidationResult detailed = cv.RunDetailed(X, y);

        // Convergence curve should have one entry per iteration
        Assert.AreEqual(iterations, detailed.ConvergenceCurve.Length);

        // Last value in convergence curve should equal the mean score
        double lastConvergence = detailed.ConvergenceCurve[iterations - 1];
        Assert.AreEqual(detailed.BestMeanScore, lastConvergence, 1e-10,
            "Last convergence value should equal the mean");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multiple pipelines
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_MultiplePipelines_ShouldRankCorrectly()
    {
        var (X, y) = ClassificationData();

        var dtPipeline = new Pipeline(
            new DecisionTree(),
            new Dictionary<string, object> { ["MaxDepth"] = 5 });

        var rfPipeline = new Pipeline(
            new RandomForest(),
            new Dictionary<string, object> { ["NumTrees"] = 30, ["MaxDepth"] = 6 });

        var cv = new MonteCarloCrossValidator(
            pipelines: new List<Pipeline> { dtPipeline, rfPipeline },
            iterations: 30,
            testSize: 0.2,
            seed: 42);

        MonteCarloValidationResult detailed = cv.RunDetailed(X, y);

        // Should have scores for both pipelines
        Assert.AreEqual(2, detailed.DetailedScores.Count);

        // Best pipeline should exist and have a valid score
        Assert.IsNotNull(detailed.BestPipeline);
        Assert.IsTrue(detailed.BestMeanScore > 0.5,
            $"Best accuracy too low: {detailed.BestMeanScore}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  PipelineGrid constructor
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_PipelineGrid_ShouldWork()
    {
        var (X, y) = ClassificationData();

        var grid = new PipelineGrid()
            .AddModel<DecisionTree>(g => g
                .Add("MaxDepth", 3, 5));

        var cv = new MonteCarloCrossValidator(
            pipelineGrid: grid,
            iterations: 20,
            testSize: 0.2,
            seed: 42);

        CrossValidationResult result = cv.Run(X, y);

        // Grid expands to 2 pipelines (MaxDepth=3 and MaxDepth=5)
        Assert.AreEqual(2, result.Scores.Count);
        Assert.IsNotNull(result.BestPipeline);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Reproducibility
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_SameSeed_ShouldProduceIdenticalResults()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cv1 = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 30, seed: 99);
        var cv2 = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 30, seed: 99);

        var r1 = cv1.RunDetailed(X, y);
        var r2 = cv2.RunDetailed(X, y);

        Assert.AreEqual(r1.BestMeanScore, r2.BestMeanScore, 1e-10,
            "Same seed should produce identical mean scores");
    }

    [TestMethod]
    public void MCCV_DifferentSeed_ShouldProduceDifferentResults()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cv1 = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 30, seed: 1);
        var cv2 = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 30, seed: 2);

        var r1 = cv1.RunDetailed(X, y);
        var r2 = cv2.RunDetailed(X, y);

        // Means will likely differ (not guaranteed but overwhelmingly likely)
        Assert.AreNotEqual(r1.BestMeanScore, r2.BestMeanScore,
            "Different seeds should (almost certainly) produce different mean scores");
    }

    // ═══════════════════════════════════════════════════════════════
    //  DetailedResult property
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_Run_ShouldPopulateDetailedResultProperty()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cv = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 20, seed: 42);

        // Before run, DetailedResult should be null
        Assert.IsNull(cv.DetailedResult);

        cv.Run(X, y);

        // After run, it should be populated
        Assert.IsNotNull(cv.DetailedResult);
        Assert.AreEqual(20, cv.DetailedResult.Iterations);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Histogram access
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_MonteCarloResult_ShouldSupportHistogram()
    {
        var (X, y) = ClassificationData();

        var pipeline = new Pipeline(
            new DecisionTree(),
            new Dictionary<string, object> { ["MaxDepth"] = 5 });

        var cv = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 50, seed: 42);

        var detailed = cv.RunDetailed(X, y);
        var mcResult = detailed.DetailedScores[pipeline];

        // Histogram should work and return bins
        var histogram = mcResult.Histogram(10);
        Assert.IsTrue(histogram.Length > 0, "Histogram should have bins");

        // Sum of all bin counts should equal iteration count
        int totalCount = histogram.Sum(h => h.count);
        Assert.AreEqual(50, totalCount, "Histogram bin counts should sum to iterations");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Standard error
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MCCV_MoreIterations_ShouldReduceStandardError()
    {
        var (X, y) = RegressionData();

        var pipeline = new Pipeline(
            new Linear(),
            new Dictionary<string, object> { ["LearningRate"] = 0.01 },
            scaler: new StandardScaler());

        var cvSmall = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 20, seed: 42);
        var cvLarge = new MonteCarloCrossValidator(
            new List<Pipeline> { pipeline }, iterations: 200, seed: 42);

        var rSmall = cvSmall.RunDetailed(X, y);
        var rLarge = cvLarge.RunDetailed(X, y);

        double seSmall = rSmall.DetailedScores[pipeline].StandardError;
        double seLarge = rLarge.DetailedScores[pipeline].StandardError;

        Assert.IsTrue(seLarge < seSmall,
            $"More iterations should reduce SE: SE(20)={seSmall:F6}, SE(200)={seLarge:F6}");
    }
}
