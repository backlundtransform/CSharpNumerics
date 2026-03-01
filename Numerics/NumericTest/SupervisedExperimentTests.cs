using CSharpNumerics.ML;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Experiment;
using CSharpNumerics.ML.Experiment.Results;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.ML.Selector;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class SupervisedExperimentTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Test data helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Perfect linear data: y = 2x + 1</summary>
    private static (Matrix X, VectorN y) LinearData()
    {
        double[,] Xdata =
        {
            { 0 }, { 1 }, { 2 }, { 3 }, { 4 },
            { 5 }, { 6 }, { 7 }, { 8 }, { 9 }
        };
        double[] ydata = { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19 };
        return (new Matrix(Xdata), new VectorN(ydata));
    }

    /// <summary>Binary classification: y = x > 4.5</summary>
    private static (Matrix X, VectorN y) BinaryClassificationData()
    {
        double[,] Xdata =
        {
            { 0 }, { 1 }, { 2 }, { 3 }, { 4 },
            { 5 }, { 6 }, { 7 }, { 8 }, { 9 }
        };
        double[] ydata = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };
        return (new Matrix(Xdata), new VectorN(ydata));
    }

    /// <summary>Multi-class data: 3 classes based on x1 + x2 thresholds.</summary>
    private static (Matrix X, VectorN y) MultiClassData(int n = 120)
    {
        var X = new Matrix(n, 2);
        var y = new VectorN(n);
        var rnd = new Random(42);

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
    //  Simple pipeline mode
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void SimpleRegression_SinglePipeline_KFold_ShouldWork()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithPipeline(new Pipeline(
                model: new Linear(),
                modelParams: new Dictionary<string, object>()))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        Assert.IsNotNull(result);
        Assert.IsNotNull(result.Best);
        Assert.AreEqual(1, result.Rankings.Count);
        Assert.AreEqual("Linear", result.BestModelName);
        Assert.AreEqual("KFold", result.PrimaryCrossValidator);
        Assert.IsTrue(result.BestScore > -1.0,
            $"Expected near-zero MSE for linear data, got {result.BestScore}");
    }

    [TestMethod]
    public void SimpleClassification_SinglePipeline_KFold_ShouldWork()
    {
        var (X, y) = BinaryClassificationData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithPipeline(new Pipeline(
                model: new Logistic(),
                modelParams: new Dictionary<string, object>
                {
                    { "LearningRate", 0.1 },
                    { "MaxIterations", 2000 }
                }))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        Assert.IsNotNull(result.Best);
        Assert.AreEqual("Logistic", result.BestModelName);
    }

    [TestMethod]
    public void MultiplePipelines_ShouldRankCorrectly()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithPipelines(
                new Pipeline(new Linear(), new Dictionary<string, object>()),
                new Pipeline(new Ridge(), new Dictionary<string, object> { { "Alpha", 0.1 } }))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        Assert.AreEqual(2, result.Rankings.Count);
        Assert.IsNotNull(result.Best);
        // First in ranking should have the highest (least negative) score
        Assert.IsTrue(result.Rankings[0].Scores["KFold"] >= result.Rankings[1].Scores["KFold"]);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Grid mode
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void GridRegression_ShouldExpandAndRank()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<Linear>(g => { })
                .AddModel<Ridge>(g => g
                    .Add("Alpha", 0.01, 0.1, 1.0)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        // Linear: 1 + Ridge: 3 Alpha values = 4 pipelines
        Assert.AreEqual(4, result.Rankings.Count);
        Assert.IsNotNull(result.Best);
        Assert.IsTrue(result.BestScore > -10.0);
    }

    [TestMethod]
    public void GridClassification_MultiModel_ShouldWork()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)
                    .AddScaler<StandardScaler>(s => { }))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        // KNN: 3 K values = 3, DecisionTree: 3 MaxDepth = 3 → 6 total
        Assert.AreEqual(6, result.Rankings.Count);

        var names = result.Rankings.Select(r => r.ModelName).Distinct().ToList();
        Assert.IsTrue(names.Contains("KNearestNeighbors"));
        Assert.IsTrue(names.Contains("DecisionTree"));

        Assert.IsTrue(result.BestScore > 0.5,
            $"Expected reasonable accuracy, got {result.BestScore}");
    }

    [TestMethod]
    public void GridWithScalerAndSelector_ShouldExpandAll()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<Logistic>(g => g
                    .Add("LearningRate", 0.05, 0.1)
                    .Add("MaxIterations", 1000)
                    .AddScaler<StandardScaler>(s => { })
                    .AddSelector<SelectKBest>(s => s
                        .Add("K", 1, 2))))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        // 2 LearningRate × 1 MaxIter × 1 Scaler × 2 SelectorK = 4
        Assert.AreEqual(4, result.Rankings.Count);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multiple cross-validators
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MultipleCrossValidators_ShouldScoreAll()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 5)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        // Each result should have scores from both CVs
        foreach (var r in result.Rankings)
        {
            Assert.IsTrue(r.Scores.ContainsKey("KFold"),
                $"Missing KFold score for {r.PipelineDescription}");
            Assert.IsTrue(r.Scores.ContainsKey("ShuffleSplit"),
                $"Missing ShuffleSplit score for {r.PipelineDescription}");
        }

        // Primary is KFold (first added)
        Assert.AreEqual("KFold", result.PrimaryCrossValidator);
    }

    [TestMethod]
    public void MultipleCrossValidators_Classification_KFoldAndStratified()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.StratifiedKFold(folds: 5))
            .Run();

        Assert.AreEqual(2, result.Rankings.Count);

        foreach (var r in result.Rankings)
        {
            Assert.IsTrue(r.Scores.ContainsKey("KFold"));
            Assert.IsTrue(r.Scores.ContainsKey("StratifiedKFold"));
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  BestBy cross-validator
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void BestBy_ShouldReturnBestForSpecificCV()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        var bestKFold = result.BestBy("KFold");
        var bestShuffle = result.BestBy("ShuffleSplit");

        Assert.IsNotNull(bestKFold);
        Assert.IsNotNull(bestShuffle);

        // Best by KFold should have the highest KFold score
        foreach (var r in result.Rankings)
        {
            if (r.Scores.ContainsKey("KFold"))
                Assert.IsTrue(bestKFold.Scores["KFold"] >= r.Scores["KFold"]);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  CVResults property
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CVResults_ShouldContainDetailedResults()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 5)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5))
            .Run();

        Assert.IsTrue(result.CVResults.ContainsKey("KFold"));
        var cvResult = result.CVResults["KFold"];

        Assert.IsNotNull(cvResult.BestPipeline);
        Assert.IsTrue(cvResult.BestScore > 0);
        Assert.IsNotNull(cvResult.ConfusionMatrix);
    }

    [TestMethod]
    public void BestConfusionMatrix_Classification_ShouldBePopulated()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 5)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        Assert.IsNotNull(result.BestConfusionMatrix);
        Assert.AreEqual(3, result.BestConfusionMatrix.rowLength);
        Assert.AreEqual(3, result.BestConfusionMatrix.columnLength);
    }

    [TestMethod]
    public void BestR2_Regression_ShouldBePopulated()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithPipeline(new Pipeline(new Linear(), new Dictionary<string, object>()))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        // Perfect linear data should give high R²
        Assert.IsTrue(result.BestR2 > 0.9,
            $"Expected high R² for linear data, got {result.BestR2}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Pipeline description
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PipelineDescription_ShouldIncludeModelAndParams()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 5)
                    .AddScaler<StandardScaler>(s => { })))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        var desc = result.Best.PipelineDescription;
        Assert.IsTrue(desc.Contains("KNearestNeighbors"), $"Description should contain model name: {desc}");
        Assert.IsTrue(desc.Contains("StandardScaler"), $"Description should contain scaler name: {desc}");
        Assert.IsTrue(desc.Contains("K=5"), $"Description should contain params: {desc}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Duration tracking
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Duration_ShouldBeTracked()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithPipeline(new Pipeline(new Linear(), new Dictionary<string, object>()))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        Assert.IsTrue(result.TotalDuration.TotalMilliseconds > 0);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Error handling
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void NoCrossValidator_ShouldThrow()
    {
        var (X, y) = LinearData();

        SupervisedExperiment
            .For(X, y)
            .WithPipeline(new Pipeline(new Linear(), new Dictionary<string, object>()))
            .Run();
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void NoPipelines_ShouldThrow()
    {
        var (X, y) = LinearData();

        SupervisedExperiment
            .For(X, y)
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();
    }

    // ═══════════════════════════════════════════════════════════════
    //  Full kitchen-sink experiment
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FullExperiment_MultiModel_MultiCV_ShouldWork()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)
                    .AddScaler<StandardScaler>(s => { }))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10))
                .AddModel<RandomForest>(g => g
                    .Add("NumTrees", 50)
                    .Add("MaxDepth", 5, 8)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        // KNN: 3, DT: 3, RF: 2 = 8 pipelines
        Assert.AreEqual(8, result.Rankings.Count);

        // All pipelines should have scores from both CVs
        foreach (var r in result.Rankings)
        {
            Assert.AreEqual(2, r.Scores.Count, $"Expected 2 CV scores for {r.PipelineDescription}");
            Assert.IsTrue(r.Scores.ContainsKey("KFold"));
            Assert.IsTrue(r.Scores.ContainsKey("ShuffleSplit"));
        }

        // Model names should include all three types
        var modelNames = result.Rankings.Select(r => r.ModelName).Distinct().ToList();
        Assert.IsTrue(modelNames.Contains("KNearestNeighbors"));
        Assert.IsTrue(modelNames.Contains("DecisionTree"));
        Assert.IsTrue(modelNames.Contains("RandomForest"));

        // Best should be reasonable
        Assert.IsTrue(result.BestScore > 0.5);
        Assert.IsNotNull(result.BestPipeline);

        // CVResults should have full details
        Assert.AreEqual(2, result.CVResults.Count);
        Assert.IsNotNull(result.CVResults["KFold"].BestPipeline);
        Assert.IsNotNull(result.CVResults["ShuffleSplit"].BestPipeline);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Mixed regression grid
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FullExperiment_Regression_MultiModel_ShouldWork()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<Linear>(g => { })
                .AddModel<Ridge>(g => g
                    .Add("Alpha", 0.01, 0.1, 1.0))
                .AddModel<Lasso>(g => g
                    .Add("Alpha", 0.01, 0.1)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        // Linear: 1, Ridge: 3, Lasso: 2 = 6
        Assert.AreEqual(6, result.Rankings.Count);
        Assert.IsTrue(result.BestScore > -10.0,
            $"Expected reasonable regression score, got {result.BestScore}");

        var modelNames = result.Rankings.Select(r => r.ModelName).Distinct().ToList();
        Assert.IsTrue(modelNames.Contains("Linear"));
        Assert.IsTrue(modelNames.Contains("Ridge"));
        Assert.IsTrue(modelNames.Contains("Lasso"));
    }

    // ═══════════════════════════════════════════════════════════════
    //  CrossValidatorConfig.Custom
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CustomCrossValidator_ShouldWork()
    {
        var (X, y) = LinearData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithPipeline(new Pipeline(new Linear(), new Dictionary<string, object>()))
            .WithCrossValidator(
                CrossValidatorConfig.Custom("MyKFold3",
                    p => new CSharpNumerics.ML.CrossValidators.KFoldCrossValidator(p, folds: 3)))
            .Run();

        Assert.IsNotNull(result.Best);
        Assert.IsTrue(result.Best.Scores.ContainsKey("MyKFold3"));
    }

    // ═══════════════════════════════════════════════════════════════
    //  Score Distribution Statistics
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ScoreSummary_Classification_ShouldComputeDescriptiveStats()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)
                    .AddScaler<StandardScaler>(s => { }))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        var summary = result.ScoreSummary();

        Assert.AreEqual(6, summary.Count);
        Assert.IsTrue(summary.Mean > 0, "Mean accuracy should be positive");
        Assert.IsTrue(summary.Median > 0);
        Assert.IsTrue(summary.StandardDeviation >= 0);
        Assert.IsTrue(summary.Range >= 0);
        Assert.IsTrue(summary.InterquartileRange >= 0);
        Assert.IsTrue(summary.Max >= summary.Min);
        Assert.IsTrue(!double.IsNaN(summary.Skewness));
        Assert.IsTrue(!double.IsNaN(summary.Kurtosis));
    }

    [TestMethod]
    public void ScoreSummary_SpecificCV_ShouldWork()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        var kfSummary = result.ScoreSummary("KFold");
        var ssSummary = result.ScoreSummary("ShuffleSplit");

        Assert.AreEqual(3, kfSummary.Count);
        Assert.AreEqual(3, ssSummary.Count);
    }

    [TestMethod]
    public void ScorePercentile_BestShouldBeHighest()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)
                    .AddScaler<StandardScaler>(s => { }))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        double bestPct = result.ScorePercentile(result.Best);
        double worstPct = result.ScorePercentile(result.Rankings.Last());

        Assert.IsTrue(bestPct >= worstPct,
            $"Best ({bestPct:F1}) should be >= worst ({worstPct:F1})");
    }

    [TestMethod]
    public void RankCorrelation_TwoCVs_ShouldReturnValidRho()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)
                    .AddScaler<StandardScaler>(s => { }))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        var (rho, pValue) = result.RankCorrelation("KFold", "ShuffleSplit");

        Assert.IsTrue(rho >= -1.0 && rho <= 1.0,
            $"Spearman rho should be in [-1, 1], got {rho:F4}");
        Assert.IsTrue(pValue >= 0 && pValue <= 1.0,
            $"p-value should be in [0, 1], got {pValue:F4}");
        // For similar CV strategies on same data, ranking should be somewhat correlated
        Assert.IsTrue(rho > -0.5,
            $"Expected at least moderate agreement between KFold and ShuffleSplit, got rho={rho:F4}");
    }

    [TestMethod]
    public void ScoreConsistency_SingleCV_ShouldBeZero()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 5)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        // With only one CV, consistency stddev should be 0
        double consistency = result.ScoreConsistency(result.Best);
        Assert.AreEqual(0.0, consistency, 1e-10);
    }

    [TestMethod]
    public void ScoreConsistency_MultipleCVs_ShouldBeNonNegative()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 5)))
            .WithCrossValidators(
                CrossValidatorConfig.KFold(folds: 5),
                CrossValidatorConfig.ShuffleSplit(nSplits: 5, testSize: 0.2, randomState: 42))
            .Run();

        double consistency = result.ScoreConsistency(result.Best);
        Assert.IsTrue(consistency >= 0, "Consistency should be non-negative");
    }

    [TestMethod]
    public void ScoreSummary_ConfidenceInterval_ShouldBracketMean()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7)
                    .AddScaler<StandardScaler>(s => { }))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        var summary = result.ScoreSummary();
        var ci = summary.ConfidenceInterval;

        Assert.IsTrue(ci.Lower <= summary.Mean,
            $"CI lower ({ci.Lower:F4}) should be <= mean ({summary.Mean:F4})");
        Assert.IsTrue(ci.Upper >= summary.Mean,
            $"CI upper ({ci.Upper:F4}) should be >= mean ({summary.Mean:F4})");
    }

    [TestMethod]
    public void ScoreSummary_ToString_ShouldBeReadable()
    {
        var (X, y) = MultiClassData();

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<KNearestNeighbors>(g => g
                    .Add("K", 3, 5, 7, 9))
                .AddModel<DecisionTree>(g => g
                    .Add("MaxDepth", 3, 5, 10)))
            .WithCrossValidator(CrossValidatorConfig.KFold(folds: 5))
            .Run();

        string text = result.ScoreSummary().ToString();
        Assert.IsTrue(text.Contains("N="), $"ToString should contain 'N=': {text}");
        Assert.IsTrue(text.Contains("Mean="), $"ToString should contain 'Mean=': {text}");
        Assert.IsTrue(text.Contains("IQR="), $"ToString should contain 'IQR=': {text}");
    }
}
