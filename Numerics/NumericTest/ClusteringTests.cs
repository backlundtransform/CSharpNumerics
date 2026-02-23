using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class ClusteringTests
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

    /// <summary>Two clearly separated 2D clusters.</summary>
    private static Matrix TwoClusterData()
    {
        var data = new double[,]
        {
            { 0.0, 0.0 }, { 0.5, 0.5 }, { 1.0, 0.0 }, { 0.5, -0.5 }, { 0.0, 1.0 },
            { 10.0, 10.0 }, { 10.5, 10.5 }, { 11.0, 10.0 }, { 10.5, 9.5 }, { 10.0, 11.0 },
        };
        return new Matrix(data);
    }

    // ═══════════════════════════════════════════════════════════════
    //  KMeans
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void KMeans_ThreeClusters_ShouldFindThreeGroups()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42 };

        VectorN labels = kmeans.FitPredict(X);

        Assert.AreEqual(3, kmeans.ClusterCount);
        Assert.AreEqual(X.rowLength, labels.Length);

        // All points in first 10 rows should have same label
        int label0 = (int)labels[0];
        for (int i = 1; i < 10; i++)
            Assert.AreEqual(label0, (int)labels[i], $"Point {i} should be in same cluster as point 0");

        // All points in rows 10-19 should share label, different from cluster 0
        int label1 = (int)labels[10];
        Assert.AreNotEqual(label0, label1);
        for (int i = 11; i < 20; i++)
            Assert.AreEqual(label1, (int)labels[i], $"Point {i} should be in same cluster as point 10");

        // All points in rows 20-29 should share label, different from both
        int label2 = (int)labels[20];
        Assert.AreNotEqual(label0, label2);
        Assert.AreNotEqual(label1, label2);
        for (int i = 21; i < 30; i++)
            Assert.AreEqual(label2, (int)labels[i], $"Point {i} should be in same cluster as point 20");
    }

    [TestMethod]
    public void KMeans_PlusPlus_ShouldConverge()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42, InitMethod = KMeansInit.PlusPlus };

        kmeans.FitPredict(X);

        Assert.IsTrue(kmeans.Iterations > 0);
        Assert.IsTrue(kmeans.Iterations <= kmeans.MaxIterations);
        Assert.IsTrue(kmeans.Inertia >= 0);
        Assert.AreEqual(3, kmeans.Centroids.rowLength);
        Assert.AreEqual(2, kmeans.Centroids.columnLength);
    }

    [TestMethod]
    public void KMeans_Predict_AfterFit_ShouldAssignNewPoints()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42 };
        kmeans.Fit(X);

        // New point near cluster 1 center (10,10)
        var newPoint = new Matrix(new double[,] { { 10.1, 10.1 } });
        VectorN pred = kmeans.Predict(newPoint);

        // Should get same label as training point 10
        VectorN trainLabels = kmeans.Predict(X);
        Assert.AreEqual((int)trainLabels[10], (int)pred[0]);
    }

    [TestMethod]
    public void KMeans_Clone_ShouldPreserveHyperparameters()
    {
        var original = new KMeans { K = 5, MaxIterations = 100, Tolerance = 1e-6, Seed = 7 };
        var clone = (KMeans)original.Clone();

        Assert.AreEqual(5, clone.K);
        Assert.AreEqual(100, clone.MaxIterations);
        Assert.AreEqual(1e-6, clone.Tolerance);
        Assert.AreEqual(7, clone.Seed);
    }

    [TestMethod]
    public void KMeans_SetHyperParameters_ShouldApply()
    {
        var kmeans = new KMeans();
        kmeans.SetHyperParameters(new Dictionary<string, object>
        {
            ["K"] = 7,
            ["MaxIterations"] = 50
        });

        Assert.AreEqual(7, kmeans.K);
        Assert.AreEqual(50, kmeans.MaxIterations);
    }

    // ═══════════════════════════════════════════════════════════════
    //  DBSCAN
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void DBSCAN_ThreeClusters_ShouldDiscoverThreeGroups()
    {
        var X = ThreeClusterData();
        var dbscan = new DBSCAN { Epsilon = 2.0, MinPoints = 3 };

        VectorN labels = dbscan.FitPredict(X);

        Assert.AreEqual(3, dbscan.ClusterCount);
        Assert.AreEqual(0, dbscan.NoiseCount);
    }

    [TestMethod]
    public void DBSCAN_TightEpsilon_ShouldProduceNoise()
    {
        var X = ThreeClusterData();
        var dbscan = new DBSCAN { Epsilon = 0.01, MinPoints = 5 };

        VectorN labels = dbscan.FitPredict(X);

        // With very tight epsilon, most points become noise
        Assert.IsTrue(dbscan.NoiseCount > 0, "Expected some noise points with tight epsilon");
    }

    [TestMethod]
    public void DBSCAN_Clone_ShouldPreserveHyperparameters()
    {
        var original = new DBSCAN { Epsilon = 1.5, MinPoints = 10 };
        var clone = (DBSCAN)original.Clone();

        Assert.AreEqual(1.5, clone.Epsilon);
        Assert.AreEqual(10, clone.MinPoints);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Agglomerative Clustering
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Agglomerative_Ward_ThreeClusters_ShouldFindThreeGroups()
    {
        var X = ThreeClusterData();
        var agg = new AgglomerativeClustering { K = 3, Linkage = LinkageType.Ward };

        VectorN labels = agg.FitPredict(X);

        Assert.AreEqual(3, agg.ClusterCount);
        Assert.AreEqual(X.rowLength, labels.Length);

        // Verify cluster consistency within known groups
        int label0 = (int)labels[0];
        for (int i = 1; i < 10; i++)
            Assert.AreEqual(label0, (int)labels[i]);

        int label1 = (int)labels[10];
        Assert.AreNotEqual(label0, label1);
        for (int i = 11; i < 20; i++)
            Assert.AreEqual(label1, (int)labels[i]);
    }

    [TestMethod]
    public void Agglomerative_SingleLinkage_ShouldWork()
    {
        var X = TwoClusterData();
        var agg = new AgglomerativeClustering { K = 2, Linkage = LinkageType.Single };

        VectorN labels = agg.FitPredict(X);

        Assert.AreEqual(2, agg.ClusterCount);
        Assert.IsNotNull(agg.Dendrogram);
        Assert.IsTrue(agg.Dendrogram.Count > 0);
    }

    [TestMethod]
    public void Agglomerative_CompleteLinkage_ShouldWork()
    {
        var X = TwoClusterData();
        var agg = new AgglomerativeClustering { K = 2, Linkage = LinkageType.Complete };

        VectorN labels = agg.FitPredict(X);
        Assert.AreEqual(2, agg.ClusterCount);
    }

    [TestMethod]
    public void Agglomerative_AverageLinkage_ShouldWork()
    {
        var X = TwoClusterData();
        var agg = new AgglomerativeClustering { K = 2, Linkage = LinkageType.Average };

        VectorN labels = agg.FitPredict(X);
        Assert.AreEqual(2, agg.ClusterCount);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Evaluators
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void SilhouetteEvaluator_WellSeparatedClusters_ShouldBeHigh()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42 };
        VectorN labels = kmeans.FitPredict(X);

        var evaluator = new SilhouetteEvaluator();
        double score = evaluator.Score(X, labels);

        Assert.IsTrue(score > 0.7, $"Expected silhouette > 0.7 for well-separated clusters, got {score:F3}");
    }

    [TestMethod]
    public void SilhouetteEvaluator_SingleCluster_ShouldReturnNegative()
    {
        var X = ThreeClusterData();
        var labels = new VectorN(X.rowLength);
        // All same label
        for (int i = 0; i < X.rowLength; i++)
            labels[i] = 0;

        var evaluator = new SilhouetteEvaluator();
        double score = evaluator.Score(X, labels);

        Assert.AreEqual(-1, score, "Single cluster should return -1");
    }

    [TestMethod]
    public void InertiaEvaluator_ShouldReturnNegativeValue()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42 };
        VectorN labels = kmeans.FitPredict(X);

        var evaluator = new InertiaEvaluator();
        double score = evaluator.Score(X, labels);
        double raw = evaluator.RawInertia(X, labels);

        Assert.IsTrue(score < 0, "Inertia score should be negative (negated convention)");
        Assert.IsTrue(raw > 0, "Raw inertia should be positive");
        Assert.AreEqual(-raw, score, 1e-10, "Score should be negated raw inertia");
    }

    [TestMethod]
    public void InertiaEvaluator_MoreClusters_ShouldHaveLowerRawInertia()
    {
        var X = ThreeClusterData();
        var evaluator = new InertiaEvaluator();

        var kmeans2 = new KMeans { K = 2, Seed = 42 };
        var labels2 = kmeans2.FitPredict(X);
        double inertia2 = evaluator.RawInertia(X, labels2);

        var kmeans3 = new KMeans { K = 3, Seed = 42 };
        var labels3 = kmeans3.FitPredict(X);
        double inertia3 = evaluator.RawInertia(X, labels3);

        Assert.IsTrue(inertia3 < inertia2, "More clusters should reduce inertia");
    }

    [TestMethod]
    public void DaviesBouldinEvaluator_WellSeparated_ShouldBeNearZero()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42 };
        VectorN labels = kmeans.FitPredict(X);

        var evaluator = new DaviesBouldinEvaluator();
        double score = evaluator.Score(X, labels);

        // Negated: well-separated → low DB → high score (close to 0 but negative)
        Assert.IsTrue(score < 0, "DB score should be negative (negated)");
        Assert.IsTrue(score > -1.0, $"Expected negated DB close to 0 for separated clusters, got {score:F3}");
    }

    [TestMethod]
    public void CalinskiHarabaszEvaluator_WellSeparated_ShouldBeHigh()
    {
        var X = ThreeClusterData();
        var kmeans = new KMeans { K = 3, Seed = 42 };
        VectorN labels = kmeans.FitPredict(X);

        var evaluator = new CalinskiHarabaszEvaluator();
        double score = evaluator.Score(X, labels);

        Assert.IsTrue(score > 100, $"Expected high CH for well-separated clusters, got {score:F1}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  ClusteringPipeline
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ClusteringPipeline_FitPredict_ShouldWork()
    {
        var X = ThreeClusterData();
        var pipeline = new ClusteringPipeline(
            model: new KMeans(),
            modelParams: new Dictionary<string, object> { ["K"] = 3, ["Seed"] = 42 }
        );

        VectorN labels = pipeline.FitPredict(X);

        Assert.AreEqual(X.rowLength, labels.Length);
    }

    [TestMethod]
    public void ClusteringPipeline_WithScaler_ShouldWork()
    {
        var X = ThreeClusterData();
        var pipeline = new ClusteringPipeline(
            model: new KMeans(),
            modelParams: new Dictionary<string, object> { ["K"] = 3, ["Seed"] = 42 },
            scaler: new StandardScaler()
        );

        VectorN labels = pipeline.FitPredict(X);

        Assert.AreEqual(X.rowLength, labels.Length);
    }

    [TestMethod]
    public void ClusteringPipeline_Clone_ShouldBeIndependent()
    {
        var pipeline = new ClusteringPipeline(
            model: new KMeans(),
            modelParams: new Dictionary<string, object> { ["K"] = 3 },
            scaler: new StandardScaler()
        );

        var clone = pipeline.Clone();

        Assert.IsNotNull(clone);
        Assert.AreNotSame(pipeline.Model, clone.Model);
    }

    // ═══════════════════════════════════════════════════════════════
    //  ClusteringGrid
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ClusteringGrid_Expand_ShouldProduceCorrectCount()
    {
        var grid = new ClusteringGrid()
            .AddModel<KMeans>(g => g
                .Add("K", 2, 3, 4)
                .Add("Seed", 42))
            .AddModel<AgglomerativeClustering>(g => g
                .Add("K", 2, 3));

        var pipelines = grid.Expand().ToList();

        // KMeans: 3 K values × 1 Seed = 3
        // Agglomerative: 2 K values = 2
        Assert.AreEqual(5, pipelines.Count);
    }

    [TestMethod]
    public void ClusteringGrid_WithScaler_ShouldExpand()
    {
        var grid = new ClusteringGrid()
            .AddModel<KMeans>(g => g
                .Add("K", 2, 3)
                .AddScaler<StandardScaler>(s => { }));

        var pipelines = grid.Expand().ToList();

        Assert.AreEqual(2, pipelines.Count);
        Assert.IsNotNull(pipelines[0].Scaler);
    }

    // ═══════════════════════════════════════════════════════════════
    //  ClusteringExperiment (Fluent API)
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ClusteringExperiment_SimpleKMeans_ShouldFindBestK()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 5)
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        Assert.IsNotNull(result);
        Assert.IsNotNull(result.Best);
        Assert.AreEqual(3, result.BestClusterCount, "Should find K=3 as best for 3 well-separated clusters");
        Assert.IsTrue(result.BestScore > 0.5);
        Assert.AreEqual(4, result.Rankings.Count, "Should have results for K=2,3,4,5");
    }

    [TestMethod]
    public void ClusteringExperiment_WithScaler_ShouldWork()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 4)
            .WithEvaluator(new SilhouetteEvaluator())
            .WithScaler(new StandardScaler())
            .Run();

        Assert.IsNotNull(result.Best);
        Assert.IsTrue(result.Rankings.Count == 3);
    }

    [TestMethod]
    public void ClusteringExperiment_MultipleEvaluators_ShouldScoreAll()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 4)
            .WithEvaluators(
                new SilhouetteEvaluator(),
                new CalinskiHarabaszEvaluator(),
                new DaviesBouldinEvaluator())
            .Run();

        // Each result should have all 3 evaluator scores
        foreach (var r in result.Rankings)
        {
            Assert.IsTrue(r.Scores.ContainsKey("Silhouette"));
            Assert.IsTrue(r.Scores.ContainsKey("CalinskiHarabasz"));
            Assert.IsTrue(r.Scores.ContainsKey("DaviesBouldin"));
        }

        // Primary is Silhouette (first added)
        Assert.AreEqual("Silhouette", result.PrimaryEvaluator);
    }

    [TestMethod]
    public void ClusteringExperiment_BestByEvaluator_ShouldWork()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 5)
            .WithEvaluators(new SilhouetteEvaluator(), new CalinskiHarabaszEvaluator())
            .Run();

        var bestSilhouette = result.BestBy<SilhouetteEvaluator>();
        var bestCH = result.BestBy<CalinskiHarabaszEvaluator>();

        Assert.IsNotNull(bestSilhouette);
        Assert.IsNotNull(bestCH);
    }

    [TestMethod]
    public void ClusteringExperiment_WithInertia_ShouldPopulateElbowCurve()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 6)
            .WithEvaluators(new SilhouetteEvaluator(), new InertiaEvaluator())
            .Run();

        Assert.IsTrue(result.ElbowCurve.Count == 5, "Should have elbow data for K=2..6");
        // Inertia should decrease as K increases
        for (int i = 1; i < result.ElbowCurve.Count; i++)
            Assert.IsTrue(result.ElbowCurve[i].Inertia <= result.ElbowCurve[i - 1].Inertia,
                "Inertia should decrease with more clusters");
    }

    [TestMethod]
    public void ClusteringExperiment_MultiAlgorithm_ShouldCompareAll()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithms(
                new KMeans { Seed = 42 },
                new AgglomerativeClustering())
            .TryClusterCounts(2, 4)
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        // 3 K values × 2 algorithms = 6 results
        Assert.AreEqual(6, result.Rankings.Count);

        // Should have both algorithm names
        var names = result.Rankings.Select(r => r.AlgorithmName).Distinct().ToList();
        Assert.IsTrue(names.Contains("KMeans"));
        Assert.IsTrue(names.Contains("AgglomerativeClustering"));
    }

    [TestMethod]
    public void ClusteringExperiment_DBSCAN_ShouldIgnoreClusterRange()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new DBSCAN { Epsilon = 2.0, MinPoints = 3 })
            .TryClusterCounts(2, 10)  // should be ignored for DBSCAN
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        // DBSCAN should produce exactly 1 result (not 9)
        Assert.AreEqual(1, result.Rankings.Count);
        Assert.AreEqual(3, result.BestClusterCount);
    }

    [TestMethod]
    public void ClusteringExperiment_WithGrid_ShouldExpandAll()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithGrid(new ClusteringGrid()
                .AddModel<KMeans>(g => g
                    .Add("K", 2, 3)
                    .Add("Seed", 42)))
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        Assert.AreEqual(2, result.Rankings.Count);
    }

    [TestMethod]
    public void ClusteringExperiment_ResultHasDuration()
    {
        var X = ThreeClusterData();

        var result = ClusteringExperiment
            .For(X)
            .WithAlgorithm(new KMeans { Seed = 42 })
            .TryClusterCounts(2, 3)
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        Assert.IsTrue(result.TotalDuration.TotalMilliseconds > 0);
        foreach (var r in result.Rankings)
            Assert.IsTrue(r.Duration.TotalMilliseconds >= 0);
    }
}
