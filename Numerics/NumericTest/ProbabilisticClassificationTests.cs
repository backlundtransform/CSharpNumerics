using CSharpNumerics.ML;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Evaluation;
using CSharpNumerics.ML.Models;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.ML.Sequence.Models.Classification;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericTest;

[TestClass]
public class ProbabilisticClassificationTests
{
    #region Data

    private static double Gaussian(Random rng) =>
        Math.Sqrt(-2.0 * Math.Log(1.0 - rng.NextDouble())) * Math.Cos(2.0 * Math.PI * rng.NextDouble());

    /// <summary>Two Gaussian blobs in 2-D: class 0 around (0, 0), class 1 around (d, d).</summary>
    private static (Matrix X, VectorN y) Blobs(int negatives, int positives, double distance, int seed)
    {
        var rng = new Random(seed);
        int n = negatives + positives;
        var X = new double[n, 2];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            bool positive = i >= negatives;
            double centre = positive ? distance : 0.0;
            X[i, 0] = centre + Gaussian(rng);
            X[i, 1] = centre + Gaussian(rng);
            y[i] = positive ? 1 : 0;
        }
        return (new Matrix(X), new VectorN(y));
    }

    private static void AssertRowsAreProbabilities(Matrix p)
    {
        for (int i = 0; i < p.rowLength; i++)
        {
            double sum = 0;
            for (int c = 0; c < p.columnLength; c++)
            {
                Assert.IsTrue(p.values[i, c] >= 0 && p.values[i, c] <= 1, $"p[{i},{c}] = {p.values[i, c]}");
                sum += p.values[i, c];
            }
            Assert.AreEqual(1.0, sum, 1e-9, $"row {i}");
        }
    }

    private static void AssertArgmaxMatchesPredict(IProbabilisticClassifier model, Matrix X)
    {
        var p = model.PredictProbabilities(X);
        var predicted = model.Predict(X);
        for (int i = 0; i < X.rowLength; i++)
        {
            int best = 0;
            for (int c = 1; c < p.columnLength; c++)
                if (p.values[i, c] > p.values[i, best]) best = c;
            Assert.AreEqual(best, predicted[i], $"row {i}");
        }
    }

    private static double Recall(VectorN y, VectorN predicted)
    {
        int positives = 0, found = 0;
        for (int i = 0; i < y.Length; i++)
        {
            if (y[i] != 1) continue;
            positives++;
            if (predicted[i] == 1) found++;
        }
        return (double)found / positives;
    }

    #endregion

    #region DecisionTree

    [TestMethod]
    public void DecisionTree_Probabilities_SumToOneAndAgreeWithPredict()
    {
        var (X, y) = Blobs(80, 80, 2.0, seed: 1);
        var tree = new DecisionTree();
        tree.SetHyperParameters(new() { ["MaxDepth"] = 4 });
        tree.Fit(X, y);

        var (Xtest, _) = Blobs(40, 40, 2.0, seed: 2);
        var p = tree.PredictProbabilities(Xtest);
        Assert.AreEqual(2, p.columnLength);
        AssertRowsAreProbabilities(p);
        AssertArgmaxMatchesPredict(tree, Xtest);
    }

    [TestMethod]
    public void DecisionTree_EqualWeights_GiveTheSameTreeAsNoWeights()
    {
        var (X, y) = Blobs(60, 60, 1.5, seed: 3);
        var (Xtest, _) = Blobs(30, 30, 1.5, seed: 4);

        var plain = new DecisionTree();
        plain.Fit(X, y);
        var weighted = new DecisionTree();
        // Doubling is exact in floating point, so both trees see bit-identical impurities.
        weighted.Fit(X, y, new VectorN(Enumerable.Repeat(2.0, y.Length).ToArray()));

        var a = plain.PredictProbabilities(Xtest);
        var b = weighted.PredictProbabilities(Xtest);
        for (int i = 0; i < Xtest.rowLength; i++)
            Assert.AreEqual(a.values[i, 1], b.values[i, 1], 1e-12, $"row {i}");
    }

    [TestMethod]
    public void DecisionTree_WeightTwo_IsTheSameAsADuplicatedRow()
    {
        var (X, y) = Blobs(30, 30, 1.0, seed: 5);
        var (Xtest, _) = Blobs(30, 30, 1.0, seed: 6);
        var parameters = new Dictionary<string, object> { ["MaxDepth"] = 3 };

        // Duplicate every third row, or give it weight 2.
        var rows = Enumerable.Range(0, y.Length).SelectMany(i => i % 3 == 0 ? new[] { i, i } : new[] { i }).ToArray();
        var duplicated = new DecisionTree();
        duplicated.SetHyperParameters(parameters);
        duplicated.Fit(X.SubMatrix(rows), new VectorN(rows.Select(i => y[i]).ToArray()));

        var weighted = new DecisionTree();
        weighted.SetHyperParameters(parameters);
        weighted.Fit(X, y, new VectorN(Enumerable.Range(0, y.Length).Select(i => i % 3 == 0 ? 2.0 : 1.0).ToArray()));

        var a = duplicated.PredictProbabilities(Xtest);
        var b = weighted.PredictProbabilities(Xtest);
        for (int i = 0; i < Xtest.rowLength; i++)
            Assert.AreEqual(a.values[i, 1], b.values[i, 1], 1e-9, $"row {i}");
    }

    [TestMethod]
    public void DecisionTree_MinSamplesLeaf_GivesUncertainLeaves()
    {
        var (X, y) = Blobs(50, 50, 1.0, seed: 7);

        var deep = new DecisionTree();
        deep.SetHyperParameters(new() { ["MaxDepth"] = 30 });
        deep.Fit(X, y);
        var pure = deep.PredictProbabilities(X);
        for (int i = 0; i < X.rowLength; i++)
            Assert.IsTrue(pure.values[i, 1] == 0 || pure.values[i, 1] == 1, "a fully grown tree has pure leaves");

        var smooth = new DecisionTree();
        smooth.SetHyperParameters(new() { ["MaxDepth"] = 30, ["MinSamplesLeaf"] = 10 });
        smooth.Fit(X, y);
        var p = smooth.PredictProbabilities(X);
        Assert.IsTrue(Enumerable.Range(0, X.rowLength).Any(i => p.values[i, 1] > 0 && p.values[i, 1] < 1));
    }

    [TestMethod]
    public void DecisionTree_InvalidWeights_Throw()
    {
        var (X, y) = Blobs(5, 5, 2.0, seed: 8);
        var tree = new DecisionTree();
        Assert.ThrowsException<ArgumentException>(() => tree.Fit(X, y, new VectorN(new double[3])));
        Assert.ThrowsException<ArgumentException>(() => tree.Fit(X, y, new VectorN(Enumerable.Repeat(-1.0, 10).ToArray())));
        Assert.ThrowsException<ArgumentException>(() => tree.Fit(X, y, new VectorN(new double[10])));
    }

    #endregion

    #region RandomForest

    [TestMethod]
    public void RandomForest_SeparatesBlobsWithGoodProbabilities()
    {
        var (X, y) = Blobs(150, 150, 3.0, seed: 10);
        var (Xtest, ytest) = Blobs(100, 100, 3.0, seed: 11);
        var forest = new RandomForest();
        forest.SetHyperParameters(new() { ["NumTrees"] = 40, ["MaxDepth"] = 6, ["MaxFeatures"] = 1 });
        forest.Fit(X, y);

        var p = forest.PredictProbabilities(Xtest);
        AssertRowsAreProbabilities(p);
        AssertArgmaxMatchesPredict(forest, Xtest);

        double accuracy = Enumerable.Range(0, ytest.Length).Count(i => forest.Predict(Xtest)[i] == ytest[i]) / (double)ytest.Length;
        double auc = ClassificationMetrics.RocAuc(ytest, ClassificationMetrics.ProbabilityOf(p));
        Assert.IsTrue(accuracy > 0.9, $"accuracy {accuracy:F3}");
        Assert.IsTrue(auc > 0.95, $"AUC {auc:F3}");
    }

    [TestMethod]
    public void RandomForest_IsRepeatableForASeed()
    {
        var (X, y) = Blobs(60, 60, 1.0, seed: 12);
        Matrix Train(int seed)
        {
            var forest = new RandomForest();
            forest.SetHyperParameters(new() { ["NumTrees"] = 20, ["MaxFeatures"] = 1, ["Seed"] = seed });
            forest.Fit(X, y);
            return forest.PredictProbabilities(X);
        }

        var a = Train(1);
        var b = Train(1);
        var c = Train(2);
        bool differs = false;
        for (int i = 0; i < X.rowLength; i++)
        {
            Assert.AreEqual(a.values[i, 1], b.values[i, 1], 0.0, $"row {i}");
            differs |= a.values[i, 1] != c.values[i, 1];
        }
        Assert.IsTrue(differs, "another seed should give another forest");
    }

    [TestMethod]
    public void RandomForest_BalancedWeights_FindMoreOfARareClass()
    {
        var (X, y) = Blobs(270, 30, 1.2, seed: 13);
        var (Xtest, ytest) = Blobs(270, 30, 1.2, seed: 14);
        var parameters = new Dictionary<string, object> { ["NumTrees"] = 30, ["MaxDepth"] = 4 };

        var plain = new RandomForest();
        plain.SetHyperParameters(parameters);
        plain.Fit(X, y);
        var balanced = new RandomForest();
        balanced.SetHyperParameters(parameters);
        balanced.Fit(X, y, ClassWeights.Balanced(y));

        double plainRecall = Recall(ytest, plain.Predict(Xtest));
        double balancedRecall = Recall(ytest, balanced.Predict(Xtest));
        Assert.IsTrue(balancedRecall > plainRecall, $"balanced {balancedRecall:F2}, plain {plainRecall:F2}");
    }

    [TestMethod]
    public void RandomForest_Clone_KeepsTheNewHyperparameters()
    {
        var forest = new RandomForest();
        forest.SetHyperParameters(new() { ["MinSamplesLeaf"] = 4, ["MaxFeatures"] = 2, ["Seed"] = 9 });

        var clone = (RandomForest)forest.Clone();

        Assert.AreEqual(4, clone.MinSamplesLeaf);
        Assert.AreEqual(2, clone.MaxFeatures);
        Assert.AreEqual(9, clone.Seed);
    }

    #endregion

    #region Logistic, MLP, CNN1D

    [TestMethod]
    public void Logistic_Probabilities_RiseWithTheFeatureAndAgreeWithPredict()
    {
        var X = new Matrix(Enumerable.Range(0, 20).Select(i => (double)i).ToArray().ToColumn());
        var y = new VectorN(Enumerable.Range(0, 20).Select(i => i >= 10 ? 1.0 : 0.0).ToArray());
        var model = new Logistic();
        model.SetHyperParameters(new() { ["LearningRate"] = 0.1, ["MaxIterations"] = 3000 });
        model.Fit(X, y);

        var p = model.PredictProbabilities(X);
        AssertRowsAreProbabilities(p);
        AssertArgmaxMatchesPredict(model, X);
        for (int i = 1; i < X.rowLength; i++)
            Assert.IsTrue(p.values[i, 1] >= p.values[i - 1, 1], $"row {i}");
    }

    [TestMethod]
    public void Logistic_EqualWeights_GiveTheSameModelAsNoWeights()
    {
        var (X, y) = Blobs(30, 30, 1.0, seed: 15);
        var plain = new Logistic();
        plain.Fit(X, y);
        var weighted = new Logistic();
        weighted.Fit(X, y, new VectorN(Enumerable.Repeat(1.0, y.Length).ToArray()));

        var a = plain.PredictProbabilities(X);
        var b = weighted.PredictProbabilities(X);
        for (int i = 0; i < X.rowLength; i++)
            Assert.AreEqual(a.values[i, 1], b.values[i, 1], 1e-12);
    }

    [TestMethod]
    public void MLPClassifier_Probabilities_SumToOneAndAgreeWithPredict()
    {
        var (X, y) = Blobs(40, 40, 3.0, seed: 16);
        var model = new MLPClassifier { HiddenLayers = new[] { 8 }, Epochs = 40, ValidationSplit = 0.0 };
        model.Fit(X, y, ClassWeights.Balanced(y));

        var p = model.PredictProbabilities(X);
        Assert.AreEqual(2, p.columnLength);
        AssertRowsAreProbabilities(p);
        AssertArgmaxMatchesPredict(model, X);
    }

    [TestMethod]
    public void CNN1DClassifier_Probabilities_SumToOneAndAgreeWithPredict()
    {
        // A dip in the middle of the sequence for class 1.
        var rng = new Random(17);
        int n = 24, steps = 10;
        var data = new double[n, steps];
        var labels = new double[n];
        for (int i = 0; i < n; i++)
        {
            labels[i] = i % 2;
            for (int t = 0; t < steps; t++)
                data[i, t] = 0.1 * Gaussian(rng) - (labels[i] == 1 && t >= 4 && t <= 5 ? 1.0 : 0.0);
        }
        var X = new Matrix(data);
        var model = new CNN1DClassifier();
        model.SetHyperParameters(new()
        {
            ["TimeSteps"] = steps, ["Features"] = 1, ["Filters"] = 3, ["KernelSize"] = 3,
            ["HiddenUnits"] = 4, ["Epochs"] = 10, ["ValidationSplit"] = 0.0,
        });
        model.Fit(X, new VectorN(labels));

        var p = model.PredictProbabilities(X);
        AssertRowsAreProbabilities(p);
        AssertArgmaxMatchesPredict(model, X);
    }

    #endregion

    #region Pipeline and class weights

    [TestMethod]
    public void Pipeline_PredictProbabilities_AppliesTheScaler()
    {
        var (X, y) = Blobs(50, 50, 3.0, seed: 18);
        var pipeline = new Pipeline(new RandomForest(), new() { ["NumTrees"] = 10 }, new StandardScaler());
        pipeline.Fit(X, y, ClassWeights.Balanced(y));

        var p = pipeline.PredictProbabilities(X);
        Assert.AreEqual(X.rowLength, p.rowLength);
        AssertRowsAreProbabilities(p);
    }

    [TestMethod]
    public void Pipeline_ModelWithoutWeightsOrProbabilities_Throws()
    {
        var (X, y) = Blobs(10, 10, 3.0, seed: 19);
        var pipeline = new Pipeline(new KNearestNeighbors(), new());

        Assert.ThrowsException<InvalidOperationException>(() => pipeline.Fit(X, y, ClassWeights.Balanced(y)));
        pipeline.Fit(X, y);
        Assert.ThrowsException<InvalidOperationException>(() => pipeline.PredictProbabilities(X));
    }

    [TestMethod]
    public void ClassWeights_Balanced_GiveEveryClassTheSameTotal()
    {
        var y = new VectorN(new double[] { 0, 0, 0, 0, 0, 0, 1, 1, 2, 2 });

        var w = ClassWeights.Balanced(y);

        Assert.AreEqual(10.0 / 3, Enumerable.Range(0, 6).Sum(i => w[i]), 1e-12);
        Assert.AreEqual(10.0 / 3, w[6] + w[7], 1e-12);
        Assert.AreEqual(10.0 / 3, w[8] + w[9], 1e-12);
        Assert.AreEqual(1.0, w.Values.Average(), 1e-12);
    }

    #endregion
}

internal static class TestArrayExtensions
{
    public static double[,] ToColumn(this double[] values)
    {
        var column = new double[values.Length, 1];
        for (int i = 0; i < values.Length; i++) column[i, 0] = values[i];
        return column;
    }
}
