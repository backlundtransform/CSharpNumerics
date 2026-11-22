using CSharpNumerics.ML;
using CSharpNumerics.ML.CrossValidators;
using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Experiment;
using CSharpNumerics.ML.Models;
using CSharpNumerics.ML.Models.Classification;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericTest;

[TestClass]
public class GroupedCrossValidationTests
{
    #region Data

    private static double Gaussian(Random rng) =>
        Math.Sqrt(-2.0 * Math.Log(1.0 - rng.NextDouble())) * Math.Cos(2.0 * Math.PI * rng.NextDouble());

    /// <summary>
    /// <paramref name="groupCount"/> groups of 1–5 rows. Every row of a group has the group's label
    /// (1 with probability <paramref name="positiveShare"/>), and its features are the group's
    /// centre (shifted by the label) plus a little noise.
    /// </summary>
    private static (Matrix X, VectorN y, long[] Groups) GroupedData(int groupCount, double positiveShare,
        double separation, int seed)
    {
        var rng = new Random(seed);
        var rows = new List<(double A, double B, double Label, long Group)>();
        for (int g = 0; g < groupCount; g++)
        {
            double label = rng.NextDouble() < positiveShare ? 1 : 0;
            double a = Gaussian(rng) + separation * label, b = Gaussian(rng) + separation * label;
            int size = 1 + rng.Next(5);
            for (int k = 0; k < size; k++)
                rows.Add((a + 0.05 * Gaussian(rng), b + 0.05 * Gaussian(rng), label, 1000L + g * 7));
        }

        var X = new double[rows.Count, 2];
        for (int i = 0; i < rows.Count; i++) { X[i, 0] = rows[i].A; X[i, 1] = rows[i].B; }
        return (new Matrix(X), new VectorN(rows.Select(r => r.Label).ToArray()), rows.Select(r => r.Group).ToArray());
    }

    private static List<Pipeline> Forest(int trees = 15) =>
        [new Pipeline(new RandomForest(), new() { ["NumTrees"] = trees, ["MaxDepth"] = 5 })];

    private static void AssertGroupsStayTogether(long[] groups, int[] folds)
    {
        foreach (var group in Enumerable.Range(0, groups.Length).GroupBy(i => groups[i]))
            Assert.AreEqual(1, group.Select(i => folds[i]).Distinct().Count(), $"group {group.Key} is split");
    }

    #endregion

    [TestMethod]
    public void GroupKFold_KeepsGroupsTogetherAndFoldsBalanced()
    {
        var (_, y, groups) = GroupedData(120, 0.3, 2.0, seed: 1);
        var cv = new GroupKFoldCrossValidator(Forest(), groups, folds: 5);

        int[] folds = cv.FoldOfEachRow(y.Length);

        AssertGroupsStayTogether(groups, folds);
        var sizes = Enumerable.Range(0, 5).Select(f => folds.Count(x => x == f)).ToArray();
        Assert.IsTrue(sizes.Max() - sizes.Min() <= 5, $"fold sizes {string.Join(", ", sizes)}"); // at most one group apart
    }

    [TestMethod]
    public void StratifiedGroupKFold_KeepsGroupsTogetherAndClassSharesEven()
    {
        var (_, y, groups) = GroupedData(300, 0.2, 2.0, seed: 2);
        var cv = new StratifiedGroupKFoldCrossValidator(Forest(), groups, folds: 5);

        int[] folds = cv.FoldOfEachRow(y);

        AssertGroupsStayTogether(groups, folds);
        double overall = y.Values.Average();
        for (int f = 0; f < 5; f++)
        {
            double share = Enumerable.Range(0, y.Length).Where(i => folds[i] == f).Average(i => y[i]);
            Assert.AreEqual(overall, share, 0.05, $"fold {f}");
        }
    }

    [TestMethod]
    public void GroupKFold_TooFewGroups_Throws()
    {
        var y = new VectorN(new double[] { 0, 1, 0, 1 });
        var cv = new GroupKFoldCrossValidator(Forest(), new long[] { 1, 1, 2, 2 }, folds: 3);
        Assert.ThrowsException<InvalidOperationException>(() => cv.FoldOfEachRow(y.Length));
    }

    [TestMethod]
    public void GroupKFold_DoesNotLeakGroupsIntoTheScore()
    {
        // Labels are random per group and the features only identify the group: there is nothing
        // to learn. Plain k-fold rewards memorising groups; group k-fold does not.
        var (ordered, labels, orderedGroups) = GroupedData(150, 0.5, separation: 0.0, seed: 3);
        // Shuffle the rows: k-fold splits in blocks, and group members would otherwise sit together.
        var rng = new Random(3);
        var order = Enumerable.Range(0, labels.Length).OrderBy(_ => rng.Next()).ToArray();
        var X = ordered.SubMatrix(order);
        var y = new VectorN(order.Select(i => labels[i]).ToArray());
        var groups = order.Select(i => orderedGroups[i]).ToArray();
        var tree = new List<Pipeline> { new Pipeline(new DecisionTree(), new() { ["MaxDepth"] = 30 }) };

        double leaky = new KFoldCrossValidator(tree, 5).Run(X, y).BestScore;
        double grouped = new GroupKFoldCrossValidator(tree, groups, folds: 5).Run(X, y).BestScore;

        Assert.IsTrue(leaky > 0.7, $"plain k-fold accuracy {leaky:F2}");
        Assert.IsTrue(grouped < 0.65, $"group k-fold accuracy {grouped:F2}");
    }

    [TestMethod]
    public void Run_ScoresByRocAucAndReturnsOutOfFoldProbabilities()
    {
        var (X, y, groups) = GroupedData(150, 0.4, 3.0, seed: 4);
        var cv = new StratifiedGroupKFoldCrossValidator(Forest(), groups, folds: 4, scoring: ClassificationScoring.RocAuc);

        var result = cv.Run(X, y);

        Assert.IsTrue(result.BestScore > 0.9, $"AUC {result.BestScore:F3}");
        Assert.IsTrue(result.PredictedProbabilities.HasValue);
        Assert.AreEqual(y.Length, result.PredictedProbabilities.Value.rowLength);
        Assert.AreEqual(2, result.PredictedProbabilities.Value.columnLength);
        CollectionAssert.AreEqual(Enumerable.Range(0, y.Length).ToArray(), result.RowIndices);
        CollectionAssert.AreEqual(y.Values, result.ActualValues.Values);
    }

    [TestMethod]
    public void Run_WithSampleWeights_TrainsWeightedModels()
    {
        var (X, y, groups) = GroupedData(120, 0.15, 1.0, seed: 5);
        var cv = new GroupKFoldCrossValidator(Forest(), groups, folds: 4, sampleWeights: ClassWeights.Balanced(y),
            scoring: ClassificationScoring.AveragePrecision);

        var result = cv.Run(X, y);

        Assert.IsTrue(result.BestScore > 0 && result.BestScore <= 1, $"AP {result.BestScore:F3}");
    }

    [TestMethod]
    public void Run_ProbabilityScoringWithoutProbabilities_Throws()
    {
        var (X, y, groups) = GroupedData(40, 0.5, 2.0, seed: 6);
        var knn = new List<Pipeline> { new Pipeline(new KNearestNeighbors(), new()) };

        var byAuc = new GroupKFoldCrossValidator(knn, groups, folds: 3, scoring: ClassificationScoring.RocAuc);
        Assert.ThrowsException<InvalidOperationException>(() => byAuc.Run(X, y));

        var weighted = new GroupKFoldCrossValidator(knn, groups, folds: 3, sampleWeights: ClassWeights.Balanced(y));
        Assert.ThrowsException<InvalidOperationException>(() => weighted.Run(X, y));
    }

    [TestMethod]
    public void SupervisedExperiment_RunsWithStratifiedGroupKFold()
    {
        var (X, y, groups) = GroupedData(100, 0.3, 3.0, seed: 7);

        var result = SupervisedExperiment
            .For(X, y)
            .WithGrid(new PipelineGrid()
                .AddModel<RandomForest>(g => g
                    .Add("NumTrees", 10)
                    .Add("MaxDepth", 3, 6)))
            .WithCrossValidator(CrossValidatorConfig.StratifiedGroupKFold(groups, folds: 4, scoring: ClassificationScoring.RocAuc))
            .Run();

        Assert.AreEqual(2, result.Rankings.Count);
        Assert.IsTrue(result.BestScore > 0.9, $"AUC {result.BestScore:F3}");
    }
}
