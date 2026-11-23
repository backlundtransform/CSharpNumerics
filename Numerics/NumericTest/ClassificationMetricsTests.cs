using CSharpNumerics.ML.Evaluation;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Linq;

namespace NumericTest;

/// <summary>Expected values are worked by hand; several are the examples in scikit-learn's documentation.</summary>
[TestClass]
public class ClassificationMetricsTests
{
    private static VectorN V(params double[] values) => new VectorN(values);

    private static readonly VectorN Labels = V(0, 0, 1, 1);
    private static readonly VectorN Scores = V(0.1, 0.4, 0.35, 0.8);

    [TestMethod]
    public void Roc_FollowsTheThresholdsDown()
    {
        var roc = ClassificationMetrics.Roc(Labels, Scores);

        CollectionAssert.AreEqual(new[] { double.PositiveInfinity, 0.8, 0.4, 0.35, 0.1 }, roc.Thresholds);
        CollectionAssert.AreEqual(new[] { 0.0, 0.0, 0.5, 0.5, 1.0 }, roc.FalsePositiveRate);
        CollectionAssert.AreEqual(new[] { 0.0, 0.5, 0.5, 1.0, 1.0 }, roc.TruePositiveRate);
        Assert.AreEqual(0.75, roc.Auc, 1e-12);
    }

    [TestMethod]
    public void RocAuc_PerfectReversedAndTied()
    {
        Assert.AreEqual(1.0, ClassificationMetrics.RocAuc(V(0, 0, 1, 1), V(0.1, 0.2, 0.8, 0.9)), 1e-12);
        Assert.AreEqual(0.0, ClassificationMetrics.RocAuc(V(0, 0, 1, 1), V(0.9, 0.8, 0.2, 0.1)), 1e-12);
        Assert.AreEqual(0.5, ClassificationMetrics.RocAuc(V(0, 1, 0, 1), V(0.5, 0.5, 0.5, 0.5)), 1e-12);
    }

    [TestMethod]
    public void RocAuc_WeightTwo_IsTheSameAsADuplicatedRow()
    {
        var y = V(0, 1, 0, 1, 1, 0);
        var s = V(0.2, 0.7, 0.6, 0.4, 0.9, 0.1);

        double weighted = ClassificationMetrics.RocAuc(y, s, V(1, 1, 2, 1, 1, 1));
        double duplicated = ClassificationMetrics.RocAuc(V(0, 1, 0, 0, 1, 1, 0), V(0.2, 0.7, 0.6, 0.6, 0.4, 0.9, 0.1));

        Assert.AreEqual(duplicated, weighted, 1e-12);
    }

    [TestMethod]
    public void AveragePrecision_MatchesTheWorkedExample()
    {
        Assert.AreEqual(0.5 * 1.0 + 0.5 * (2.0 / 3.0), ClassificationMetrics.AveragePrecision(Labels, Scores), 1e-12);
    }

    [TestMethod]
    public void PrecisionRecall_GivesPrecisionAndRecallPerThreshold()
    {
        var curve = ClassificationMetrics.PrecisionRecall(Labels, Scores);

        CollectionAssert.AreEqual(new[] { 0.8, 0.4, 0.35, 0.1 }, curve.Thresholds);
        CollectionAssert.AreEqual(new[] { 0.5, 0.5, 1.0, 1.0 }, curve.Recall);
        Assert.AreEqual(1.0, curve.Precision[0], 1e-12);
        Assert.AreEqual(0.5, curve.Precision[1], 1e-12);
        Assert.AreEqual(2.0 / 3.0, curve.Precision[2], 1e-12);
        Assert.AreEqual(0.5, curve.Precision[3], 1e-12);
    }

    [TestMethod]
    public void PrecisionAtRecall_PicksTheBestThresholdThatReachesTheRecall()
    {
        var (allPrecision, allThreshold) = ClassificationMetrics.PrecisionAtRecall(Labels, Scores, minRecall: 1.0);
        Assert.AreEqual(2.0 / 3.0, allPrecision, 1e-12);
        Assert.AreEqual(0.35, allThreshold);

        var (halfPrecision, halfThreshold) = ClassificationMetrics.PrecisionAtRecall(Labels, Scores, minRecall: 0.5);
        Assert.AreEqual(1.0, halfPrecision, 1e-12);
        Assert.AreEqual(0.8, halfThreshold);
    }

    [TestMethod]
    public void BrierScore_MatchesTheWorkedExample()
    {
        Assert.AreEqual(0.0375, ClassificationMetrics.BrierScore(V(0, 1, 1, 0), V(0.1, 0.9, 0.8, 0.3)), 1e-12);
    }

    [TestMethod]
    public void LogLoss_MatchesTheWorkedExample()
    {
        double expected = -(Math.Log(0.9) + Math.Log(0.9) + Math.Log(0.8) + Math.Log(0.65)) / 4;
        Assert.AreEqual(expected, ClassificationMetrics.LogLoss(V(1, 0, 0, 1), V(0.9, 0.1, 0.2, 0.65)), 1e-12);
        Assert.AreEqual(0.21616, expected, 1e-5);
    }

    [TestMethod]
    public void LogLoss_ConfidentMistakeIsLargeButFinite()
    {
        double loss = ClassificationMetrics.LogLoss(V(1), V(0.0));
        Assert.IsTrue(loss > 30 && !double.IsInfinity(loss), $"loss {loss}");
    }

    [TestMethod]
    public void Calibration_BinsPredictionsAndLeavesOutEmptyBins()
    {
        var curve = ClassificationMetrics.Calibration(V(0, 1, 1, 1), V(0.1, 0.15, 0.8, 0.85), bins: 4);

        CollectionAssert.AreEqual(new[] { 0.0, 0.75 }, curve.BinLower);
        Assert.AreEqual(0.125, curve.MeanPredicted[0], 1e-12);
        Assert.AreEqual(0.5, curve.FractionPositive[0], 1e-12);
        Assert.AreEqual(0.825, curve.MeanPredicted[1], 1e-12);
        Assert.AreEqual(1.0, curve.FractionPositive[1], 1e-12);
        CollectionAssert.AreEqual(new[] { 2, 2 }, curve.Counts);
    }

    [TestMethod]
    public void ProbabilityOf_ReturnsAColumn()
    {
        var p = new Matrix(new double[,] { { 0.3, 0.7 }, { 0.9, 0.1 } });
        CollectionAssert.AreEqual(new[] { 0.7, 0.1 }, ClassificationMetrics.ProbabilityOf(p).Values);
        CollectionAssert.AreEqual(new[] { 0.3, 0.9 }, ClassificationMetrics.ProbabilityOf(p, 0).Values);
    }

    [TestMethod]
    public void InvalidInput_Throws()
    {
        Assert.ThrowsException<ArgumentException>(() => ClassificationMetrics.RocAuc(V(0, 2), V(0.1, 0.2)), "label 2");
        Assert.ThrowsException<ArgumentException>(() => ClassificationMetrics.RocAuc(V(0, 0), V(0.1, 0.2)), "no positive");
        Assert.ThrowsException<ArgumentException>(() => ClassificationMetrics.RocAuc(V(1, 1), V(0.1, 0.2)), "no negative");
        Assert.ThrowsException<ArgumentException>(() => ClassificationMetrics.RocAuc(V(0, 1), V(0.1)), "length");
        Assert.ThrowsException<ArgumentException>(() => ClassificationMetrics.BrierScore(V(0, 1), V(0.1, 1.2)), "not a probability");
        Assert.ThrowsException<ArgumentException>(() => ClassificationMetrics.RocAuc(V(0, 1), V(double.NaN, 0.2)), "NaN");
    }
}
