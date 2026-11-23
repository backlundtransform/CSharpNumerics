using CSharpNumerics.ML.Models;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Evaluation;

/// <summary>
/// Metrics for binary classifiers that give a score or probability, such as
/// <see cref="Models.Interfaces.IProbabilisticClassifier"/>.
/// </summary>
/// <remarks>
/// Labels are 0 (negative) and 1 (positive). Scores can be any finite numbers where higher means
/// "more likely positive"; <see cref="BrierScore"/>, <see cref="LogLoss"/> and
/// <see cref="Calibration"/> need probabilities in [0, 1]. Every metric takes optional sample
/// weights: a weight of 2 counts a row twice.
/// </remarks>
public static class ClassificationMetrics
{
    /// <summary>Column <paramref name="classIndex"/> of a probability matrix, e.g. P(class 1) for a binary classifier.</summary>
    public static VectorN ProbabilityOf(Matrix probabilities, int classIndex = 1)
    {
        if (classIndex < 0 || classIndex >= probabilities.columnLength)
            throw new ArgumentOutOfRangeException(nameof(classIndex));
        var result = new double[probabilities.rowLength];
        for (int i = 0; i < result.Length; i++)
            result[i] = probabilities.values[i, classIndex];
        return new VectorN(result);
    }

    /// <summary>ROC curve and its area. Needs at least one positive and one negative.</summary>
    public static RocCurve Roc(VectorN yTrue, VectorN scores, VectorN? sampleWeights = null)
    {
        var sweep = Sweep(yTrue, scores, sampleWeights);
        if (sweep.Negatives <= 0)
            throw new ArgumentException("ROC needs at least one negative (label 0).", nameof(yTrue));

        int k = sweep.Thresholds.Length;
        var fpr = new double[k + 1];
        var tpr = new double[k + 1];
        var thresholds = new double[k + 1];
        thresholds[0] = double.PositiveInfinity;
        double auc = 0;
        for (int i = 0; i < k; i++)
        {
            fpr[i + 1] = sweep.FalsePositives[i] / sweep.Negatives;
            tpr[i + 1] = sweep.TruePositives[i] / sweep.Positives;
            thresholds[i + 1] = sweep.Thresholds[i];
            // Trapezoid: a group of tied scores moves diagonally, which counts ties as half.
            auc += (fpr[i + 1] - fpr[i]) * (tpr[i + 1] + tpr[i]) / 2.0;
        }

        return new RocCurve(fpr, tpr, thresholds, auc);
    }

    /// <summary>Area under the ROC curve. See <see cref="Roc"/>.</summary>
    public static double RocAuc(VectorN yTrue, VectorN scores, VectorN? sampleWeights = null) =>
        Roc(yTrue, scores, sampleWeights).Auc;

    /// <summary>Precision–recall curve and average precision. Needs at least one positive.</summary>
    public static PrecisionRecallCurve PrecisionRecall(VectorN yTrue, VectorN scores, VectorN? sampleWeights = null)
    {
        var sweep = Sweep(yTrue, scores, sampleWeights);
        int k = sweep.Thresholds.Length;
        var precision = new double[k];
        var recall = new double[k];
        double averagePrecision = 0, previousRecall = 0;
        for (int i = 0; i < k; i++)
        {
            double predicted = sweep.TruePositives[i] + sweep.FalsePositives[i];
            precision[i] = predicted > 0 ? sweep.TruePositives[i] / predicted : 1.0;
            recall[i] = sweep.TruePositives[i] / sweep.Positives;
            averagePrecision += (recall[i] - previousRecall) * precision[i];
            previousRecall = recall[i];
        }

        return new PrecisionRecallCurve(precision, recall, (double[])sweep.Thresholds.Clone(), averagePrecision);
    }

    /// <summary>Average precision (area under the precision–recall curve). See <see cref="PrecisionRecall"/>.</summary>
    public static double AveragePrecision(VectorN yTrue, VectorN scores, VectorN? sampleWeights = null) =>
        PrecisionRecall(yTrue, scores, sampleWeights).AveragePrecision;

    /// <summary>
    /// The highest precision among thresholds that reach at least <paramref name="minRecall"/>, and
    /// the highest such threshold: "if we must find 90 % of the planets, how pure is the list?"
    /// </summary>
    public static (double Precision, double Threshold) PrecisionAtRecall(VectorN yTrue, VectorN scores,
        double minRecall, VectorN? sampleWeights = null)
    {
        if (!(minRecall >= 0 && minRecall <= 1))
            throw new ArgumentOutOfRangeException(nameof(minRecall), "Recall must be between 0 and 1.");

        var curve = PrecisionRecall(yTrue, scores, sampleWeights);
        double bestPrecision = -1, bestThreshold = double.NaN;
        for (int i = 0; i < curve.Recall.Length; i++)
        {
            if (curve.Recall[i] < minRecall - 1e-12) continue;
            if (curve.Precision[i] > bestPrecision)
            {
                bestPrecision = curve.Precision[i];
                bestThreshold = curve.Thresholds[i];
            }
        }

        return (bestPrecision, bestThreshold);
    }

    /// <summary>Mean squared difference between probability and label; 0 is perfect, 0.25 is a constant 0.5.</summary>
    public static double BrierScore(VectorN yTrue, VectorN probabilities, VectorN? sampleWeights = null)
    {
        var (labels, p, w) = CheckProbabilities(yTrue, probabilities, sampleWeights);
        double sum = 0, total = 0;
        for (int i = 0; i < p.Length; i++)
        {
            double d = p[i] - labels[i];
            sum += w[i] * d * d;
            total += w[i];
        }
        return sum / total;
    }

    /// <summary>
    /// Mean negative log-likelihood (cross-entropy) of the labels. Probabilities are clipped to
    /// [<paramref name="epsilon"/>, 1 − <paramref name="epsilon"/>] so a confident mistake is
    /// large but finite.
    /// </summary>
    public static double LogLoss(VectorN yTrue, VectorN probabilities, VectorN? sampleWeights = null, double epsilon = 1e-15)
    {
        var (labels, p, w) = CheckProbabilities(yTrue, probabilities, sampleWeights);
        double sum = 0, total = 0;
        for (int i = 0; i < p.Length; i++)
        {
            double q = Math.Min(Math.Max(p[i], epsilon), 1 - epsilon);
            sum -= w[i] * (labels[i] == 1 ? Math.Log(q) : Math.Log(1 - q));
            total += w[i];
        }
        return sum / total;
    }

    /// <summary>
    /// Calibration (reliability) curve over <paramref name="bins"/> equal-width probability bins.
    /// Empty bins are left out.
    /// </summary>
    public static CalibrationCurve Calibration(VectorN yTrue, VectorN probabilities, int bins = 10, VectorN? sampleWeights = null)
    {
        if (bins < 1) throw new ArgumentOutOfRangeException(nameof(bins));
        var (labels, p, w) = CheckProbabilities(yTrue, probabilities, sampleWeights);

        var predicted = new double[bins];
        var positive = new double[bins];
        var weight = new double[bins];
        var counts = new int[bins];
        for (int i = 0; i < p.Length; i++)
        {
            int b = Math.Min((int)(p[i] * bins), bins - 1);
            predicted[b] += w[i] * p[i];
            positive[b] += w[i] * labels[i];
            weight[b] += w[i];
            counts[b]++;
        }

        var lower = new List<double>();
        var mean = new List<double>();
        var fraction = new List<double>();
        var count = new List<int>();
        for (int b = 0; b < bins; b++)
        {
            if (counts[b] == 0 || weight[b] <= 0) continue;
            lower.Add((double)b / bins);
            mean.Add(predicted[b] / weight[b]);
            fraction.Add(positive[b] / weight[b]);
            count.Add(counts[b]);
        }

        return new CalibrationCurve(lower.ToArray(), mean.ToArray(), fraction.ToArray(), count.ToArray());
    }

    private sealed class SweepResult
    {
        public double[] Thresholds;
        public double[] TruePositives;
        public double[] FalsePositives;
        public double Positives;
        public double Negatives;
    }

    /// <summary>Cumulative weighted true and false positives at each distinct score, from the highest down.</summary>
    private static SweepResult Sweep(VectorN yTrue, VectorN scores, VectorN? sampleWeights)
    {
        int[] labels = Labels(yTrue);
        if (scores.Values == null || scores.Length != labels.Length)
            throw new ArgumentException("Scores must have one value per label.", nameof(scores));
        double[] w = ClassWeights.Validate(sampleWeights, labels.Length);

        int n = labels.Length;
        var keys = new double[n];
        var order = new int[n];
        for (int i = 0; i < n; i++)
        {
            double s = scores[i];
            if (double.IsNaN(s) || double.IsInfinity(s))
                throw new ArgumentException("Scores must be finite.", nameof(scores));
            keys[i] = -s; // ascending sort of −score = descending score
            order[i] = i;
        }
        Array.Sort(keys, order);

        var thresholds = new List<double>();
        var tps = new List<double>();
        var fps = new List<double>();
        double tp = 0, fp = 0;
        for (int k = 0; k < n; k++)
        {
            int i = order[k];
            if (labels[i] == 1) tp += w[i]; else fp += w[i];
            if (k + 1 < n && keys[k + 1] == keys[k]) continue; // finish the group of tied scores
            thresholds.Add(-keys[k]);
            tps.Add(tp);
            fps.Add(fp);
        }

        if (tp <= 0)
            throw new ArgumentException("At least one positive (label 1) with a positive weight is required.", nameof(yTrue));

        return new SweepResult
        {
            Thresholds = thresholds.ToArray(),
            TruePositives = tps.ToArray(),
            FalsePositives = fps.ToArray(),
            Positives = tp,
            Negatives = fp,
        };
    }

    private static (int[] Labels, double[] Probabilities, double[] Weights) CheckProbabilities(
        VectorN yTrue, VectorN probabilities, VectorN? sampleWeights)
    {
        int[] labels = Labels(yTrue);
        if (probabilities.Values == null || probabilities.Length != labels.Length)
            throw new ArgumentException("Probabilities must have one value per label.", nameof(probabilities));
        foreach (double p in probabilities.Values)
            if (!(p >= 0 && p <= 1))
                throw new ArgumentException("Probabilities must be between 0 and 1.", nameof(probabilities));
        return (labels, probabilities.Values, ClassWeights.Validate(sampleWeights, labels.Length));
    }

    private static int[] Labels(VectorN yTrue)
    {
        if (yTrue.Values == null || yTrue.Length == 0)
            throw new ArgumentException("No labels.", nameof(yTrue));
        var labels = new int[yTrue.Length];
        for (int i = 0; i < labels.Length; i++)
        {
            if (yTrue[i] == 1) labels[i] = 1;
            else if (yTrue[i] != 0) throw new ArgumentException($"Labels must be 0 or 1; got {yTrue[i]}.", nameof(yTrue));
        }
        return labels;
    }
}
