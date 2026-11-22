namespace CSharpNumerics.ML.Evaluation;

/// <summary>
/// Receiver operating characteristic: true- and false-positive rates as the decision threshold
/// is lowered, from (0, 0) to (1, 1). See <see cref="ClassificationMetrics.Roc"/>.
/// </summary>
public sealed class RocCurve
{
    internal RocCurve(double[] falsePositiveRate, double[] truePositiveRate, double[] thresholds, double auc)
    {
        FalsePositiveRate = falsePositiveRate;
        TruePositiveRate = truePositiveRate;
        Thresholds = thresholds;
        Auc = auc;
    }

    /// <summary>Share of negatives scored at or above each threshold.</summary>
    public double[] FalsePositiveRate { get; }

    /// <summary>Share of positives scored at or above each threshold (recall).</summary>
    public double[] TruePositiveRate { get; }

    /// <summary>
    /// Decision thresholds, decreasing; the first is +∞ (nothing predicted positive), giving the
    /// point (0, 0).
    /// </summary>
    public double[] Thresholds { get; }

    /// <summary>Area under the curve: the probability that a random positive scores above a random negative (ties count half).</summary>
    public double Auc { get; }
}
