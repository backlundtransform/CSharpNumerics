namespace CSharpNumerics.ML.Evaluation;

/// <summary>
/// Precision and recall as the decision threshold is lowered. See
/// <see cref="ClassificationMetrics.PrecisionRecall"/>.
/// </summary>
public sealed class PrecisionRecallCurve
{
    internal PrecisionRecallCurve(double[] precision, double[] recall, double[] thresholds, double averagePrecision)
    {
        Precision = precision;
        Recall = recall;
        Thresholds = thresholds;
        AveragePrecision = averagePrecision;
    }

    /// <summary>Share of rows scored at or above each threshold that are positive.</summary>
    public double[] Precision { get; }

    /// <summary>Share of positives scored at or above each threshold; increasing.</summary>
    public double[] Recall { get; }

    /// <summary>Decision thresholds, decreasing: every distinct score.</summary>
    public double[] Thresholds { get; }

    /// <summary>
    /// Average precision, Σ (Rₖ − Rₖ₋₁) · Pₖ: the area under the curve without interpolation.
    /// Equals the share of positives for a model that ranks at random.
    /// </summary>
    public double AveragePrecision { get; }
}
