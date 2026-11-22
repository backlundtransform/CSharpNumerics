namespace CSharpNumerics.ML.Evaluation;

/// <summary>
/// Predicted probability against the observed share of positives, per probability bin. A
/// calibrated model has <see cref="MeanPredicted"/> ≈ <see cref="FractionPositive"/>. See
/// <see cref="ClassificationMetrics.Calibration"/>.
/// </summary>
public sealed class CalibrationCurve
{
    internal CalibrationCurve(double[] binLower, double[] meanPredicted, double[] fractionPositive, int[] counts)
    {
        BinLower = binLower;
        MeanPredicted = meanPredicted;
        FractionPositive = fractionPositive;
        Counts = counts;
    }

    /// <summary>Lower edge of each non-empty bin.</summary>
    public double[] BinLower { get; }

    /// <summary>Mean predicted probability in each bin.</summary>
    public double[] MeanPredicted { get; }

    /// <summary>Share of positives in each bin.</summary>
    public double[] FractionPositive { get; }

    /// <summary>Rows in each bin.</summary>
    public int[] Counts { get; }
}
