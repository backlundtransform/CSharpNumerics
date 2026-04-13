using System;
using System.Linq;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Result of a MAD computation.
/// </summary>
public class MadResult
{
    /// <summary>The median of the input data.</summary>
    public double Median { get; }

    /// <summary>The Median Absolute Deviation: median(|xᵢ − median(x)|).</summary>
    public double Mad { get; }

    /// <summary>
    /// Consistency-adjusted MAD that estimates σ for normal data:
    /// σ̂ = k × MAD, where k = 1.4826 by default.
    /// </summary>
    public double ScaledMad { get; }

    public MadResult(double median, double mad, double scaledMad)
    {
        Median = median;
        Mad = mad;
        ScaledMad = scaledMad;
    }
}

/// <summary>
/// Median Absolute Deviation (MAD): a robust measure of statistical dispersion.
/// MAD = median(|xᵢ − median(x)|).
/// The scaled MAD (with consistency constant k = 1.4826) is an estimator of σ
/// for normally distributed data.
/// </summary>
public static class MedianAbsoluteDeviation
{
    /// <summary>Default consistency constant for normal distributions.</summary>
    public const double DefaultScale = 1.4826;

    /// <summary>
    /// Computes the Median Absolute Deviation and a consistency-scaled estimate.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="scale">Consistency constant (default 1.4826 for normal data).</param>
    /// <returns>A <see cref="MadResult"/> with median, MAD, and scaled MAD.</returns>
    public static MadResult Compute(double[] data, double scale = DefaultScale)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) throw new ArgumentException("Data must not be empty.", nameof(data));
        if (scale <= 0) throw new ArgumentOutOfRangeException(nameof(scale), "Must be positive.");

        double median = ComputeMedian(data);

        double[] absDeviations = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
            absDeviations[i] = Math.Abs(data[i] - median);

        double mad = ComputeMedian(absDeviations);
        double scaledMad = scale * mad;

        return new MadResult(median, mad, scaledMad);
    }

    internal static double ComputeMedian(double[] values)
    {
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        if (n % 2 == 1)
            return sorted[n / 2];
        return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }
}
