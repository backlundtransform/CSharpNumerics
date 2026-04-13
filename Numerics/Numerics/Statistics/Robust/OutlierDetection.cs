using System;
using System.Collections.Generic;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Result of outlier detection.
/// </summary>
public class OutlierResult
{
    /// <summary>Boolean mask: true = outlier, false = inlier.</summary>
    public bool[] OutlierMask { get; }

    /// <summary>Indices of detected outliers.</summary>
    public int[] OutlierIndices { get; }

    /// <summary>Number of detected outliers.</summary>
    public int OutlierCount { get; }

    /// <summary>Scores used for outlier determination (interpretation depends on method).</summary>
    public double[] Scores { get; }

    public OutlierResult(bool[] outlierMask, int[] outlierIndices, int outlierCount, double[] scores)
    {
        OutlierMask = outlierMask;
        OutlierIndices = outlierIndices;
        OutlierCount = outlierCount;
        Scores = scores;
    }
}

/// <summary>
/// Methods for detecting outliers in univariate data:
/// <list type="bullet">
///   <item><description><see cref="Iqr"/> — interquartile range fence method.</description></item>
///   <item><description><see cref="ZScore"/> — standard Z-score method.</description></item>
///   <item><description><see cref="ModifiedZScore"/> — MAD-based modified Z-score (Iglewicz &amp; Hoaglin).</description></item>
/// </list>
/// </summary>
public static class OutlierDetection
{
    /// <summary>
    /// Detects outliers using the IQR fence method.
    /// A point is an outlier if it falls below Q1 − k·IQR or above Q3 + k·IQR.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="k">Fence multiplier (default 1.5 for mild outliers, 3.0 for extreme).</param>
    /// <returns>An <see cref="OutlierResult"/> with outlier mask and deviation scores.</returns>
    public static OutlierResult Iqr(double[] data, double k = 1.5)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length < 4) throw new ArgumentException("Need at least 4 data points for IQR method.", nameof(data));
        if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k), "Must be positive.");

        double[] sorted = (double[])data.Clone();
        Array.Sort(sorted);

        double q1 = Percentile(sorted, 0.25);
        double q3 = Percentile(sorted, 0.75);
        double iqr = q3 - q1;

        double lowerFence = q1 - k * iqr;
        double upperFence = q3 + k * iqr;

        return BuildResult(data, i => data[i] < lowerFence || data[i] > upperFence,
            i => Math.Max(lowerFence - data[i], data[i] - upperFence));
    }

    /// <summary>
    /// Detects outliers using the standard Z-score method.
    /// A point is an outlier if |z| &gt; threshold.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="threshold">Z-score threshold (default 3.0).</param>
    /// <returns>An <see cref="OutlierResult"/> with outlier mask and absolute Z-scores.</returns>
    public static OutlierResult ZScore(double[] data, double threshold = 3.0)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length < 2) throw new ArgumentException("Need at least 2 data points.", nameof(data));
        if (threshold <= 0) throw new ArgumentOutOfRangeException(nameof(threshold), "Must be positive.");

        double mean = 0;
        for (int i = 0; i < data.Length; i++) mean += data[i];
        mean /= data.Length;

        double variance = 0;
        for (int i = 0; i < data.Length; i++)
        {
            double d = data[i] - mean;
            variance += d * d;
        }
        double std = Math.Sqrt(variance / data.Length);

        if (std < 1e-15)
            return BuildResult(data, _ => false, _ => 0.0);

        double[] zScores = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
            zScores[i] = Math.Abs((data[i] - mean) / std);

        return BuildResult(data, i => zScores[i] > threshold, i => zScores[i]);
    }

    /// <summary>
    /// Detects outliers using the modified Z-score based on MAD (Iglewicz &amp; Hoaglin, 1993).
    /// Modified Z = 0.6745 × (xᵢ − median) / MAD.
    /// A point is an outlier if |modified Z| &gt; threshold.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="threshold">Modified Z-score threshold (default 3.5).</param>
    /// <returns>An <see cref="OutlierResult"/> with outlier mask and absolute modified Z-scores.</returns>
    public static OutlierResult ModifiedZScore(double[] data, double threshold = 3.5)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length < 2) throw new ArgumentException("Need at least 2 data points.", nameof(data));
        if (threshold <= 0) throw new ArgumentOutOfRangeException(nameof(threshold), "Must be positive.");

        double median = MedianAbsoluteDeviation.ComputeMedian(data);

        double[] absDeviations = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
            absDeviations[i] = Math.Abs(data[i] - median);

        double mad = MedianAbsoluteDeviation.ComputeMedian(absDeviations);

        if (mad < 1e-15)
            return BuildResult(data, _ => false, _ => 0.0);

        double[] scores = new double[data.Length];
        for (int i = 0; i < data.Length; i++)
            scores[i] = Math.Abs(0.6745 * (data[i] - median) / mad);

        return BuildResult(data, i => scores[i] > threshold, i => scores[i]);
    }

    private static OutlierResult BuildResult(double[] data, Func<int, bool> isOutlier, Func<int, double> score)
    {
        bool[] mask = new bool[data.Length];
        double[] scores = new double[data.Length];
        var outlierIndices = new List<int>();

        for (int i = 0; i < data.Length; i++)
        {
            scores[i] = score(i);
            mask[i] = isOutlier(i);
            if (mask[i]) outlierIndices.Add(i);
        }

        return new OutlierResult(mask, outlierIndices.ToArray(), outlierIndices.Count, scores);
    }

    private static double Percentile(double[] sorted, double p)
    {
        double index = p * (sorted.Length - 1);
        int lower = (int)Math.Floor(index);
        int upper = (int)Math.Ceiling(index);
        if (lower == upper)
            return sorted[lower];
        double fraction = index - lower;
        return sorted[lower] * (1 - fraction) + sorted[upper] * fraction;
    }
}
