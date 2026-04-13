using System;
using System.Linq;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Trimmed mean: computes the arithmetic mean after discarding a specified
/// proportion of the smallest and largest values.
/// Trimming ratio α ∈ [0, 0.5): α = 0 gives the ordinary mean,
/// α → 0.5 converges on the median.
/// </summary>
public static class TrimmedMean
{
    /// <summary>
    /// Computes the trimmed mean of the data.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="trimRatio">
    /// Proportion of data to trim from each end, in [0, 0.5).
    /// A value of 0.1 trims the lowest 10 % and highest 10 %.
    /// </param>
    /// <returns>The trimmed mean.</returns>
    public static double Compute(double[] data, double trimRatio = 0.1)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) throw new ArgumentException("Data must not be empty.", nameof(data));
        if (trimRatio < 0 || trimRatio >= 0.5)
            throw new ArgumentOutOfRangeException(nameof(trimRatio), "Must be in [0, 0.5).");

        double[] sorted = (double[])data.Clone();
        Array.Sort(sorted);

        int trimCount = (int)Math.Floor(sorted.Length * trimRatio);
        int start = trimCount;
        int end = sorted.Length - trimCount;

        if (end <= start)
            return sorted[sorted.Length / 2];

        double sum = 0;
        for (int i = start; i < end; i++)
            sum += sorted[i];

        return sum / (end - start);
    }
}
