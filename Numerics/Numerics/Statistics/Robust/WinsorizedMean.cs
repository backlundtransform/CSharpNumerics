using System;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Winsorized mean: computes the arithmetic mean after replacing a specified
/// proportion of the smallest values with the next smallest retained value,
/// and the largest values with the next largest retained value.
/// Unlike trimming, winsorizing preserves the original sample size.
/// </summary>
public static class WinsorizedMean
{
    /// <summary>
    /// Computes the winsorized mean of the data.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="trimRatio">
    /// Proportion of data to winsorize from each end, in [0, 0.5).
    /// A value of 0.1 replaces the lowest 10 % and highest 10 % with boundary values.
    /// </param>
    /// <returns>The winsorized mean.</returns>
    public static double Compute(double[] data, double trimRatio = 0.1)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) throw new ArgumentException("Data must not be empty.", nameof(data));
        if (trimRatio < 0 || trimRatio >= 0.5)
            throw new ArgumentOutOfRangeException(nameof(trimRatio), "Must be in [0, 0.5).");

        double[] sorted = (double[])data.Clone();
        Array.Sort(sorted);

        int trimCount = (int)Math.Floor(sorted.Length * trimRatio);

        if (trimCount > 0)
        {
            double lowerBound = sorted[trimCount];
            double upperBound = sorted[sorted.Length - 1 - trimCount];

            for (int i = 0; i < trimCount; i++)
                sorted[i] = lowerBound;

            for (int i = sorted.Length - trimCount; i < sorted.Length; i++)
                sorted[i] = upperBound;
        }

        double sum = 0;
        for (int i = 0; i < sorted.Length; i++)
            sum += sorted[i];

        return sum / sorted.Length;
    }
}
