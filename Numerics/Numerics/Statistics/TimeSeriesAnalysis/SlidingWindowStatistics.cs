using System;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Computes running (sliding window) statistics: mean, standard deviation,
/// and median absolute deviation (MAD) over a configurable window size.
/// </summary>
public static class SlidingWindowStatistics
{
    /// <summary>
    /// Computes the running mean with a centered window of given size.
    /// Edge windows are computed with fewer points (minimum 1).
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="windowSize">Window size (must be odd for symmetric centering; even values are incremented by 1).</param>
    /// <returns>Array of running means, same length as input.</returns>
    public static double[] RunningMean(double[] data, int windowSize)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) return Array.Empty<double>();
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize), "Must be at least 1.");

        // Make window odd for symmetric centering
        if (windowSize % 2 == 0) windowSize++;
        int halfWin = windowSize / 2;

        double[] result = new double[data.Length];

        for (int i = 0; i < data.Length; i++)
        {
            int start = Math.Max(0, i - halfWin);
            int end = Math.Min(data.Length - 1, i + halfWin);
            double sum = 0;
            int count = 0;
            for (int j = start; j <= end; j++)
            {
                sum += data[j];
                count++;
            }
            result[i] = sum / count;
        }

        return result;
    }

    /// <summary>
    /// Computes the running standard deviation with a centered window.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="windowSize">Window size.</param>
    /// <returns>Array of running standard deviations, same length as input.</returns>
    public static double[] RunningStd(double[] data, int windowSize)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) return Array.Empty<double>();
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize), "Must be at least 1.");

        if (windowSize % 2 == 0) windowSize++;
        int halfWin = windowSize / 2;

        double[] result = new double[data.Length];

        for (int i = 0; i < data.Length; i++)
        {
            int start = Math.Max(0, i - halfWin);
            int end = Math.Min(data.Length - 1, i + halfWin);

            double sum = 0;
            int count = 0;
            for (int j = start; j <= end; j++)
            {
                sum += data[j];
                count++;
            }
            double mean = sum / count;

            double variance = 0;
            for (int j = start; j <= end; j++)
            {
                double d = data[j] - mean;
                variance += d * d;
            }

            result[i] = count > 1 ? Math.Sqrt(variance / count) : 0.0;
        }

        return result;
    }

    /// <summary>
    /// Computes the running median absolute deviation (MAD) with a centered window.
    /// MAD = median(|x_i - median(x)|) within the window.
    /// The MAD is a robust measure of spread, less sensitive to outliers than std.
    /// To convert to a std-equivalent, multiply by 1.4826.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="windowSize">Window size.</param>
    /// <returns>Array of running MAD values, same length as input.</returns>
    public static double[] RunningMAD(double[] data, int windowSize)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) return Array.Empty<double>();
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize), "Must be at least 1.");

        if (windowSize % 2 == 0) windowSize++;
        int halfWin = windowSize / 2;

        double[] result = new double[data.Length];

        for (int i = 0; i < data.Length; i++)
        {
            int start = Math.Max(0, i - halfWin);
            int end = Math.Min(data.Length - 1, i + halfWin);
            int count = end - start + 1;

            // Extract window values
            double[] window = new double[count];
            for (int j = 0; j < count; j++)
                window[j] = data[start + j];

            // Compute median
            double median = Median(window);

            // Compute absolute deviations from median
            double[] absDevs = new double[count];
            for (int j = 0; j < count; j++)
                absDevs[j] = Math.Abs(window[j] - median);

            // MAD = median of absolute deviations
            result[i] = Median(absDevs);
        }

        return result;
    }

    /// <summary>
    /// Computes the running median with a centered window.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="windowSize">Window size.</param>
    /// <returns>Array of running medians, same length as input.</returns>
    public static double[] RunningMedian(double[] data, int windowSize)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) return Array.Empty<double>();
        if (windowSize < 1) throw new ArgumentOutOfRangeException(nameof(windowSize), "Must be at least 1.");

        if (windowSize % 2 == 0) windowSize++;
        int halfWin = windowSize / 2;

        double[] result = new double[data.Length];

        for (int i = 0; i < data.Length; i++)
        {
            int start = Math.Max(0, i - halfWin);
            int end = Math.Min(data.Length - 1, i + halfWin);
            int count = end - start + 1;

            double[] window = new double[count];
            for (int j = 0; j < count; j++)
                window[j] = data[start + j];

            result[i] = Median(window);
        }

        return result;
    }

    /// <summary>
    /// Computes the median of an array (sorts in-place).
    /// </summary>
    private static double Median(double[] sorted)
    {
        Array.Sort(sorted);
        int n = sorted.Length;
        if (n == 0) return 0;
        return n % 2 == 1
            ? sorted[n / 2]
            : (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }
}
