using System;
using System.Collections.Generic;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Result of iterative sigma-clipping.
/// </summary>
public class ClipResult
{
    /// <summary>Boolean mask: true = retained, false = clipped.</summary>
    public bool[] Mask { get; }

    /// <summary>Mean of the retained data after final iteration.</summary>
    public double Mean { get; }

    /// <summary>Standard deviation of the retained data after final iteration.</summary>
    public double Std { get; }

    /// <summary>Number of retained data points.</summary>
    public int RetainedCount { get; }

    /// <summary>Number of clipped data points.</summary>
    public int ClippedCount { get; }

    /// <summary>Number of iterations performed.</summary>
    public int Iterations { get; }

    public ClipResult(bool[] mask, double mean, double std, int retainedCount, int clippedCount, int iterations)
    {
        Mask = mask;
        Mean = mean;
        Std = std;
        RetainedCount = retainedCount;
        ClippedCount = clippedCount;
        Iterations = iterations;
    }
}

/// <summary>
/// Iterative sigma-clipping: removes data points that fall outside
/// [mean − sigmaLow × σ, mean + sigmaHigh × σ] and recomputes statistics
/// until convergence or maxIter is reached.
/// </summary>
public static class SigmaClipping
{
    /// <summary>
    /// Performs iterative sigma-clipping on the input data.
    /// </summary>
    /// <param name="data">Input data array.</param>
    /// <param name="sigmaLow">Lower clipping threshold in standard deviations.</param>
    /// <param name="sigmaHigh">Upper clipping threshold in standard deviations.</param>
    /// <param name="maxIter">Maximum number of clipping iterations.</param>
    /// <returns>A <see cref="ClipResult"/> with the mask, mean, and std of retained points.</returns>
    public static ClipResult Clip(double[] data, double sigmaLow = 3.0, double sigmaHigh = 3.0, int maxIter = 10)
    {
        if (data == null) throw new ArgumentNullException(nameof(data));
        if (data.Length == 0) throw new ArgumentException("Data must not be empty.", nameof(data));
        if (sigmaLow <= 0) throw new ArgumentOutOfRangeException(nameof(sigmaLow), "Must be positive.");
        if (sigmaHigh <= 0) throw new ArgumentOutOfRangeException(nameof(sigmaHigh), "Must be positive.");
        if (maxIter < 1) throw new ArgumentOutOfRangeException(nameof(maxIter), "Must be at least 1.");

        bool[] mask = new bool[data.Length];
        for (int i = 0; i < data.Length; i++)
            mask[i] = true;

        double mean = 0;
        double std = 0;
        int retained = data.Length;
        int iteration = 0;

        for (iteration = 0; iteration < maxIter; iteration++)
        {
            // Compute mean of retained points
            mean = 0;
            retained = 0;
            for (int i = 0; i < data.Length; i++)
            {
                if (mask[i])
                {
                    mean += data[i];
                    retained++;
                }
            }

            if (retained < 2) break;
            mean /= retained;

            // Compute std of retained points
            double variance = 0;
            for (int i = 0; i < data.Length; i++)
            {
                if (mask[i])
                {
                    double d = data[i] - mean;
                    variance += d * d;
                }
            }
            std = Math.Sqrt(variance / retained);

            if (std < 1e-15) break;

            // Clip points outside [mean - sigmaLow*std, mean + sigmaHigh*std]
            double lower = mean - sigmaLow * std;
            double upper = mean + sigmaHigh * std;

            int clippedThisRound = 0;
            for (int i = 0; i < data.Length; i++)
            {
                if (mask[i] && (data[i] < lower || data[i] > upper))
                {
                    mask[i] = false;
                    clippedThisRound++;
                }
            }

            // Converged: no points clipped this iteration
            if (clippedThisRound == 0)
            {
                iteration++;
                break;
            }
        }

        // Final count
        retained = 0;
        for (int i = 0; i < data.Length; i++)
            if (mask[i]) retained++;

        int clipped = data.Length - retained;

        return new ClipResult(mask, mean, std, retained, clipped, iteration);
    }

    /// <summary>
    /// Convenience: symmetric sigma-clipping with the same threshold for both tails.
    /// </summary>
    public static ClipResult Clip(double[] data, double sigma, int maxIter = 10)
    {
        return Clip(data, sigma, sigma, maxIter);
    }

    /// <summary>
    /// Returns only the retained values (not clipped) from the data.
    /// </summary>
    public static double[] Apply(double[] data, double sigmaLow = 3.0, double sigmaHigh = 3.0, int maxIter = 10)
    {
        var result = Clip(data, sigmaLow, sigmaHigh, maxIter);
        var retained = new List<double>();
        for (int i = 0; i < data.Length; i++)
        {
            if (result.Mask[i])
                retained.Add(data[i]);
        }
        return retained.ToArray();
    }
}
