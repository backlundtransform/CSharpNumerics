using System;

namespace CSharpNumerics.Statistics.Robust;

/// <summary>
/// Result of a RANSAC fitting procedure.
/// </summary>
public class RansacResult
{
    /// <summary>Best-fit model parameters (e.g. [intercept, slope] for a line).</summary>
    public double[] BestModel { get; }

    /// <summary>Boolean mask: true = inlier, false = outlier.</summary>
    public bool[] InlierMask { get; }

    /// <summary>Number of inliers for the best model.</summary>
    public int InlierCount { get; }

    /// <summary>Number of iterations performed.</summary>
    public int Iterations { get; }

    public RansacResult(double[] bestModel, bool[] inlierMask, int inlierCount, int iterations)
    {
        BestModel = bestModel;
        InlierMask = inlierMask;
        InlierCount = inlierCount;
        Iterations = iterations;
    }
}

/// <summary>
/// RANSAC (RANdom SAmple Consensus): an iterative algorithm for robust model
/// estimation in the presence of outliers.
///
/// At each iteration a minimal random subset is drawn, a model is fitted,
/// and all points within the residual threshold are counted as inliers.
/// The model with the most inliers wins.
/// </summary>
public static class Ransac
{
    /// <summary>
    /// Fits a line y = a + b·x using RANSAC.
    /// </summary>
    /// <param name="x">X-coordinates.</param>
    /// <param name="y">Y-coordinates.</param>
    /// <param name="residualThreshold">Maximum absolute residual for a point to be considered an inlier.</param>
    /// <param name="maxIterations">Maximum number of random sampling iterations.</param>
    /// <param name="seed">Random seed for reproducibility. Use -1 for non-deterministic.</param>
    /// <returns>A <see cref="RansacResult"/> containing the best line [intercept, slope] and the inlier mask.</returns>
    public static RansacResult FitLine(double[] x, double[] y, double residualThreshold = 1.0, int maxIterations = 1000, int seed = 42)
    {
        if (x == null) throw new ArgumentNullException(nameof(x));
        if (y == null) throw new ArgumentNullException(nameof(y));
        if (x.Length != y.Length) throw new ArgumentException("x and y must have the same length.");
        if (x.Length < 2) throw new ArgumentException("Need at least 2 data points.");
        if (residualThreshold <= 0) throw new ArgumentOutOfRangeException(nameof(residualThreshold), "Must be positive.");
        if (maxIterations < 1) throw new ArgumentOutOfRangeException(nameof(maxIterations), "Must be at least 1.");

        System.Random rng = seed >= 0 ? new System.Random(seed) : new System.Random();

        double[] bestModel = null;
        bool[] bestMask = null;
        int bestInlierCount = 0;
        int iterations = 0;

        for (iterations = 0; iterations < maxIterations; iterations++)
        {
            // Pick 2 distinct random indices
            int i1 = rng.Next(x.Length);
            int i2;
            do { i2 = rng.Next(x.Length); } while (i2 == i1);

            // Fit line through the two points
            double dx = x[i2] - x[i1];
            if (Math.Abs(dx) < 1e-15) continue;

            double slope = (y[i2] - y[i1]) / dx;
            double intercept = y[i1] - slope * x[i1];

            // Count inliers
            bool[] mask = new bool[x.Length];
            int inlierCount = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double residual = Math.Abs(y[i] - (intercept + slope * x[i]));
                if (residual <= residualThreshold)
                {
                    mask[i] = true;
                    inlierCount++;
                }
            }

            if (inlierCount > bestInlierCount)
            {
                bestInlierCount = inlierCount;
                bestModel = new[] { intercept, slope };
                bestMask = mask;
            }

            // Early exit if all points are inliers
            if (bestInlierCount == x.Length) break;
        }

        // Refit model on all inliers for a better estimate
        if (bestMask != null && bestInlierCount >= 2)
        {
            double sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
            int n = 0;
            for (int i = 0; i < x.Length; i++)
            {
                if (bestMask[i])
                {
                    sumX += x[i];
                    sumY += y[i];
                    sumXX += x[i] * x[i];
                    sumXY += x[i] * y[i];
                    n++;
                }
            }
            double denom = n * sumXX - sumX * sumX;
            if (Math.Abs(denom) > 1e-15)
            {
                double slope = (n * sumXY - sumX * sumY) / denom;
                double intercept = (sumY - slope * sumX) / n;
                bestModel = new[] { intercept, slope };

                // Recompute inlier mask with refined model
                bestInlierCount = 0;
                for (int i = 0; i < x.Length; i++)
                {
                    double residual = Math.Abs(y[i] - (intercept + slope * x[i]));
                    bestMask[i] = residual <= residualThreshold;
                    if (bestMask[i]) bestInlierCount++;
                }
            }
        }

        bestModel ??= new[] { 0.0, 0.0 };
        bestMask ??= new bool[x.Length];

        return new RansacResult(bestModel, bestMask, bestInlierCount, iterations + 1);
    }
}
