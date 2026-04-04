using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Phase folding utilities for periodic time-series analysis.
/// Folds observation times onto a single cycle of a trial period.
/// </summary>
public static class PhaseFolding
{
    /// <summary>
    /// Fold a time series on a given period, returning phases in [0, 1).
    /// </summary>
    /// <param name="times">Observation times.</param>
    /// <param name="values">Observed values.</param>
    /// <param name="period">Trial period to fold on.</param>
    /// <param name="epoch">Reference epoch (time of phase 0). Default 0.</param>
    /// <returns>Tuple of (phases, values) sorted by phase.</returns>
    public static (VectorN Phases, VectorN Values) Fold(
        VectorN times,
        VectorN values,
        double period,
        double epoch = 0)
    {
        if (times.Length != values.Length)
            throw new ArgumentException("times and values must have the same length.");
        if (period <= 0)
            throw new ArgumentException("period must be positive.");

        int n = times.Length;
        double[] t = times.Values;
        double[] y = values.Values;

        double[] phases = new double[n];
        int[] indices = new int[n];
        for (int i = 0; i < n; i++)
        {
            phases[i] = (((t[i] - epoch) / period) % 1.0 + 1.0) % 1.0;
            indices[i] = i;
        }

        // Sort by phase
        Array.Sort(phases, indices);

        double[] sortedValues = new double[n];
        for (int i = 0; i < n; i++)
            sortedValues[i] = y[indices[i]];

        return (new VectorN(phases), new VectorN(sortedValues));
    }

    /// <summary>
    /// Fold and bin a time series, averaging values within each phase bin.
    /// </summary>
    /// <param name="times">Observation times.</param>
    /// <param name="values">Observed values.</param>
    /// <param name="period">Trial period to fold on.</param>
    /// <param name="numBins">Number of phase bins.</param>
    /// <param name="epoch">Reference epoch. Default 0.</param>
    /// <returns>Tuple of (binCenters, binMeans) — bins with no data get NaN.</returns>
    public static (VectorN BinCenters, VectorN BinMeans) FoldAndBin(
        VectorN times,
        VectorN values,
        double period,
        int numBins,
        double epoch = 0)
    {
        if (times.Length != values.Length)
            throw new ArgumentException("times and values must have the same length.");
        if (period <= 0)
            throw new ArgumentException("period must be positive.");
        if (numBins < 1)
            throw new ArgumentException("numBins must be >= 1.");

        int n = times.Length;
        double[] t = times.Values;
        double[] y = values.Values;

        double[] binSum = new double[numBins];
        int[] binCount = new int[numBins];

        for (int i = 0; i < n; i++)
        {
            double phase = (((t[i] - epoch) / period) % 1.0 + 1.0) % 1.0;
            int bin = (int)(phase * numBins);
            if (bin >= numBins) bin = numBins - 1;
            binSum[bin] += y[i];
            binCount[bin]++;
        }

        double[] centers = new double[numBins];
        double[] means = new double[numBins];
        for (int i = 0; i < numBins; i++)
        {
            centers[i] = (i + 0.5) / numBins;
            means[i] = binCount[i] > 0 ? binSum[i] / binCount[i] : double.NaN;
        }

        return (new VectorN(centers), new VectorN(means));
    }

    /// <summary>
    /// Fold and bin, also returning standard deviation within each bin.
    /// </summary>
    public static (VectorN BinCenters, VectorN BinMeans, VectorN BinStdDevs) FoldAndBinWithErrors(
        VectorN times,
        VectorN values,
        double period,
        int numBins,
        double epoch = 0)
    {
        if (times.Length != values.Length)
            throw new ArgumentException("times and values must have the same length.");
        if (period <= 0)
            throw new ArgumentException("period must be positive.");
        if (numBins < 1)
            throw new ArgumentException("numBins must be >= 1.");

        int n = times.Length;
        double[] t = times.Values;
        double[] y = values.Values;

        double[] binSum = new double[numBins];
        double[] binSumSq = new double[numBins];
        int[] binCount = new int[numBins];

        for (int i = 0; i < n; i++)
        {
            double phase = (((t[i] - epoch) / period) % 1.0 + 1.0) % 1.0;
            int bin = (int)(phase * numBins);
            if (bin >= numBins) bin = numBins - 1;
            binSum[bin] += y[i];
            binSumSq[bin] += y[i] * y[i];
            binCount[bin]++;
        }

        double[] centers = new double[numBins];
        double[] means = new double[numBins];
        double[] stds = new double[numBins];
        for (int i = 0; i < numBins; i++)
        {
            centers[i] = (i + 0.5) / numBins;
            if (binCount[i] > 0)
            {
                double m = binSum[i] / binCount[i];
                means[i] = m;
                stds[i] = binCount[i] > 1
                    ? Math.Sqrt((binSumSq[i] / binCount[i]) - m * m)
                    : 0;
            }
            else
            {
                means[i] = double.NaN;
                stds[i] = double.NaN;
            }
        }

        return (new VectorN(centers), new VectorN(means), new VectorN(stds));
    }
}
