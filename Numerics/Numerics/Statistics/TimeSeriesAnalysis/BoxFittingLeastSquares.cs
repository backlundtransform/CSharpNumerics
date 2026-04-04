using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Box-fitting Least Squares (BLS) algorithm for periodic transit detection
/// in time-series photometry. Implements the Kovács, Zucker &amp; Mazeh (2002) method.
/// </summary>
public static class BoxFittingLeastSquares
{
    /// <summary>
    /// Search for periodic box-shaped (transit) signals.
    /// </summary>
    /// <param name="times">Observation times.</param>
    /// <param name="values">Observed flux/magnitude values.</param>
    /// <param name="periods">Array of trial periods to test.</param>
    /// <param name="minDurationFraction">Minimum transit duration as a fraction of the period (e.g. 0.01).</param>
    /// <param name="maxDurationFraction">Maximum transit duration as a fraction of the period (e.g. 0.15).</param>
    /// <param name="numBins">Number of phase bins for folding (default 200).</param>
    public static BLSResult Compute(
        VectorN times,
        VectorN values,
        VectorN periods,
        double minDurationFraction = 0.01,
        double maxDurationFraction = 0.15,
        int numBins = 200)
    {
        if (times.Length != values.Length)
            throw new ArgumentException("times and values must have the same length.");
        if (times.Length < 3)
            throw new ArgumentException("At least 3 data points are required.");
        if (minDurationFraction <= 0 || maxDurationFraction <= minDurationFraction || maxDurationFraction >= 1)
            throw new ArgumentException("Duration fractions must satisfy 0 < min < max < 1.");
        if (numBins < 5)
            throw new ArgumentException("numBins must be >= 5.");

        int n = times.Length;
        int np = periods.Length;
        double[] t = times.Values;
        double[] y = values.Values;

        double[] sr = new double[np];
        int bestIdx = 0;
        double bestSR = double.NegativeInfinity;
        double bestEpoch = 0;
        double bestDurFrac = 0;
        double bestDepth = 0;

        // Precompute mean
        double yMean = 0;
        for (int i = 0; i < n; i++) yMean += y[i];
        yMean /= n;

        int minBinWidth = Math.Max(1, (int)(minDurationFraction * numBins));
        int maxBinWidth = Math.Max(minBinWidth + 1, (int)(maxDurationFraction * numBins));
        if (maxBinWidth > numBins) maxBinWidth = numBins;

        for (int ip = 0; ip < np; ip++)
        {
            double period = periods[ip];
            if (period <= 0) continue;

            // Phase-fold and bin
            double[] binSum = new double[numBins];
            int[] binCount = new int[numBins];

            for (int i = 0; i < n; i++)
            {
                double phase = ((t[i] / period) % 1.0 + 1.0) % 1.0;
                int bin = (int)(phase * numBins);
                if (bin >= numBins) bin = numBins - 1;
                binSum[bin] += y[i];
                binCount[bin]++;
            }

            // Try all contiguous box positions and widths
            double periodBestSR = double.NegativeInfinity;
            int periodBestStart = 0;
            int periodBestWidth = minBinWidth;

            for (int width = minBinWidth; width <= maxBinWidth; width++)
            {
                // Initialize the first window
                double s = 0;
                int r = 0;
                for (int b = 0; b < width; b++)
                {
                    s += binSum[b];
                    r += binCount[b];
                }

                double srVal = ComputeSR(s, r, n, yMean);
                if (srVal > periodBestSR)
                {
                    periodBestSR = srVal;
                    periodBestStart = 0;
                    periodBestWidth = width;
                }

                // Slide the window
                for (int start = 1; start < numBins; start++)
                {
                    int removeIdx = start - 1;
                    int addIdx = (start + width - 1) % numBins;

                    s -= binSum[removeIdx];
                    r -= binCount[removeIdx];
                    s += binSum[addIdx];
                    r += binCount[addIdx];

                    srVal = ComputeSR(s, r, n, yMean);
                    if (srVal > periodBestSR)
                    {
                        periodBestSR = srVal;
                        periodBestStart = start;
                        periodBestWidth = width;
                    }
                }
            }

            sr[ip] = periodBestSR;

            if (periodBestSR > bestSR)
            {
                bestSR = periodBestSR;
                bestIdx = ip;
                bestDurFrac = (double)periodBestWidth / numBins;
                // Epoch: mid-point of the best box in time
                double phaseCenter = ((double)periodBestStart + periodBestWidth * 0.5) / numBins;
                bestEpoch = phaseCenter * period;

                // Depth: average in-transit minus average out-of-transit
                double inSum = 0;
                int inCount = 0;
                for (int b = 0; b < periodBestWidth; b++)
                {
                    int idx = (periodBestStart + b) % numBins;
                    inSum += binSum[idx];
                    inCount += binCount[idx];
                }
                double inMean = inCount > 0 ? inSum / inCount : yMean;
                double outSum = 0;
                int outCount = 0;
                for (int b = periodBestWidth; b < numBins; b++)
                {
                    int idx = (periodBestStart + b) % numBins;
                    outSum += binSum[idx];
                    outCount += binCount[idx];
                }
                double outMean = outCount > 0 ? outSum / outCount : yMean;
                bestDepth = outMean - inMean;
            }
        }

        return new BLSResult(
            periods,
            new VectorN(sr),
            bestIdx,
            bestEpoch,
            bestDurFrac,
            bestDepth);
    }

    /// <summary>
    /// Convenience overload: linearly spaced trial periods.
    /// </summary>
    public static BLSResult Compute(
        VectorN times,
        VectorN values,
        double minPeriod,
        double maxPeriod,
        int numPeriods,
        double minDurationFraction = 0.01,
        double maxDurationFraction = 0.15,
        int numBins = 200)
    {
        if (minPeriod <= 0 || maxPeriod <= minPeriod)
            throw new ArgumentException("Periods must satisfy 0 < minPeriod < maxPeriod.");

        double step = (maxPeriod - minPeriod) / Math.Max(numPeriods - 1, 1);
        double[] p = new double[numPeriods];
        for (int i = 0; i < numPeriods; i++)
            p[i] = minPeriod + i * step;

        return Compute(times, values, new VectorN(p), minDurationFraction, maxDurationFraction, numBins);
    }

    /// <summary>
    /// Compute Signal Residue: SR = sqrt(s² / (r · (n - r)))
    /// where s = in-transit sum minus r·mean, r = in-transit count, n = total.
    /// </summary>
    private static double ComputeSR(double s, int r, int n, double mean)
    {
        if (r <= 0 || r >= n) return 0;
        double sAdj = s - r * mean;
        return Math.Sqrt((sAdj * sAdj) / ((double)r * (n - r)));
    }
}
