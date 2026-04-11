using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// False alarm probability (FAP) estimation for periodogram peaks.
/// Provides both analytical (Baluev 2008) and bootstrap-based methods.
/// </summary>
public static class FalseAlarmProbability
{
    /// <summary>
    /// Analytical FAP using the Baluev (2008) approximation for Lomb-Scargle periodograms.
    /// FAP ≈ 1 − (1 − p_single)^M where p_single = exp(−z) and M ≈ numFrequencies.
    /// </summary>
    /// <param name="peakPower">The normalized power at the peak frequency.</param>
    /// <param name="numFrequencies">Number of independent frequencies tested.</param>
    /// <param name="numDataPoints">Number of observations in the time series.</param>
    /// <returns>False alarm probability in [0, 1].</returns>
    public static double AnalyticalFAP(double peakPower, int numFrequencies, int numDataPoints)
    {
        if (numFrequencies <= 0) throw new ArgumentOutOfRangeException(nameof(numFrequencies));
        if (numDataPoints <= 0) throw new ArgumentOutOfRangeException(nameof(numDataPoints));
        if (peakPower < 0) return 1.0;

        // Baluev (2008): effective number of independent frequencies
        // For well-sampled data, M_eff ≈ numFrequencies
        // Single-frequency survival probability
        double tau = CalculateBaluevTau(peakPower, numDataPoints);

        // Probability of exceeding z at a single frequency
        double pSingle = tau;

        if (pSingle >= 1.0) return 1.0;
        if (pSingle <= 0.0) return 0.0;

        // Multi-trial correction
        // FAP = 1 - (1 - p_single)^M
        // For numerical stability when pSingle is small:
        if (pSingle < 1e-10)
            return Math.Min(1.0, pSingle * numFrequencies);

        double fap = 1.0 - Math.Pow(1.0 - pSingle, numFrequencies);
        return Math.Max(0.0, Math.Min(1.0, fap));
    }

    /// <summary>
    /// Bootstrap-based FAP estimation. Shuffles the observed values randomly,
    /// recomputes the periodogram peak each time, and counts how often the
    /// shuffled peak exceeds the observed peak power.
    /// </summary>
    /// <param name="times">Observation times.</param>
    /// <param name="values">Observed values.</param>
    /// <param name="peakPower">The observed peak power to test.</param>
    /// <param name="nBootstrap">Number of bootstrap iterations.</param>
    /// <param name="minPeriod">Minimum period for the periodogram (default: auto).</param>
    /// <param name="maxPeriod">Maximum period for the periodogram (default: auto).</param>
    /// <param name="numFrequencies">Number of frequencies to evaluate (default: 1000).</param>
    /// <param name="seed">Random seed for reproducibility (null = random).</param>
    /// <returns>Bootstrap FAP: fraction of shuffled datasets with peak ≥ observed peak.</returns>
    public static double BootstrapFAP(
        double[] times, double[] values, double peakPower,
        int nBootstrap = 1000,
        double minPeriod = 0, double maxPeriod = 0,
        int numFrequencies = 500,
        int? seed = null)
    {
        if (times == null) throw new ArgumentNullException(nameof(times));
        if (values == null) throw new ArgumentNullException(nameof(values));
        if (times.Length != values.Length)
            throw new ArgumentException("Times and values must have the same length.");
        if (times.Length < 4) throw new ArgumentException("Need at least 4 data points.");
        if (nBootstrap < 1) throw new ArgumentOutOfRangeException(nameof(nBootstrap));

        int n = times.Length;
        var timesVec = new VectorN(times);

        // Auto-determine period range if not specified
        if (minPeriod <= 0)
        {
            double[] sortedTimes = (double[])times.Clone();
            Array.Sort(sortedTimes);
            double minDt = double.MaxValue;
            for (int i = 1; i < sortedTimes.Length; i++)
            {
                double dt = sortedTimes[i] - sortedTimes[i - 1];
                if (dt > 0 && dt < minDt) minDt = dt;
            }
            minPeriod = 2.0 * minDt;
        }
        if (maxPeriod <= 0)
            maxPeriod = (times[times.Length - 1] - times[0]) / 2.0;

        if (maxPeriod <= minPeriod)
            maxPeriod = minPeriod * 10;

        var rng = seed.HasValue ? new System.Random(seed.Value) : new System.Random();

        int exceedCount = 0;
        double[] shuffled = new double[n];

        for (int b = 0; b < nBootstrap; b++)
        {
            // Shuffle values (keep times fixed — destroys any periodic signal)
            Array.Copy(values, shuffled, n);
            FisherYatesShuffle(shuffled, rng);

            var shuffledVec = new VectorN(shuffled);
            var result = LombScarglePeriodogram.ComputeByPeriod(
                timesVec, shuffledVec, minPeriod, maxPeriod, numFrequencies);

            if (result.BestPower >= peakPower)
                exceedCount++;
        }

        return (double)exceedCount / nBootstrap;
    }

    /// <summary>
    /// Baluev (2008) single-frequency false alarm probability.
    /// τ(z) = (1 - 2z/N)^((N-3)/2) for the Lomb-Scargle statistic.
    /// </summary>
    private static double CalculateBaluevTau(double z, int N)
    {
        if (N <= 3) return Math.Exp(-z);

        double base_val = 1.0 - 2.0 * z / N;
        if (base_val <= 0) return 0.0;

        double exponent = (N - 3.0) / 2.0;
        return Math.Pow(base_val, exponent);
    }

    /// <summary>
    /// Fisher-Yates shuffle in-place.
    /// </summary>
    private static void FisherYatesShuffle(double[] array, System.Random rng)
    {
        for (int i = array.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            double tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }
    }
}
