using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Lomb–Scargle periodogram for spectral analysis of unevenly sampled time series.
/// Implements the classical Lomb (1976) / Scargle (1982) algorithm with the
/// Townsend (2010) time-offset τ for numerical stability.
/// </summary>
public static class LombScarglePeriodogram
{
    /// <summary>
    /// Compute the Lomb–Scargle periodogram at specified angular frequencies.
    /// </summary>
    /// <param name="times">Observation times (not necessarily evenly spaced).</param>
    /// <param name="values">Observed values at each time.</param>
    /// <param name="angularFrequencies">Angular frequencies ω at which to evaluate power.</param>
    /// <param name="centerData">If true (default), subtract the mean before computing.</param>
    public static PeriodogramResult Compute(
        VectorN times,
        VectorN values,
        VectorN angularFrequencies,
        bool centerData = true)
    {
        if (times.Length != values.Length)
            throw new ArgumentException("times and values must have the same length.");
        if (times.Length < 2)
            throw new ArgumentException("At least 2 data points are required.");

        int n = times.Length;
        int nf = angularFrequencies.Length;
        double[] t = times.Values;
        double[] y = (double[])values.Values.Clone();

        if (centerData)
        {
            double mean = 0;
            for (int i = 0; i < n; i++) mean += y[i];
            mean /= n;
            for (int i = 0; i < n; i++) y[i] -= mean;
        }

        double[] power = new double[nf];
        double[] freqs = new double[nf];

        for (int k = 0; k < nf; k++)
        {
            double omega = angularFrequencies[k];
            freqs[k] = omega / (2.0 * Math.PI);

            // Compute τ: tan(2ωτ) = Σsin(2ωtᵢ) / Σcos(2ωtᵢ)
            double sinSum = 0, cosSum = 0;
            for (int i = 0; i < n; i++)
            {
                double arg = 2.0 * omega * t[i];
                sinSum += Math.Sin(arg);
                cosSum += Math.Cos(arg);
            }
            double tau = Math.Atan2(sinSum, cosSum) / (2.0 * omega);

            // Compute power using shifted times
            double ycNum = 0, ycDen = 0;
            double ysNum = 0, ysDen = 0;

            for (int i = 0; i < n; i++)
            {
                double phase = omega * (t[i] - tau);
                double c = Math.Cos(phase);
                double s = Math.Sin(phase);

                ycNum += y[i] * c;
                ycDen += c * c;
                ysNum += y[i] * s;
                ysDen += s * s;
            }

            power[k] = 0.5 * (
                (ycDen > 0 ? (ycNum * ycNum) / ycDen : 0) +
                (ysDen > 0 ? (ysNum * ysNum) / ysDen : 0));
        }

        return new PeriodogramResult(new VectorN(freqs), new VectorN(power));
    }

    /// <summary>
    /// Compute the periodogram over a linearly spaced frequency range.
    /// </summary>
    /// <param name="times">Observation times.</param>
    /// <param name="values">Observed values.</param>
    /// <param name="minFrequency">Minimum frequency (Hz, cycles per time unit).</param>
    /// <param name="maxFrequency">Maximum frequency (Hz).</param>
    /// <param name="numFrequencies">Number of frequency bins.</param>
    /// <param name="centerData">If true (default), subtract the mean.</param>
    public static PeriodogramResult Compute(
        VectorN times,
        VectorN values,
        double minFrequency,
        double maxFrequency,
        int numFrequencies,
        bool centerData = true)
    {
        if (minFrequency <= 0 || maxFrequency <= minFrequency)
            throw new ArgumentException("Frequencies must satisfy 0 < minFrequency < maxFrequency.");
        if (numFrequencies < 1)
            throw new ArgumentException("numFrequencies must be >= 1.");

        double[] omega = new double[numFrequencies];
        double step = (maxFrequency - minFrequency) / Math.Max(numFrequencies - 1, 1);
        for (int i = 0; i < numFrequencies; i++)
            omega[i] = 2.0 * Math.PI * (minFrequency + i * step);

        return Compute(times, values, new VectorN(omega), centerData);
    }

    /// <summary>
    /// Compute the periodogram over a range of trial periods.
    /// </summary>
    /// <param name="times">Observation times.</param>
    /// <param name="values">Observed values.</param>
    /// <param name="minPeriod">Minimum trial period.</param>
    /// <param name="maxPeriod">Maximum trial period.</param>
    /// <param name="numPeriods">Number of trial periods (linearly spaced in period).</param>
    /// <param name="centerData">If true, subtract the mean.</param>
    public static PeriodogramResult ComputeByPeriod(
        VectorN times,
        VectorN values,
        double minPeriod,
        double maxPeriod,
        int numPeriods,
        bool centerData = true)
    {
        if (minPeriod <= 0 || maxPeriod <= minPeriod)
            throw new ArgumentException("Periods must satisfy 0 < minPeriod < maxPeriod.");
        if (numPeriods < 1)
            throw new ArgumentException("numPeriods must be >= 1.");

        double step = (maxPeriod - minPeriod) / Math.Max(numPeriods - 1, 1);
        double[] omega = new double[numPeriods];
        for (int i = 0; i < numPeriods; i++)
        {
            double period = minPeriod + i * step;
            omega[i] = 2.0 * Math.PI / period;
        }

        return Compute(times, values, new VectorN(omega), centerData);
    }

    /// <summary>
    /// Estimate false-alarm probability for a given power level using the
    /// Baluev (2008) approximation: FAP ≈ 1 − (1 − exp(−z))^M, where M = N_eff.
    /// </summary>
    /// <param name="power">The observed spectral power.</param>
    /// <param name="numObservations">Number of data points.</param>
    /// <param name="numFrequencies">Number of independent frequencies tested.</param>
    public static double FalseAlarmProbability(double power, int numObservations, int numFrequencies)
    {
        double pSingle = Math.Exp(-power);
        if (pSingle >= 1.0) return 1.0;
        return 1.0 - Math.Pow(1.0 - pSingle, numFrequencies);
    }
}
