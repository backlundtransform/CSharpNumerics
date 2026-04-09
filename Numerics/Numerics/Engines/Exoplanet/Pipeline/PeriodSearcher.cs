using System;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public static class PeriodSearcher
{
    public static PeriodSearchResult Search(LightCurve lc, TransitDetectionConfig config)
    {
        return Search(lc, config.MinPeriodDays, config.MaxPeriodDays, config.PeriodSearchMethod, config.NumTrialPeriods);
    }

    public static PeriodSearchResult Search(LightCurve lc, double minPeriod, double maxPeriod,
        PeriodSearchMethod method, int numTrialPeriods = 5000)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));

        var times = new VectorN(lc.Time);
        var values = new VectorN(lc.Flux);

        switch (method)
        {
            case PeriodSearchMethod.BLS:
                return RunBLS(times, values, minPeriod, maxPeriod, numTrialPeriods);

            case PeriodSearchMethod.LombScargle:
                return RunLombScargle(times, values, minPeriod, maxPeriod, numTrialPeriods);

            case PeriodSearchMethod.Both:
                var bls = RunBLS(times, values, minPeriod, maxPeriod, numTrialPeriods);
                var ls = RunLombScargle(times, values, minPeriod, maxPeriod, numTrialPeriods);
                // Return BLS if it has a stronger detection, otherwise LS
                return bls.BestPower >= ls.BestPower ? bls : ls;

            default:
                return RunBLS(times, values, minPeriod, maxPeriod, numTrialPeriods);
        }
    }

    private static PeriodSearchResult RunBLS(VectorN times, VectorN values,
        double minPeriod, double maxPeriod, int numTrialPeriods)
    {
        var result = BoxFittingLeastSquares.Compute(times, values, minPeriod, maxPeriod, numTrialPeriods);

        double fap = EstimateBLSFAP(result.BestPower, result.Power);

        return new PeriodSearchResult(
            result.BestPeriod,
            result.BestPower,
            result.BestEpoch,
            result.BestDepth,
            result.BestDurationFraction,
            fap,
            result.Periods,
            result.Power);
    }

    private static PeriodSearchResult RunLombScargle(VectorN times, VectorN values,
        double minPeriod, double maxPeriod, int numTrialPeriods)
    {
        var result = LombScarglePeriodogram.ComputeByPeriod(times, values, minPeriod, maxPeriod, numTrialPeriods);

        double fap = LombScarglePeriodogram.FalseAlarmProbability(result.BestPower, times.Length, numTrialPeriods);

        return new PeriodSearchResult(
            result.BestPeriod,
            result.BestPower,
            0.0, // Lomb-Scargle doesn't directly give epoch
            0.0, // nor depth
            0.0,
            fap,
            result.Periods,
            result.Power);
    }

    /// <summary>
    /// Estimate BLS false alarm probability via bootstrap percentile (simple approximation):
    /// fraction of spectrum values that exceed the peak.
    /// </summary>
    private static double EstimateBLSFAP(double peakPower, VectorN spectrum)
    {
        if (spectrum.Length == 0) return 1.0;

        double mean = 0;
        double std = 0;
        for (int i = 0; i < spectrum.Length; i++)
            mean += spectrum[i];
        mean /= spectrum.Length;

        for (int i = 0; i < spectrum.Length; i++)
        {
            double d = spectrum[i] - mean;
            std += d * d;
        }
        std = Math.Sqrt(std / spectrum.Length);

        if (std < 1e-15) return 1.0;

        double snr = (peakPower - mean) / std;
        // Approximate FAP from Gaussian tail
        return Math.Exp(-0.5 * snr * snr) * spectrum.Length;
    }
}
