using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public static class TransitDetectionPipeline
{
    /// <summary>
    /// Detects transit candidates in a light curve using an iterative pipeline:
    /// preprocess → period search → phase fold → fit → validate → subtract → repeat.
    /// </summary>
    public static TransitCandidate[] Detect(LightCurve lc, TransitDetectionConfig config)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));
        if (config == null) throw new ArgumentNullException(nameof(config));

        // Step 1: Preprocess
        var processed = LightCurvePreprocessor.Preprocess(lc, config);

        var candidates = new List<TransitCandidate>();
        var residualLc = processed;

        for (int planet = 0; planet < config.MaxPlanets; planet++)
        {
            // Step 2: Period search
            PeriodSearchResult searchResult;
            try
            {
                searchResult = PeriodSearcher.Search(residualLc, config);
            }
            catch
            {
                break;
            }

            if (searchResult.BestPower <= 0 || searchResult.BestDepth <= 0)
                break;

            // Step 3: Fit transit model
            TransitFitResult fitResult;
            try
            {
                fitResult = TransitFitter.Fit(residualLc, searchResult.BestPeriod, searchResult.BestEpoch);
            }
            catch
            {
                break;
            }

            // Step 4: Phase-fold for the candidate
            var times = new VectorN(residualLc.Time);
            var values = new VectorN(residualLc.Flux);
            var (binCenters, binMeans) = PhaseFolding.FoldAndBin(
                times, values, fitResult.Parameters.Period, config.PhaseBins, fitResult.Parameters.Epoch);

            double[] phaseTime = new double[binCenters.Length];
            double[] phaseFlux = new double[binMeans.Length];
            for (int i = 0; i < binCenters.Length; i++)
            {
                phaseTime[i] = binCenters[i];
                phaseFlux[i] = binMeans[i];
            }
            var phaseFolded = LightCurve.FromArrays(phaseTime, phaseFlux);

            // Create candidate
            var candidate = new TransitCandidate(
                fitResult.Parameters,
                0.0,
                TransitDisposition.Unknown,
                phaseFolded,
                new TransitFeatureSet());

            // Step 5: Validate
            var validation = TransitValidator.Validate(candidate, residualLc, config);
            candidate.Score = validation.Score;

            if (!validation.IsValid)
                break;

            candidate.Disposition = validation.Score > 0.7
                ? TransitDisposition.Candidate
                : TransitDisposition.Unknown;

            candidates.Add(candidate);

            // Step 6: Subtract transit signal for multi-planet search
            residualLc = SubtractTransit(residualLc, fitResult);
        }

        return candidates.ToArray();
    }

    /// <summary>
    /// Subtracts the fitted transit model from the light curve to reveal additional planets.
    /// </summary>
    private static LightCurve SubtractTransit(LightCurve lc, TransitFitResult fit)
    {
        double period = fit.Parameters.Period;
        double epoch = fit.Parameters.Epoch;
        double depth = fit.Parameters.Depth;
        double duration = fit.Parameters.Duration; // in days
        double ingressFrac = duration > 0
            ? fit.Parameters.IngressDuration / duration
            : 0.2;
        if (double.IsNaN(ingressFrac) || double.IsInfinity(ingressFrac))
            ingressFrac = 0.2;

        double halfDur = duration / 2.0;
        double ingressWidth = halfDur * Math.Max(ingressFrac, 0.01);

        double[] correctedFlux = new double[lc.Length];
        for (int i = 0; i < lc.Length; i++)
        {
            double phase = ((lc.Time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;
            double absPhase = Math.Abs(phase); // in days

            double model;
            if (absPhase > halfDur)
                model = 1.0;
            else if (absPhase < halfDur - ingressWidth)
                model = 1.0 - depth;
            else
            {
                double t = (halfDur - absPhase) / ingressWidth;
                model = 1.0 - depth * t;
            }

            // Divide out the transit model to restore baseline
            correctedFlux[i] = model > 1e-10 ? lc.Flux[i] / model : lc.Flux[i];
        }

        return new LightCurve(
            (double[])lc.Time.Clone(), correctedFlux,
            (double[])lc.FluxError.Clone(), (int[])lc.QualityFlags.Clone(), lc.Metadata);
    }

    /// <summary>
    /// Computes the scatter (standard deviation) of the light curve flux.
    /// Uses MAD (median absolute deviation) for robustness against outliers.
    /// </summary>
    private static double ComputeScatter(LightCurve lc)
    {
        if (lc.Length < 5) return double.MaxValue;

        double[] sorted = new double[lc.Length];
        Array.Copy(lc.Flux, sorted, lc.Length);
        Array.Sort(sorted);

        double median = sorted.Length % 2 == 1
            ? sorted[sorted.Length / 2]
            : (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]) / 2.0;

        double[] absDevs = new double[lc.Length];
        for (int i = 0; i < lc.Length; i++)
            absDevs[i] = Math.Abs(sorted[i] - median);
        Array.Sort(absDevs);

        double mad = absDevs.Length % 2 == 1
            ? absDevs[absDevs.Length / 2]
            : (absDevs[absDevs.Length / 2 - 1] + absDevs[absDevs.Length / 2]) / 2.0;

        return 1.4826 * mad;
    }
}
