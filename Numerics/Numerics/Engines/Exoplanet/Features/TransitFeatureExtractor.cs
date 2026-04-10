using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Pipeline;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace CSharpNumerics.Engines.Exoplanet.Features;

public static class TransitFeatureExtractor
{
    /// <summary>
    /// Extracts a full set of transit-specific features from a candidate and its parent light curve.
    /// </summary>
    public static TransitFeatureSet Extract(TransitCandidate candidate, LightCurve lc)
    {
        if (candidate == null) throw new ArgumentNullException(nameof(candidate));
        if (lc == null) throw new ArgumentNullException(nameof(lc));

        var features = new TransitFeatureSet();
        var p = candidate.Parameters;

        // Basic transit parameters
        features[TransitFeatureSet.FeatureNames.Depth] = p.Depth;
        features[TransitFeatureSet.FeatureNames.Duration] = p.Duration;
        features[TransitFeatureSet.FeatureNames.Period] = p.Period;

        // Ingress/egress ratio
        double ingressEgressRatio = p.Duration > 0 && p.IngressDuration > 0
            ? p.IngressDuration / (p.Duration / 2.0)
            : 0.0;
        features[TransitFeatureSet.FeatureNames.IngressEgressRatio] = ingressEgressRatio;

        // V-shape metric: ratio of ingress+egress to total duration.
        // U-shaped transits (flat bottom) have low values; V-shaped (grazing) have high values.
        double vShape = p.Duration > 0
            ? Math.Min(1.0, 2.0 * p.IngressDuration / p.Duration)
            : 0.0;
        features[TransitFeatureSet.FeatureNames.VShapeMetric] = vShape;

        // BLS SNR
        double blsSnr = ComputeBlsSnr(lc, p.Period, p.Epoch, p.Duration);
        features[TransitFeatureSet.FeatureNames.SnrBls] = blsSnr;

        // Odd/even depth ratio
        double oddEvenRatio = ComputeOddEvenDepthRatio(lc, p.Period, p.Epoch, p.Duration);
        features[TransitFeatureSet.FeatureNames.OddEvenRatio] = oddEvenRatio;

        // Secondary eclipse depth
        double secondaryDepth = ComputeSecondaryEclipseDepth(lc, p.Period, p.Epoch);
        features[TransitFeatureSet.FeatureNames.SecondaryDepth] = secondaryDepth;

        // In-transit and out-of-transit scatter
        var (scatterIn, scatterOut) = ComputeScatters(lc, p.Period, p.Epoch, p.Duration);
        features[TransitFeatureSet.FeatureNames.ScatterInTransit] = scatterIn;
        features[TransitFeatureSet.FeatureNames.ScatterOutTransit] = scatterOut;

        // Limb darkening coefficients from quadratic fit to phase-folded transit shape
        var (u1, u2) = FitLimbDarkeningCoefficients(candidate);
        features[TransitFeatureSet.FeatureNames.LimbDarkeningU1] = u1;
        features[TransitFeatureSet.FeatureNames.LimbDarkeningU2] = u2;

        return features;
    }

    /// <summary>
    /// Computes the BLS signal-to-noise ratio for the transit.
    /// SNR = depth * sqrt(numInTransit) / scatter_out_of_transit
    /// </summary>
    private static double ComputeBlsSnr(LightCurve lc, double period, double epoch, double duration)
    {
        double halfDur = duration / 2.0;
        var inTransit = new List<double>();
        var outTransit = new List<double>();

        for (int i = 0; i < lc.Length; i++)
        {
            double phase = ((lc.Time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;

            if (Math.Abs(phase) <= halfDur)
                inTransit.Add(lc.Flux[i]);
            else
                outTransit.Add(lc.Flux[i]);
        }

        if (inTransit.Count < 2 || outTransit.Count < 5) return 0;

        double meanIn = Mean(inTransit);
        double meanOut = Mean(outTransit);
        double depth = meanOut - meanIn;

        double scatter = Std(outTransit);
        if (scatter < 1e-15) return 0;

        return depth * Math.Sqrt(inTransit.Count) / scatter;
    }

    private static double ComputeOddEvenDepthRatio(LightCurve lc, double period, double epoch, double duration)
    {
        var oddDepths = new List<double>();
        var evenDepths = new List<double>();
        double halfDur = duration / 2.0;

        for (int i = 0; i < lc.Length; i++)
        {
            double timeSinceEpoch = lc.Time[i] - epoch;
            double phase = (timeSinceEpoch % period + period) % period;
            if (phase > period / 2.0) phase -= period;

            if (Math.Abs(phase) <= halfDur)
            {
                int transitNumber = (int)Math.Round(timeSinceEpoch / period);
                if (transitNumber % 2 == 0)
                    evenDepths.Add(1.0 - lc.Flux[i]);
                else
                    oddDepths.Add(1.0 - lc.Flux[i]);
            }
        }

        if (oddDepths.Count < 2 || evenDepths.Count < 2) return 1.0;

        double oddMean = Mean(oddDepths);
        double evenMean = Mean(evenDepths);

        if (evenMean < 1e-10) return 1.0;
        return oddMean / evenMean;
    }

    private static double ComputeSecondaryEclipseDepth(LightCurve lc, double period, double epoch)
    {
        double secondaryEpoch = epoch + period / 2.0;
        var inSecondary = new List<double>();
        var outSecondary = new List<double>();
        double halfDur = period * 0.025;

        for (int i = 0; i < lc.Length; i++)
        {
            double phase = ((lc.Time[i] - secondaryEpoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;

            if (Math.Abs(phase) <= halfDur)
                inSecondary.Add(lc.Flux[i]);
            else if (Math.Abs(phase) > halfDur * 3)
                outSecondary.Add(lc.Flux[i]);
        }

        if (inSecondary.Count < 3 || outSecondary.Count < 3) return 0;

        double meanIn = Mean(inSecondary);
        double meanOut = Mean(outSecondary);

        return Math.Max(0, meanOut - meanIn);
    }

    private static (double scatterIn, double scatterOut) ComputeScatters(
        LightCurve lc, double period, double epoch, double duration)
    {
        var inTransit = new List<double>();
        var outTransit = new List<double>();
        double halfDur = duration / 2.0;

        for (int i = 0; i < lc.Length; i++)
        {
            double phase = ((lc.Time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;

            if (Math.Abs(phase) <= halfDur)
                inTransit.Add(lc.Flux[i]);
            else if (Math.Abs(phase) > halfDur * 1.5)
                outTransit.Add(lc.Flux[i]);
        }

        double scatterIn = inTransit.Count >= 3 ? Std(inTransit) : 0;
        double scatterOut = outTransit.Count >= 3 ? Std(outTransit) : 0;
        return (scatterIn, scatterOut);
    }

    /// <summary>
    /// Estimates quadratic limb darkening coefficients (u1, u2) from the phase-folded transit shape.
    /// Uses a simple least-squares fit of the transit profile curvature.
    /// </summary>
    private static (double u1, double u2) FitLimbDarkeningCoefficients(TransitCandidate candidate)
    {
        var pf = candidate.PhaseFoldedCurve;
        if (pf == null || pf.Length < 10) return (0.0, 0.0);

        var p = candidate.Parameters;
        double depth = p.Depth;
        if (depth < 1e-10) return (0.0, 0.0);

        // Map in-transit points to mu (cosine of angle from disk center)
        // and fit quadratic LD from the transit shape deviation from a flat-bottom model.
        double halfDurPhase = p.Duration / (2.0 * p.Period);
        if (halfDurPhase <= 0) return (0.0, 0.0);

        // Collect in-transit data points
        var muList = new List<double>();
        var normFluxList = new List<double>();

        for (int i = 0; i < pf.Length; i++)
        {
            double phase = pf.Time[i];
            // Normalize phase to [-0.5, 0.5]
            if (phase > 0.5) phase -= 1.0;
            double absPhase = Math.Abs(phase);

            if (absPhase < halfDurPhase && absPhase > 0)
            {
                // Map position within transit to mu = sqrt(1 - (r/R*)^2)
                // r/R* ≈ absPhase / halfDurPhase (linear approximation for small planets)
                double rNorm = absPhase / halfDurPhase;
                if (rNorm < 1.0)
                {
                    double mu = Math.Sqrt(1.0 - rNorm * rNorm);
                    muList.Add(mu);
                    // Observed flux drop normalized by depth
                    normFluxList.Add((1.0 - pf.Flux[i]) / depth);
                }
            }
        }

        if (muList.Count < 4) return (0.4, 0.2); // Default solar-like values

        // Fit: f(mu) ≈ 1 - u1*(1-mu) - u2*(1-mu)^2
        // Using normal equations for [u1, u2]
        double s00 = 0, s01 = 0, s11 = 0, r0 = 0, r1 = 0;
        for (int i = 0; i < muList.Count; i++)
        {
            double x0 = 1.0 - muList[i];
            double x1 = x0 * x0;
            double y = 1.0 - normFluxList[i]; // observed intensity relative to center
            s00 += x0 * x0;
            s01 += x0 * x1;
            s11 += x1 * x1;
            r0 += x0 * y;
            r1 += x1 * y;
        }

        double det = s00 * s11 - s01 * s01;
        if (Math.Abs(det) < 1e-20) return (0.4, 0.2);

        double u1 = (s11 * r0 - s01 * r1) / det;
        double u2 = (s00 * r1 - s01 * r0) / det;

        // Guard against NaN/Inf from ill-conditioned fit
        if (double.IsNaN(u1) || double.IsInfinity(u1) ||
            double.IsNaN(u2) || double.IsInfinity(u2))
            return (0.4, 0.2);

        // Clamp to physically reasonable range
        u1 = Math.Max(-1.0, Math.Min(2.0, u1));
        u2 = Math.Max(-1.0, Math.Min(2.0, u2));

        return (u1, u2);
    }

    private static double Mean(List<double> data)
    {
        double sum = 0;
        for (int i = 0; i < data.Count; i++)
            sum += data[i];
        return sum / data.Count;
    }

    private static double Std(List<double> data)
    {
        double mean = Mean(data);
        double variance = 0;
        for (int i = 0; i < data.Count; i++)
        {
            double d = data[i] - mean;
            variance += d * d;
        }
        return Math.Sqrt(variance / data.Count);
    }
}
