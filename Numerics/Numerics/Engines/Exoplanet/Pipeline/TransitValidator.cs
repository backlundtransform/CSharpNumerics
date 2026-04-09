using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Exoplanet.Data;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public class ValidationResult
{
    public bool IsValid { get; }
    public string[] Warnings { get; }
    public double Score { get; }
    public double Snr { get; }
    public double OddEvenDepthRatio { get; }

    public ValidationResult(bool isValid, string[] warnings, double score, double snr, double oddEvenDepthRatio)
    {
        IsValid = isValid;
        Warnings = warnings;
        Score = score;
        Snr = snr;
        OddEvenDepthRatio = oddEvenDepthRatio;
    }
}

public static class TransitValidator
{
    public static ValidationResult Validate(TransitCandidate candidate, LightCurve lc, TransitDetectionConfig config)
    {
        if (candidate == null) throw new ArgumentNullException(nameof(candidate));
        if (lc == null) throw new ArgumentNullException(nameof(lc));
        if (config == null) throw new ArgumentNullException(nameof(config));

        var warnings = new List<string>();
        bool isValid = true;
        double score = 1.0;

        double depth = candidate.Parameters.Depth;
        double period = candidate.Parameters.Period;
        double epoch = candidate.Parameters.Epoch;

        // 1. SNR check: depth / scatter
        double scatter = ComputeOutOfTransitScatter(lc, period, epoch, candidate.Parameters.Duration);
        double snr = scatter > 0 ? depth / scatter : 0;

        if (snr < config.SnrThreshold)
        {
            isValid = false;
            warnings.Add($"SNR ({snr:F1}) below threshold ({config.SnrThreshold}).");
            score *= 0.3;
        }
        else
        {
            score *= Math.Min(1.0, snr / (config.SnrThreshold * 3));
        }

        // 2. Depth check
        double depthPpm = depth * 1e6;
        if (depthPpm < config.MinTransitDepthPpm)
        {
            isValid = false;
            warnings.Add($"Transit depth ({depthPpm:F0} ppm) below minimum ({config.MinTransitDepthPpm} ppm).");
            score *= 0.3;
        }

        // 3. Odd/even transit depth test
        double oddEvenRatio = ComputeOddEvenDepthRatio(lc, period, epoch, candidate.Parameters.Duration);

        if (oddEvenRatio < 0.5 || oddEvenRatio > 2.0)
        {
            warnings.Add($"Odd/even depth ratio ({oddEvenRatio:F2}) suggests eclipsing binary.");
            score *= 0.5;
        }

        // 4. Duration reasonability: transit duration / period should be < ~0.15
        double durFrac = candidate.Parameters.Duration / period;
        if (durFrac > 0.2)
        {
            warnings.Add($"Transit duration fraction ({durFrac:F3}) is unusually large.");
            score *= 0.7;
        }

        // 5. Period boundaries
        if (period < config.MinPeriodDays || period > config.MaxPeriodDays)
        {
            isValid = false;
            warnings.Add($"Period ({period:F4} days) outside search range.");
            score *= 0.1;
        }

        // 6. Secondary eclipse check
        double secondaryDepth = ComputeSecondaryEclipseDepth(lc, period, epoch);
        if (secondaryDepth > depth * 0.3)
        {
            warnings.Add($"Significant secondary eclipse detected (depth ratio {secondaryDepth / depth:F2}), possible eclipsing binary.");
            score *= 0.5;
        }

        score = Math.Max(0, Math.Min(1, score));

        return new ValidationResult(isValid, warnings.ToArray(), score, snr, oddEvenRatio);
    }

    private static double ComputeOutOfTransitScatter(LightCurve lc, double period, double epoch, double duration)
    {
        var outOfTransit = new List<double>();
        double halfDur = duration / 2.0;

        for (int i = 0; i < lc.Length; i++)
        {
            double phase = ((lc.Time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;

            if (Math.Abs(phase) > halfDur * 1.5)
                outOfTransit.Add(lc.Flux[i]);
        }

        if (outOfTransit.Count < 5) return double.MaxValue;

        double mean = 0;
        for (int i = 0; i < outOfTransit.Count; i++)
            mean += outOfTransit[i];
        mean /= outOfTransit.Count;

        double variance = 0;
        for (int i = 0; i < outOfTransit.Count; i++)
        {
            double d = outOfTransit[i] - mean;
            variance += d * d;
        }

        return Math.Sqrt(variance / outOfTransit.Count);
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

        double oddMean = 0, evenMean = 0;
        for (int i = 0; i < oddDepths.Count; i++) oddMean += oddDepths[i];
        for (int i = 0; i < evenDepths.Count; i++) evenMean += evenDepths[i];
        oddMean /= oddDepths.Count;
        evenMean /= evenDepths.Count;

        if (evenMean < 1e-10) return 1.0;
        return oddMean / evenMean;
    }

    private static double ComputeSecondaryEclipseDepth(LightCurve lc, double period, double epoch)
    {
        // Look at phase 0.5 (opposite side of orbit)
        double secondaryEpoch = epoch + period / 2.0;
        var inSecondary = new List<double>();
        var outSecondary = new List<double>();

        double halfDur = period * 0.025; // assume secondary ~ same duration

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

        double meanIn = 0, meanOut = 0;
        for (int i = 0; i < inSecondary.Count; i++) meanIn += inSecondary[i];
        for (int i = 0; i < outSecondary.Count; i++) meanOut += outSecondary[i];
        meanIn /= inSecondary.Count;
        meanOut /= outSecondary.Count;

        return Math.Max(0, meanOut - meanIn);
    }
}
