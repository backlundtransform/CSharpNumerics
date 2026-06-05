using System;
using System.Linq;

namespace CSharpNumerics.Numerics.SignalProcessing.Wavelets;

/// <summary>Thresholding rule applied to wavelet detail coefficients.</summary>
public enum ThresholdType
{
    /// <summary>Shrink towards zero: <c>sign(x)·max(|x|-λ, 0)</c>. Smoother result.</summary>
    Soft,

    /// <summary>Keep or kill: <c>x</c> if <c>|x| &gt; λ</c>, else 0.</summary>
    Hard
}

/// <summary>Strategy for choosing the denoising threshold.</summary>
public enum ThresholdRule
{
    /// <summary>Universal threshold <c>λ = σ√(2 ln N)</c> (Donoho &amp; Johnstone).</summary>
    VisuShrink,

    /// <summary>Per-level data-adaptive threshold <c>λ = σ² / σ_x</c> (Chang et al.).</summary>
    BayesShrink
}

/// <summary>
/// Wavelet shrinkage de-noising. The signal is decomposed, the detail coefficients are
/// thresholded to suppress noise (which spreads thinly across coefficients while the signal
/// concentrates into a few large ones), and the signal is reconstructed. The noise level is
/// estimated robustly from the finest detail band via the median absolute deviation.
/// </summary>
public static class WaveletDenoising
{
    /// <summary>
    /// De-noises <paramref name="signal"/> using the given wavelet, decomposition depth,
    /// thresholding type and rule.
    /// </summary>
    public static double[] Denoise(
        double[] signal,
        WaveletFamily family,
        int levels,
        ThresholdType type = ThresholdType.Soft,
        ThresholdRule rule = ThresholdRule.VisuShrink)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (family == null) throw new ArgumentNullException(nameof(family));

        var dwt = new DiscreteWaveletTransform(family);
        var decomposition = dwt.Decompose(signal, levels);

        // Robust noise estimate σ from the finest detail band (MAD / 0.6745).
        double sigma = RobustNoiseSigma(decomposition.Details[0]);
        int n = signal.Length;
        double universal = sigma * Math.Sqrt(2.0 * Math.Log(n));

        foreach (var band in decomposition.Details)
        {
            double threshold = rule == ThresholdRule.VisuShrink
                ? universal
                : BayesThreshold(band, sigma);

            for (int i = 0; i < band.Length; i++)
                band[i] = ApplyThreshold(band[i], threshold, type);
        }

        var inverse = new InverseWaveletTransform(family);
        return inverse.Reconstruct(decomposition);
    }

    private static double ApplyThreshold(double x, double threshold, ThresholdType type)
    {
        if (type == ThresholdType.Hard)
            return Math.Abs(x) > threshold ? x : 0.0;

        // Soft.
        double magnitude = Math.Abs(x) - threshold;
        if (magnitude <= 0.0) return 0.0;
        return Math.Sign(x) * magnitude;
    }

    private static double RobustNoiseSigma(double[] finestDetail)
    {
        if (finestDetail.Length == 0) return 0.0;
        var absValues = finestDetail.Select(Math.Abs).ToArray();
        double median = Median(absValues);
        return median / 0.6745;   // MAD-to-σ scaling for Gaussian noise
    }

    private static double BayesThreshold(double[] band, double sigma)
    {
        double noiseVar = sigma * sigma;

        double mean = band.Average();
        double observedVar = 0.0;
        foreach (double v in band) observedVar += (v - mean) * (v - mean);
        observedVar /= band.Length;

        double signalVar = observedVar - noiseVar;
        if (signalVar <= 0.0)
            return band.Select(Math.Abs).DefaultIfEmpty(0.0).Max();   // pure noise → kill the band

        return noiseVar / Math.Sqrt(signalVar);
    }

    private static double Median(double[] values)
    {
        var sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int mid = sorted.Length / 2;
        return sorted.Length % 2 == 0
            ? 0.5 * (sorted[mid - 1] + sorted[mid])
            : sorted[mid];
    }
}
