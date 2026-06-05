using System;

namespace CSharpNumerics.Numerics.SignalProcessing.Wavelets;

/// <summary>
/// Critically-sampled discrete wavelet transform (DWT) with periodic (circular) boundary
/// handling. Each level convolves with the analysis filters and downsamples by two, splitting a
/// signal into an approximation (low-frequency) and a detail (high-frequency) band. Because the
/// filter bank is orthonormal, the transform is exactly invertible (see
/// <see cref="InverseWaveletTransform"/>).
/// </summary>
public class DiscreteWaveletTransform
{
    public DiscreteWaveletTransform(WaveletFamily family)
    {
        Family = family ?? throw new ArgumentNullException(nameof(family));
    }

    public WaveletFamily Family { get; }

    /// <summary>
    /// Single-level DWT: returns the approximation and detail coefficients, each of length N/2.
    /// </summary>
    public (double[] approximation, double[] detail) SingleLevel(double[] signal)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        int n = signal.Length;
        if (n < 2 || n % 2 != 0)
            throw new ArgumentException("Signal length must be even and at least 2.", nameof(signal));

        double[] lo = Family.LowPass;
        double[] hi = Family.HighPass;
        int l = lo.Length;
        int half = n / 2;

        var approx = new double[half];
        var detail = new double[half];

        for (int i = 0; i < half; i++)
        {
            double a = 0.0, d = 0.0;
            for (int k = 0; k < l; k++)
            {
                int idx = Mod((2 * i) + k, n);
                a += lo[k] * signal[idx];
                d += hi[k] * signal[idx];
            }
            approx[i] = a;
            detail[i] = d;
        }

        return (approx, detail);
    }

    /// <summary>
    /// N-level DWT. The signal length must be divisible by 2^levels. Detail bands are returned
    /// finest-first.
    /// </summary>
    public WaveletDecomposition Decompose(double[] signal, int levels)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (levels < 1) throw new ArgumentOutOfRangeException(nameof(levels), "Levels must be at least 1.");
        if (signal.Length % (1 << levels) != 0)
            throw new ArgumentException($"Signal length must be divisible by 2^levels ({1 << levels}).", nameof(signal));

        var details = new double[levels][];
        double[] current = (double[])signal.Clone();

        for (int level = 0; level < levels; level++)
        {
            var (approx, detail) = SingleLevel(current);
            details[level] = detail;
            current = approx;
        }

        return new WaveletDecomposition(current, details, Family);
    }

    internal static int Mod(int a, int n) => ((a % n) + n) % n;
}
