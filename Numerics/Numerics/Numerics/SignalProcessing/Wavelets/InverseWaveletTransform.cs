using System;

namespace CSharpNumerics.Numerics.SignalProcessing.Wavelets;

/// <summary>
/// Inverse discrete wavelet transform. Because the periodic filter bank is orthonormal, synthesis
/// is the transpose (adjoint) of analysis: each coefficient scatters its filtered contribution
/// back onto the reconstructed signal. Reconstruction is exact to machine precision.
/// </summary>
public class InverseWaveletTransform
{
    public InverseWaveletTransform(WaveletFamily family)
    {
        Family = family ?? throw new ArgumentNullException(nameof(family));
    }

    public WaveletFamily Family { get; }

    /// <summary>
    /// Single-level reconstruction from approximation and detail bands (each length N/2) back to
    /// the length-N signal.
    /// </summary>
    public double[] SingleLevel(double[] approximation, double[] detail)
    {
        if (approximation == null) throw new ArgumentNullException(nameof(approximation));
        if (detail == null) throw new ArgumentNullException(nameof(detail));
        if (approximation.Length != detail.Length)
            throw new ArgumentException("Approximation and detail bands must have equal length.");

        double[] lo = Family.LowPass;
        double[] hi = Family.HighPass;
        int l = lo.Length;
        int half = approximation.Length;
        int n = half * 2;

        var reconstructed = new double[n];
        for (int i = 0; i < half; i++)
        {
            double a = approximation[i];
            double d = detail[i];
            for (int k = 0; k < l; k++)
            {
                int idx = DiscreteWaveletTransform.Mod((2 * i) + k, n);
                reconstructed[idx] += (lo[k] * a) + (hi[k] * d);
            }
        }

        return reconstructed;
    }

    /// <summary>
    /// Reconstructs the original signal from a multi-level decomposition by combining the
    /// approximation with each detail band from coarsest to finest.
    /// </summary>
    public double[] Reconstruct(WaveletDecomposition decomposition)
    {
        if (decomposition == null) throw new ArgumentNullException(nameof(decomposition));

        double[] current = (double[])decomposition.Approximation.Clone();
        for (int level = decomposition.Levels - 1; level >= 0; level--)
            current = SingleLevel(current, decomposition.Details[level]);

        return current;
    }
}
