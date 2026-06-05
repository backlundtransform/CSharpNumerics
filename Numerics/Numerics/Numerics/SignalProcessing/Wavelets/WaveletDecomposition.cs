using System;

namespace CSharpNumerics.Numerics.SignalProcessing.Wavelets;

/// <summary>
/// Result of an N-level discrete wavelet decomposition: the coarsest approximation plus the
/// detail coefficients at each level. <see cref="Details"/>[0] holds the finest (level-1) detail
/// and <see cref="Details"/>[Levels-1] the coarsest.
/// </summary>
public class WaveletDecomposition
{
    public WaveletDecomposition(double[] approximation, double[][] details, WaveletFamily family)
    {
        Approximation = approximation ?? throw new ArgumentNullException(nameof(approximation));
        Details = details ?? throw new ArgumentNullException(nameof(details));
        Family = family ?? throw new ArgumentNullException(nameof(family));
    }

    /// <summary>Coarsest-scale approximation (scaling) coefficients.</summary>
    public double[] Approximation { get; }

    /// <summary>Detail (wavelet) coefficients per level; index 0 = finest.</summary>
    public double[][] Details { get; }

    /// <summary>Number of decomposition levels.</summary>
    public int Levels => Details.Length;

    /// <summary>The wavelet family used for the decomposition.</summary>
    public WaveletFamily Family { get; }
}
