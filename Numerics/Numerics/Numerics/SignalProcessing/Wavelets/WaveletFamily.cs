using System;

namespace CSharpNumerics.Numerics.SignalProcessing.Wavelets;

/// <summary>
/// An orthonormal wavelet filter bank. Holds the decomposition (analysis) scaling/low-pass and
/// wavelet/high-pass filters. The high-pass filter is derived from the low-pass via the standard
/// quadrature-mirror relationship <c>g[k] = (-1)^k · h[L-1-k]</c>, which makes the periodic
/// transform orthonormal — so reconstruction is simply the transpose of analysis (exact to
/// machine precision).
/// </summary>
public class WaveletFamily
{
    private readonly double[] _decLow;
    private readonly double[] _decHigh;

    private WaveletFamily(string name, double[] decompositionLowPass)
    {
        if (decompositionLowPass == null || decompositionLowPass.Length < 2)
            throw new ArgumentException("A wavelet needs at least two scaling coefficients.", nameof(decompositionLowPass));

        Name = name;
        _decLow = (double[])decompositionLowPass.Clone();

        // Normalise the scaling filter to exact unit L2 norm to improve the orthonormality of the
        // periodic filter bank (and hence reconstruction accuracy with the transpose inverse).
        double norm = 0.0;
        foreach (double c in _decLow) norm += c * c;
        norm = Math.Sqrt(norm);
        for (int i = 0; i < _decLow.Length; i++) _decLow[i] /= norm;

        int l = _decLow.Length;
        _decHigh = new double[l];
        for (int k = 0; k < l; k++)
            _decHigh[k] = ((k % 2 == 0) ? 1.0 : -1.0) * _decLow[l - 1 - k];
    }

    /// <summary>Human-readable family name.</summary>
    public string Name { get; }

    /// <summary>Number of filter taps.</summary>
    public int Length => _decLow.Length;

    /// <summary>Decomposition low-pass (scaling) filter.</summary>
    public double[] DecompositionLowPass => (double[])_decLow.Clone();

    /// <summary>Decomposition high-pass (wavelet) filter.</summary>
    public double[] DecompositionHighPass => (double[])_decHigh.Clone();

    // Internal accessors avoid per-call array cloning in the transforms.
    internal double[] LowPass => _decLow;
    internal double[] HighPass => _decHigh;

    private const double Sqrt2 = 1.4142135623730951;

    /// <summary>Haar / Daubechies-1 wavelet (2 taps). The simplest orthonormal wavelet.</summary>
    public static WaveletFamily Haar { get; } = new("Haar", new[]
    {
        1.0 / Sqrt2, 1.0 / Sqrt2
    });

    /// <summary>Daubechies wavelet with 4 coefficients (db2). Computed from the closed form
    /// (1±√3)/(4√2), (3±√3)/(4√2) for full orthonormality to machine precision.</summary>
    public static WaveletFamily Daubechies4 { get; } = new("Daubechies4", BuildDaubechies4());

    private static double[] BuildDaubechies4()
    {
        double s3 = Math.Sqrt(3.0);
        double denom = 4.0 * Sqrt2;
        return new[]
        {
            (1.0 + s3) / denom,
            (3.0 + s3) / denom,
            (3.0 - s3) / denom,
            (1.0 - s3) / denom
        };
    }

    /// <summary>Daubechies wavelet with 8 coefficients (db4).</summary>
    public static WaveletFamily Daubechies8 { get; } = new("Daubechies8", new[]
    {
        0.23037781330885523,
        0.71484657055254153,
        0.63088076792959036,
        -0.02798376941698385,
        -0.18703481171888114,
        0.030841381835986965,
        0.032883011666982945,
        -0.010597401784997278
    });

    /// <summary>Symlet wavelet with 8 coefficients (sym4) — a near-symmetric Daubechies variant.</summary>
    public static WaveletFamily Symlet4 { get; } = new("Symlet4", new[]
    {
        -0.07576571478927333,
        -0.02963552764599851,
        0.49761866763201545,
        0.80373875180591614,
        0.29785779560527736,
        -0.09921954357684722,
        -0.012603967262037833,
        0.032223100604042702
    });
}
