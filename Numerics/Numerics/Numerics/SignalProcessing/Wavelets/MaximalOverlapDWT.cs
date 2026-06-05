using System;

namespace CSharpNumerics.Numerics.SignalProcessing.Wavelets;

/// <summary>
/// Result of a MODWT decomposition: per-level detail coefficients (<see cref="Details"/>[0] =
/// level 1) and the final smooth band. Every band has the same length as the input.
/// </summary>
public class ModwtDecomposition
{
    public ModwtDecomposition(double[][] details, double[] smooth)
    {
        Details = details ?? throw new ArgumentNullException(nameof(details));
        Smooth = smooth ?? throw new ArgumentNullException(nameof(smooth));
    }

    /// <summary>Detail coefficients per level; index 0 = level 1 (finest).</summary>
    public double[][] Details { get; }

    /// <summary>Final smooth (scaling) band.</summary>
    public double[] Smooth { get; }

    /// <summary>Number of levels.</summary>
    public int Levels => Details.Length;
}

/// <summary>
/// Maximal Overlap Discrete Wavelet Transform (MODWT), an undecimated, shift-invariant wavelet
/// transform. Unlike the critically-sampled <see cref="DiscreteWaveletTransform"/>, no
/// downsampling occurs: every band keeps the full signal length and a circular shift of the input
/// produces the same circular shift of every coefficient band — ideal for time-series alignment.
/// Filters are rescaled by 1/√2 and dilated (à trous) by 2^(level-1) at each level. Works for any
/// signal length (periodic boundary).
/// </summary>
public class MaximalOverlapDWT
{
    private readonly double[] _scalingLow;   // scaling filter / √2
    private readonly double[] _scalingHigh;  // wavelet filter / √2

    public MaximalOverlapDWT(WaveletFamily family)
    {
        Family = family ?? throw new ArgumentNullException(nameof(family));

        double[] lo = family.LowPass;
        double[] hi = family.HighPass;
        int l = lo.Length;

        _scalingLow = new double[l];
        _scalingHigh = new double[l];
        const double invSqrt2 = 0.7071067811865476;
        for (int k = 0; k < l; k++)
        {
            _scalingLow[k] = lo[k] * invSqrt2;
            _scalingHigh[k] = hi[k] * invSqrt2;
        }
    }

    public WaveletFamily Family { get; }

    /// <summary>Forward MODWT to <paramref name="levels"/> levels.</summary>
    public ModwtDecomposition Forward(double[] signal, int levels)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (levels < 1) throw new ArgumentOutOfRangeException(nameof(levels), "Levels must be at least 1.");

        int n = signal.Length;
        int l = _scalingLow.Length;
        var details = new double[levels][];
        double[] v = (double[])signal.Clone();   // running smooth band (V_0 = signal)

        for (int level = 0; level < levels; level++)
        {
            int spacing = 1 << level;            // 2^level dilation
            var smooth = new double[n];
            var detail = new double[n];

            for (int t = 0; t < n; t++)
            {
                double s = 0.0, d = 0.0;
                for (int k = 0; k < l; k++)
                {
                    int idx = DiscreteWaveletTransform.Mod(t - (spacing * k), n);
                    s += _scalingLow[k] * v[idx];
                    d += _scalingHigh[k] * v[idx];
                }
                smooth[t] = s;
                detail[t] = d;
            }

            details[level] = detail;
            v = smooth;
        }

        return new ModwtDecomposition(details, v);
    }

    /// <summary>Inverse MODWT — reconstructs the signal (exact to machine precision).</summary>
    public double[] Inverse(ModwtDecomposition decomposition)
    {
        if (decomposition == null) throw new ArgumentNullException(nameof(decomposition));

        int n = decomposition.Smooth.Length;
        int l = _scalingLow.Length;
        double[] v = (double[])decomposition.Smooth.Clone();

        for (int level = decomposition.Levels - 1; level >= 0; level--)
        {
            int spacing = 1 << level;
            double[] detail = decomposition.Details[level];
            var previous = new double[n];

            // Adjoint of the forward step (note the +spacing*k shift).
            for (int t = 0; t < n; t++)
            {
                double value = 0.0;
                for (int k = 0; k < l; k++)
                {
                    int idx = DiscreteWaveletTransform.Mod(t + (spacing * k), n);
                    value += (_scalingLow[k] * v[idx]) + (_scalingHigh[k] * detail[idx]);
                }
                previous[t] = value;
            }

            v = previous;
        }

        return v;
    }
}
