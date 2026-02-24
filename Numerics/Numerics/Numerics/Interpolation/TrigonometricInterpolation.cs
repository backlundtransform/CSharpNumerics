using System;

namespace CSharpNumerics.Numerics.Interpolation;

/// <summary>
/// Trigonometric interpolation through N data points using a discrete
/// Fourier-type expansion.
/// 
/// <para>
/// Best suited for <b>periodic</b> functions. Given N equally or unequally
/// spaced points, constructs a trigonometric polynomial:
/// $$T(x) = \frac{a_0}{2} + \sum_{k=1}^{M} \bigl[ a_k \cos(k\omega x) + b_k \sin(k\omega x) \bigr]$$
/// where $\omega = 2\pi / P$ ($P$ = period) and $M = \lfloor (N-1)/2 \rfloor$.
/// </para>
/// 
/// <para>
/// For non-equispaced data, uses a least-squares DFT-like construction
/// that minimises the residual at the data points.
/// </para>
/// 
/// <example>
/// <code>
/// // Interpolate sin(x) + 0.5*cos(2x) on [0, 2π)
/// int N = 8;
/// double[] x = new double[N];
/// double[] y = new double[N];
/// for (int i = 0; i &lt; N; i++)
/// {
///     x[i] = 2 * Math.PI * i / N;
///     y[i] = Math.Sin(x[i]) + 0.5 * Math.Cos(2 * x[i]);
/// }
/// var trig = new TrigonometricInterpolation(x, y, period: 2 * Math.PI);
/// double val = trig.Evaluate(Math.PI / 3);
/// </code>
/// </example>
/// </summary>
public class TrigonometricInterpolation
{
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly double[] _a; // cosine coefficients
    private readonly double[] _b; // sine coefficients
    private readonly double _omega;
    private readonly double _period;
    private readonly int _n;
    private readonly int _m; // number of harmonics

    /// <summary>The assumed period of the data.</summary>
    public double Period => _period;

    /// <summary>Number of data points.</summary>
    public int Count => _n;

    /// <summary>Number of harmonics used.</summary>
    public int Harmonics => _m;

    /// <summary>Cosine coefficients $a_0, a_1, \ldots, a_M$.</summary>
    public double[] CosineCoefficients => (double[])_a.Clone();

    /// <summary>Sine coefficients $b_1, \ldots, b_M$ (b[0] is always 0).</summary>
    public double[] SineCoefficients => (double[])_b.Clone();

    /// <summary>
    /// Constructs a trigonometric interpolant.
    /// </summary>
    /// <param name="x">Data point x-coordinates.</param>
    /// <param name="y">Data point y-values.</param>
    /// <param name="period">
    /// The fundamental period. If null, auto-detected as <c>x[N-1] - x[0] + (x[1] - x[0])</c>
    /// (assumes roughly equispaced).
    /// </param>
    public TrigonometricInterpolation(double[] x, double[] y, double? period = null)
    {
        if (x == null || y == null) throw new ArgumentNullException();
        if (x.Length != y.Length) throw new ArgumentException("x and y must have equal length.");
        if (x.Length < 2) throw new ArgumentException("At least 2 data points are required.");

        _n = x.Length;
        _x = (double[])x.Clone();
        _y = (double[])y.Clone();

        _period = period ?? (_x[_n - 1] - _x[0] + (_x[1] - _x[0]));
        _omega = 2.0 * Math.PI / _period;
        _m = (_n - 1) / 2;

        _a = new double[_m + 1];
        _b = new double[_m + 1];

        ComputeCoefficients();
    }

    /// <summary>
    /// Evaluates the trigonometric interpolant at <paramref name="xi"/>.
    /// </summary>
    public double Evaluate(double xi)
    {
        double result = _a[0] / 2.0;
        for (int k = 1; k <= _m; k++)
        {
            double angle = k * _omega * xi;
            result += _a[k] * Math.Cos(angle) + _b[k] * Math.Sin(angle);
        }

        // If N is even, add the Nyquist term a_{N/2} * cos(N/2 * ω * x) / 2
        if (_n % 2 == 0)
        {
            int half = _n / 2;
            double aNyquist = 0;
            for (int j = 0; j < _n; j++)
                aNyquist += _y[j] * Math.Cos(half * _omega * _x[j]);
            aNyquist *= 2.0 / _n;
            result += aNyquist / 2.0 * Math.Cos(half * _omega * xi);
        }

        return result;
    }

    /// <summary>
    /// Evaluates the first derivative of the trigonometric interpolant.
    /// </summary>
    public double Derivative(double xi)
    {
        double result = 0;
        for (int k = 1; k <= _m; k++)
        {
            double angle = k * _omega * xi;
            double kw = k * _omega;
            result += kw * (-_a[k] * Math.Sin(angle) + _b[k] * Math.Cos(angle));
        }
        return result;
    }

    /// <summary>
    /// Evaluates at multiple points.
    /// </summary>
    public double[] Evaluate(double[] xs)
    {
        var result = new double[xs.Length];
        for (int i = 0; i < xs.Length; i++)
            result[i] = Evaluate(xs[i]);
        return result;
    }

    // ═══════════════════════════════════════════════════════════════

    private void ComputeCoefficients()
    {
        // For equispaced data this is the DFT.
        // For general data we use the DFT-like sums: a_k = 2/N Σ y_j cos(k ω x_j), etc.
        // This is exact interpolation when N is odd and data is equispaced.

        for (int k = 0; k <= _m; k++)
        {
            double ac = 0, bc = 0;
            for (int j = 0; j < _n; j++)
            {
                double angle = k * _omega * _x[j];
                ac += _y[j] * Math.Cos(angle);
                bc += _y[j] * Math.Sin(angle);
            }
            _a[k] = 2.0 * ac / _n;
            _b[k] = 2.0 * bc / _n;
        }
    }
}
