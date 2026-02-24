using System;

namespace CSharpNumerics.Numerics.Interpolation;

/// <summary>
/// Rational interpolation using Bulirsch–Stoer's algorithm.
/// 
/// <para>
/// Builds a rational function $R(x) = P(x)/Q(x)$ (ratio of two polynomials)
/// that passes through all N data points. Rational interpolation is superior
/// to polynomial interpolation when the function has poles or near-singularities.
/// </para>
/// 
/// <para>
/// Also provides <see cref="FloaterHormann"/> — a barycentric rational interpolant
/// that is guaranteed pole-free and requires no pole detection.
/// </para>
/// 
/// <example>
/// <code>
/// double[] x = { 0, 1, 2, 3, 4 };
/// double[] y = { 1.0, 0.5, 0.333, 0.25, 0.2 }; // ≈ 1/(1+x)
/// var interp = new RationalInterpolation(x, y);
/// double val = interp.Evaluate(1.5);   // Bulirsch–Stoer
/// double val2 = interp.EvaluateFloaterHormann(1.5, d: 3); // barycentric, blending order d
/// </code>
/// </example>
/// </summary>
public class RationalInterpolation
{
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly int _n;

    /// <summary>The x-coordinates of the data points.</summary>
    public ReadOnlySpan<double> X => _x;

    /// <summary>The y-coordinates of the data points.</summary>
    public ReadOnlySpan<double> Y => _y;

    /// <summary>Number of data points.</summary>
    public int Count => _n;

    /// <summary>
    /// Constructs a rational interpolation from the given data points.
    /// </summary>
    public RationalInterpolation(double[] x, double[] y)
    {
        if (x == null || y == null) throw new ArgumentNullException();
        if (x.Length != y.Length) throw new ArgumentException("x and y must have equal length.");
        if (x.Length < 2) throw new ArgumentException("At least 2 data points are required.");

        _n = x.Length;
        _x = (double[])x.Clone();
        _y = (double[])y.Clone();
    }

    // ═══════════════════════════════════════════════════════════════
    //  Bulirsch–Stoer algorithm
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates the rational interpolant at <paramref name="xi"/> using the
    /// Bulirsch–Stoer algorithm. Returns the interpolated value and an error estimate.
    /// </summary>
    public (double value, double errorEstimate) Evaluate(double xi)
    {
        var c = new double[_n];
        var d = new double[_n];
        double tiny = 1e-25;
        int ns = 0;
        double hh = Math.Abs(xi - _x[0]);

        for (int i = 0; i < _n; i++)
        {
            double h = Math.Abs(xi - _x[i]);
            if (h == 0.0)
                return (_y[i], 0.0);
            if (h < hh)
            {
                ns = i;
                hh = h;
            }
            c[i] = _y[i];
            d[i] = _y[i] + tiny; // tiny to prevent rare 0/0
        }

        double result = _y[ns--];
        double err = 0;

        for (int m = 1; m < _n; m++)
        {
            for (int i = 0; i < _n - m; i++)
            {
                double w = c[i + 1] - d[i];
                double h = _x[i + m] - xi; // h is never zero (checked above)
                double t = (_x[i] - xi) * d[i] / h;
                double dd = t - c[i + 1];

                if (dd == 0.0)
                    dd = tiny; // Pole encountered — use tiny offset

                dd = w / dd;
                d[i] = c[i + 1] * dd;
                c[i] = t * dd;
            }

            err = (2 * (ns + 1) < (_n - m)) ? c[ns + 1] : d[ns--];
            result += err;
        }

        return (result, Math.Abs(err));
    }

    // ═══════════════════════════════════════════════════════════════
    //  Floater–Hormann barycentric rational interpolation
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates the Floater–Hormann barycentric rational interpolant at <paramref name="xi"/>.
    /// 
    /// <para>
    /// The blending parameter <paramref name="d"/> controls smoothness:
    /// d=0 → piecewise constant, d=1 → piecewise linear, d=N−1 → polynomial.
    /// Typical choice: d = 3 (good balance).
    /// </para>
    /// 
    /// <para>
    /// This interpolant is guaranteed to have no poles in the convex hull of the data.
    /// </para>
    /// </summary>
    /// <param name="xi">The x-coordinate at which to evaluate.</param>
    /// <param name="d">Blending order (0 ≤ d ≤ N−1). Default 3.</param>
    public double EvaluateFloaterHormann(double xi, int d = 3)
    {
        if (d < 0) d = 0;
        if (d >= _n) d = _n - 1;

        // Compute barycentric weights
        double[] w = ComputeFloaterHormannWeights(d);

        // Barycentric formula: R(x) = Σ w_i / (x - x_i) * y_i  /  Σ w_i / (x - x_i)
        double numer = 0, denom = 0;
        for (int i = 0; i < _n; i++)
        {
            double diff = xi - _x[i];
            if (Math.Abs(diff) < 1e-15)
                return _y[i];

            double t = w[i] / diff;
            numer += t * _y[i];
            denom += t;
        }

        return numer / denom;
    }

    private double[] ComputeFloaterHormannWeights(int d)
    {
        var w = new double[_n];

        for (int i = 0; i < _n; i++)
        {
            int jMin = Math.Max(0, i - d);
            int jMax = Math.Min(i, _n - d - 1);
            double sum = 0;

            for (int j = jMin; j <= jMax; j++)
            {
                double prod = 1.0;
                for (int k = j; k <= j + d; k++)
                {
                    if (k != i)
                        prod *= Math.Abs(_x[i] - _x[k]);
                }
                sum += 1.0 / prod;
            }

            w[i] = ((i - d) % 2 == 0 ? 1 : -1) * sum;
            // More robust sign: (-1)^(i-d) but only when d is well-defined
            // Actually the sign should be (-1)^i for the barycentric formula
            // Floater-Hormann: w_i = (-1)^i * sum_{j in J_i} 1/prod
            w[i] = (i % 2 == 0 ? 1 : -1) * sum;
        }

        return w;
    }
}
