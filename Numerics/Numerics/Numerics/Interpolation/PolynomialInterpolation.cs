using System;
using System.Linq;

namespace CSharpNumerics.Numerics.Interpolation;

/// <summary>
/// Polynomial interpolation through N data points using Lagrange's formula
/// and Newton's divided-difference form.
/// 
/// <para>
/// Given N distinct points $(x_0 , y_0), \ldots , (x_{N-1}, y_{N-1})$,
/// there exists a unique polynomial $P(x)$ of degree ≤ N−1 passing through
/// all data points.
/// </para>
/// 
/// <list type="bullet">
///   <item><see cref="Lagrange"/> — Direct Lagrange evaluation $O(N^2)$.</item>
///   <item><see cref="Newton"/> — Newton divided-difference form $O(N^2)$ build, $O(N)$ evaluation.</item>
///   <item><see cref="Neville"/> — Neville's tableau $O(N^2)$, also returns error estimate.</item>
/// </list>
/// 
/// <example>
/// <code>
/// double[] x = { 0, 1, 2, 3 };
/// double[] y = { 1, 2, 0, 5 };
/// var poly = new PolynomialInterpolation(x, y);
/// double val = poly.Evaluate(1.5);          // Lagrange by default
/// double val2 = poly.EvaluateNewton(1.5);   // Newton form
/// var (val3, err) = poly.EvaluateNeville(1.5); // Neville with error estimate
/// </code>
/// </example>
/// </summary>
public class PolynomialInterpolation
{
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly double[] _dividedDiffs; // Newton coefficients (computed once)
    private readonly int _n;

    /// <summary>The x-coordinates of the data points.</summary>
    public ReadOnlySpan<double> X => _x;

    /// <summary>The y-coordinates of the data points.</summary>
    public ReadOnlySpan<double> Y => _y;

    /// <summary>Degree of the interpolating polynomial (N − 1).</summary>
    public int Degree => _n - 1;

    /// <summary>
    /// Constructs a polynomial interpolation from the given data points.
    /// </summary>
    /// <param name="x">Distinct x-coordinates (must not contain duplicates).</param>
    /// <param name="y">Corresponding y-values (same length as <paramref name="x"/>).</param>
    /// <exception cref="ArgumentException">Thrown if arrays differ in length or have fewer than 2 points.</exception>
    public PolynomialInterpolation(double[] x, double[] y)
    {
        if (x == null || y == null) throw new ArgumentNullException();
        if (x.Length != y.Length) throw new ArgumentException("x and y must have equal length.");
        if (x.Length < 2) throw new ArgumentException("At least 2 data points are required.");

        _n = x.Length;
        _x = (double[])x.Clone();
        _y = (double[])y.Clone();
        _dividedDiffs = BuildDividedDifferences(_x, _y);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Lagrange
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates the interpolating polynomial at <paramref name="xi"/> using
    /// the Lagrange basis form:
    /// $P(x) = \sum_{i=0}^{N-1} y_i \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$.
    /// </summary>
    public double Evaluate(double xi) => Lagrange(_x, _y, xi);

    /// <summary>Static Lagrange evaluation.</summary>
    public static double Lagrange(double[] x, double[] y, double xi)
    {
        int n = x.Length;
        double result = 0;
        for (int i = 0; i < n; i++)
        {
            double basis = 1;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                    basis *= (xi - x[j]) / (x[i] - x[j]);
            }
            result += y[i] * basis;
        }
        return result;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Newton divided-difference
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates the interpolating polynomial at <paramref name="xi"/> using
    /// Newton's divided-difference form:
    /// $P(x) = c_0 + c_1(x-x_0) + c_2(x-x_0)(x-x_1) + \ldots$
    /// </summary>
    public double EvaluateNewton(double xi)
    {
        double result = _dividedDiffs[_n - 1];
        for (int i = _n - 2; i >= 0; i--)
            result = result * (xi - _x[i]) + _dividedDiffs[i];
        return result;
    }

    /// <summary>Returns the Newton divided-difference coefficients.</summary>
    public double[] GetCoefficients() => (double[])_dividedDiffs.Clone();

    private static double[] BuildDividedDifferences(double[] x, double[] y)
    {
        int n = x.Length;
        var table = (double[])y.Clone();
        for (int j = 1; j < n; j++)
            for (int i = n - 1; i >= j; i--)
                table[i] = (table[i] - table[i - 1]) / (x[i] - x[i - j]);
        return table;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Neville (with error estimate)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates using Neville's algorithm. Returns the interpolated value
    /// together with an error estimate (difference of last two tableau entries).
    /// </summary>
    public (double value, double errorEstimate) EvaluateNeville(double xi)
    {
        var p = (double[])_y.Clone();
        double err = 0;
        for (int j = 1; j < _n; j++)
        {
            for (int i = 0; i < _n - j; i++)
            {
                double prev = p[i];
                p[i] = ((xi - _x[i + j]) * p[i] - (xi - _x[i]) * p[i + 1]) / (_x[i] - _x[i + j]);
                err = Math.Abs(p[i] - prev);
            }
        }
        return (p[0], err);
    }
}
