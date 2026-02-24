using System;
using System.Linq;

namespace CSharpNumerics.Numerics.Interpolation;

/// <summary>
/// Cubic spline interpolation through N data points.
/// 
/// Supports three boundary conditions:
/// <list type="bullet">
///   <item><see cref="SplineBoundary.Natural"/> — Second derivative = 0 at endpoints (free ends).</item>
///   <item><see cref="SplineBoundary.Clamped"/> — First derivative specified at endpoints.</item>
///   <item><see cref="SplineBoundary.NotAKnot"/> — Third derivative continuous at second and second-to-last knots.</item>
/// </list>
/// 
/// <example>
/// <code>
/// double[] x = { 0, 1, 2, 3, 4 };
/// double[] y = { 0, 1, 0, 1, 0 };
/// var spline = new CubicSplineInterpolation(x, y);
/// double val = spline.Evaluate(2.5);       // smooth interpolation
/// double dydx = spline.Derivative(2.5);    // first derivative 
/// double d2ydx2 = spline.SecondDerivative(2.5);  // curvature
/// </code>
/// </example>
/// </summary>
public class CubicSplineInterpolation
{
    private readonly double[] _x;
    private readonly double[] _y;
    private readonly double[] _a, _b, _c, _d; // coefficients per interval
    private readonly int _n; // number of data points

    /// <summary>The x-coordinates (sorted, distinct).</summary>
    public ReadOnlySpan<double> X => _x;

    /// <summary>The y-coordinates.</summary>
    public ReadOnlySpan<double> Y => _y;

    /// <summary>Number of data points.</summary>
    public int Count => _n;

    /// <summary>
    /// Constructs a natural cubic spline (second derivative = 0 at endpoints).
    /// </summary>
    public CubicSplineInterpolation(double[] x, double[] y)
        : this(x, y, SplineBoundary.Natural)
    {
    }

    /// <summary>
    /// Constructs a cubic spline with the specified boundary condition.
    /// </summary>
    /// <param name="x">Sorted, distinct x-coordinates.</param>
    /// <param name="y">Corresponding y-values.</param>
    /// <param name="boundary">Boundary condition type.</param>
    /// <param name="clampedLeft">Left endpoint first derivative (only for <see cref="SplineBoundary.Clamped"/>).</param>
    /// <param name="clampedRight">Right endpoint first derivative (only for <see cref="SplineBoundary.Clamped"/>).</param>
    public CubicSplineInterpolation(
        double[] x, double[] y,
        SplineBoundary boundary,
        double clampedLeft = 0, double clampedRight = 0)
    {
        if (x == null || y == null) throw new ArgumentNullException();
        if (x.Length != y.Length) throw new ArgumentException("x and y must have equal length.");
        if (x.Length < 3) throw new ArgumentException("At least 3 data points required for cubic spline.");

        _n = x.Length;
        _x = (double[])x.Clone();
        _y = (double[])y.Clone();

        int m = _n - 1; // number of intervals
        _a = new double[m];
        _b = new double[m];
        _c = new double[m];
        _d = new double[m];

        ComputeCoefficients(boundary, clampedLeft, clampedRight);
    }

    // ═══════════════════════════════════════════════════════════════

    /// <summary>Evaluates the spline at <paramref name="xi"/>.</summary>
    public double Evaluate(double xi)
    {
        int i = FindInterval(xi);
        double dx = xi - _x[i];
        return _a[i] + _b[i] * dx + _c[i] * dx * dx + _d[i] * dx * dx * dx;
    }

    /// <summary>Evaluates the first derivative at <paramref name="xi"/>.</summary>
    public double Derivative(double xi)
    {
        int i = FindInterval(xi);
        double dx = xi - _x[i];
        return _b[i] + 2.0 * _c[i] * dx + 3.0 * _d[i] * dx * dx;
    }

    /// <summary>Evaluates the second derivative at <paramref name="xi"/>.</summary>
    public double SecondDerivative(double xi)
    {
        int i = FindInterval(xi);
        double dx = xi - _x[i];
        return 2.0 * _c[i] + 6.0 * _d[i] * dx;
    }

    /// <summary>
    /// Evaluates the spline at many points efficiently.
    /// </summary>
    public double[] Evaluate(double[] xs)
    {
        var result = new double[xs.Length];
        for (int k = 0; k < xs.Length; k++)
            result[k] = Evaluate(xs[k]);
        return result;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Coefficient computation
    // ═══════════════════════════════════════════════════════════════

    private void ComputeCoefficients(SplineBoundary boundary, double clampedLeft, double clampedRight)
    {
        int m = _n - 1;
        var h = new double[m];
        for (int i = 0; i < m; i++)
            h[i] = _x[i + 1] - _x[i];

        // Solve for the second derivatives (M[i])
        double[] M = SolveSecondDerivatives(h, boundary, clampedLeft, clampedRight);

        // Compute cubic coefficients for each interval [x_i, x_{i+1}]
        for (int i = 0; i < m; i++)
        {
            _a[i] = _y[i];
            _b[i] = (_y[i + 1] - _y[i]) / h[i] - h[i] * (2.0 * M[i] + M[i + 1]) / 6.0;
            _c[i] = M[i] / 2.0;
            _d[i] = (M[i + 1] - M[i]) / (6.0 * h[i]);
        }
    }

    private double[] SolveSecondDerivatives(double[] h, SplineBoundary boundary, double clampedLeft, double clampedRight)
    {
        int m = _n - 1;
        var M = new double[_n];

        // Build tridiagonal system
        var lower = new double[_n];
        var diag = new double[_n];
        var upper = new double[_n];
        var rhs = new double[_n];

        // Interior equations (same for all boundary types)
        for (int i = 1; i < m; i++)
        {
            lower[i] = h[i - 1];
            diag[i] = 2.0 * (h[i - 1] + h[i]);
            upper[i] = h[i];
            rhs[i] = 6.0 * ((_y[i + 1] - _y[i]) / h[i] - (_y[i] - _y[i - 1]) / h[i - 1]);
        }

        switch (boundary)
        {
            case SplineBoundary.Natural:
                diag[0] = 1; upper[0] = 0; rhs[0] = 0;
                diag[m] = 1; lower[m] = 0; rhs[m] = 0;
                break;

            case SplineBoundary.Clamped:
                diag[0] = 2.0 * h[0];
                upper[0] = h[0];
                rhs[0] = 6.0 * ((_y[1] - _y[0]) / h[0] - clampedLeft);

                diag[m] = 2.0 * h[m - 1];
                lower[m] = h[m - 1];
                rhs[m] = 6.0 * (clampedRight - (_y[m] - _y[m - 1]) / h[m - 1]);
                break;

            case SplineBoundary.NotAKnot:
                // First equation: continuity of third derivative at x[1]
                diag[0] = h[1];
                upper[0] = -(h[0] + h[1]);
                rhs[0] = 0;
                // We need M[0], M[1], M[2] relation: h[1]*M[0] - (h[0]+h[1])*M[1] + h[0]*M[2] = 0
                // Place M[2] contribution via modifying upper row
                // Use row 0: h[1]*M[0] - (h[0]+h[1])*M[1] + h[0]*M[2] = 0
                // This is a 3-wide row: [h[1], -(h[0]+h[1]), h[0]] at positions [0,1,2]
                // For tridiagonal solver we handle the extra via reduction
                // Actually, let's use the full approach:
                // We'll fold the not-a-knot condition into the first and last rows
                // Row 0: M[0] * h[1] + M[1] * (-(h[0] + h[1])) + M[2] * h[0] = 0
                // Since upper[0] connects to M[1], we need to handle M[2] separately
                // Use a simple substitution approach instead:
                ComputeNotAKnot(h, M);
                return M;
        }

        // Solve tridiagonal system (Thomas algorithm)
        SolveTridiagonal(lower, diag, upper, rhs, M, _n);
        return M;
    }

    private void ComputeNotAKnot(double[] h, double[] M)
    {
        int m = _n - 1;

        // Build full system of n equations
        // Interior: h[i-1]*M[i-1] + 2(h[i-1]+h[i])*M[i] + h[i]*M[i+1] = rhs[i]
        // First row (not-a-knot at x[1]):    h[1]*M[0] - (h[0]+h[1])*M[1] + h[0]*M[2] = 0
        // Last row  (not-a-knot at x[m-1]):  h[m-1]*M[m-2] - (h[m-2]+h[m-1])*M[m-1] + h[m-2]*M[m] = 0

        // For n >= 4, this works directly. For n == 3, reduce to a simpler system.

        if (_n == 3)
        {
            // Not-a-knot with 3 points → M[0] = M[1] = M[2] = second derivative of the unique quadratic
            // d³S is constant → M[0]=M[1] and M[1]=M[2] → all equal
            // From interior equation: h[0]*M[0] + 2(h[0]+h[1])*M[1] + h[1]*M[2] = rhs
            double rhs1 = 6.0 * ((_y[2] - _y[1]) / h[1] - (_y[1] - _y[0]) / h[0]);
            double mVal = rhs1 / (h[0] + 2.0 * (h[0] + h[1]) + h[1]);
            M[0] = mVal; M[1] = mVal; M[2] = mVal;
            return;
        }

        // General case: n >= 4 — Gaussian elimination on dense system
        int n = _n;
        var A = new double[n, n];
        var b = new double[n];

        // Row 0: not-a-knot at x[1]
        A[0, 0] = h[1]; A[0, 1] = -(h[0] + h[1]); A[0, 2] = h[0];
        b[0] = 0;

        // Interior rows
        for (int i = 1; i < m; i++)
        {
            A[i, i - 1] = h[i - 1];
            A[i, i] = 2.0 * (h[i - 1] + h[i]);
            A[i, i + 1] = h[i];
            b[i] = 6.0 * ((_y[i + 1] - _y[i]) / h[i] - (_y[i] - _y[i - 1]) / h[i - 1]);
        }

        // Row m: not-a-knot at x[m-1]
        A[m, m - 2] = h[m - 1]; A[m, m - 1] = -(h[m - 2] + h[m - 1]); A[m, m] = h[m - 2];
        b[m] = 0;

        // Gaussian elimination with partial pivoting
        for (int col = 0; col < n; col++)
        {
            // Pivoting
            int maxRow = col;
            for (int row = col + 1; row < n; row++)
                if (Math.Abs(A[row, col]) > Math.Abs(A[maxRow, col]))
                    maxRow = row;
            if (maxRow != col)
            {
                for (int j = 0; j < n; j++)
                    (A[col, j], A[maxRow, j]) = (A[maxRow, j], A[col, j]);
                (b[col], b[maxRow]) = (b[maxRow], b[col]);
            }

            for (int row = col + 1; row < n; row++)
            {
                double factor = A[row, col] / A[col, col];
                for (int j = col; j < n; j++)
                    A[row, j] -= factor * A[col, j];
                b[row] -= factor * b[col];
            }
        }

        // Back substitution
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = b[i];
            for (int j = i + 1; j < n; j++)
                sum -= A[i, j] * M[j];
            M[i] = sum / A[i, i];
        }
    }

    private static void SolveTridiagonal(double[] lower, double[] diag, double[] upper, double[] rhs, double[] result, int n)
    {
        // Thomas algorithm — forward sweep
        var cp = new double[n];
        var dp = new double[n];

        cp[0] = upper[0] / diag[0];
        dp[0] = rhs[0] / diag[0];

        for (int i = 1; i < n; i++)
        {
            double m = 1.0 / (diag[i] - lower[i] * cp[i - 1]);
            cp[i] = upper[i] * m;
            dp[i] = (rhs[i] - lower[i] * dp[i - 1]) * m;
        }

        // Back substitution
        result[n - 1] = dp[n - 1];
        for (int i = n - 2; i >= 0; i--)
            result[i] = dp[i] - cp[i] * result[i + 1];
    }

    private int FindInterval(double xi)
    {
        int m = _n - 1;
        // Clamp to valid range
        if (xi <= _x[0]) return 0;
        if (xi >= _x[m]) return m - 1;

        // Binary search
        int lo = 0, hi = m;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (xi < _x[mid])
                hi = mid;
            else
                lo = mid + 1;
        }
        return lo - 1;
    }
}

/// <summary>Boundary condition type for cubic spline interpolation.</summary>
public enum SplineBoundary
{
    /// <summary>Second derivative = 0 at both endpoints (natural/free spline).</summary>
    Natural,

    /// <summary>First derivative specified at both endpoints.</summary>
    Clamped,

    /// <summary>Third derivative continuous at x[1] and x[N−2] (not-a-knot).</summary>
    NotAKnot
}
