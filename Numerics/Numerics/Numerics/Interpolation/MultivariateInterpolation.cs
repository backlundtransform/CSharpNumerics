using System;
using System.Linq;

namespace CSharpNumerics.Numerics.Interpolation;

/// <summary>
/// Multivariate interpolation for scattered data in $\mathbb{R}^d$.
/// 
/// Provides three methods:
/// <list type="bullet">
///   <item><b>Inverse Distance Weighting (IDW / Shepard)</b> — Simple, robust, no linear system to solve.</item>
///   <item><b>Radial Basis Function (RBF)</b> — Accurate for scattered data, supports multiple kernel functions.</item>
///   <item><b>Multilinear</b> — Fast tensor-product interpolation on a regular grid (2D/3D).</item>
/// </list>
/// 
/// <example>
/// <code>
/// // 2D scattered data: f(x,y) = x² + y²
/// double[][] points = { new[]{0.0,0.0}, new[]{1.0,0.0}, new[]{0.0,1.0}, new[]{1.0,1.0}, new[]{0.5,0.5} };
/// double[] values = points.Select(p => p[0]*p[0] + p[1]*p[1]).ToArray();
/// 
/// var interp = new MultivariateInterpolation(points, values);
/// double val = interp.EvaluateIDW(new[] { 0.25, 0.75 });             // IDW
/// double val2 = interp.EvaluateRBF(new[] { 0.25, 0.75 });            // RBF (Gaussian)
/// </code>
/// </example>
/// </summary>
public class MultivariateInterpolation
{
    private readonly double[][] _points;
    private readonly double[] _values;
    private readonly int _n;
    private readonly int _dim;

    // RBF cache
    private double[] _rbfWeights;
    private RbfKernel _cachedKernel;
    private double _cachedEpsilon;

    /// <summary>Number of data points.</summary>
    public int Count => _n;

    /// <summary>Dimensionality of input space.</summary>
    public int Dimension => _dim;

    /// <summary>
    /// Constructs a multivariate interpolation from scattered data.
    /// </summary>
    /// <param name="points">Array of N points, each a d-dimensional coordinate array.</param>
    /// <param name="values">Corresponding function values (length N).</param>
    public MultivariateInterpolation(double[][] points, double[] values)
    {
        if (points == null || values == null) throw new ArgumentNullException();
        if (points.Length != values.Length) throw new ArgumentException("points and values must have equal length.");
        if (points.Length < 2) throw new ArgumentException("At least 2 data points are required.");

        _n = points.Length;
        _dim = points[0].Length;
        _points = points.Select(p => (double[])p.Clone()).ToArray();
        _values = (double[])values.Clone();
    }

    // ═══════════════════════════════════════════════════════════════
    //  Inverse Distance Weighting (Shepard's method)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates at <paramref name="xi"/> using inverse distance weighting.
    /// $f(x) = \frac{\sum w_i \cdot f_i}{\sum w_i}$ where $w_i = 1 / d(x, x_i)^p$.
    /// </summary>
    /// <param name="xi">Query point (d-dimensional).</param>
    /// <param name="power">Distance exponent (default 2). Higher = more local influence.</param>
    public double EvaluateIDW(double[] xi, double power = 2.0)
    {
        if (xi.Length != _dim)
            throw new ArgumentException($"Query point must be {_dim}-dimensional.");

        double numer = 0, denom = 0;
        for (int i = 0; i < _n; i++)
        {
            double dist = Distance(xi, _points[i]);
            if (dist < 1e-15)
                return _values[i]; // Exact match

            double w = 1.0 / Math.Pow(dist, power);
            numer += w * _values[i];
            denom += w;
        }

        return numer / denom;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Radial Basis Function (RBF) interpolation
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Evaluates at <paramref name="xi"/> using Radial Basis Function interpolation.
    /// Solves the linear system $\Phi \mathbf{w} = \mathbf{f}$ on first call (cached).
    /// </summary>
    /// <param name="xi">Query point (d-dimensional).</param>
    /// <param name="kernel">RBF kernel function. Default <see cref="RbfKernel.Gaussian"/>.</param>
    /// <param name="epsilon">Shape parameter for the kernel. Default auto-estimated.</param>
    public double EvaluateRBF(double[] xi, RbfKernel kernel = RbfKernel.Gaussian, double epsilon = 0)
    {
        if (xi.Length != _dim)
            throw new ArgumentException($"Query point must be {_dim}-dimensional.");

        if (epsilon <= 0)
            epsilon = EstimateEpsilon();

        // Build weights if not cached or parameters changed
        if (_rbfWeights == null || _cachedKernel != kernel || Math.Abs(_cachedEpsilon - epsilon) > 1e-15)
        {
            _rbfWeights = SolveRBFWeights(kernel, epsilon);
            _cachedKernel = kernel;
            _cachedEpsilon = epsilon;
        }

        // Evaluate: f(x) = Σ w_i * φ(||x - x_i||)
        double result = 0;
        for (int i = 0; i < _n; i++)
        {
            double r = Distance(xi, _points[i]);
            result += _rbfWeights[i] * KernelFunction(r, epsilon, kernel);
        }
        return result;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Bilinear interpolation (regular 2D grid)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Bilinear interpolation on a regular 2D grid.
    /// </summary>
    /// <param name="xGrid">Sorted x-coordinates of the grid (length Nx).</param>
    /// <param name="yGrid">Sorted y-coordinates of the grid (length Ny).</param>
    /// <param name="gridValues">Grid values [ix, iy] of size Nx × Ny.</param>
    /// <param name="xi">Query x-coordinate.</param>
    /// <param name="yi">Query y-coordinate.</param>
    public static double Bilinear(double[] xGrid, double[] yGrid, double[,] gridValues, double xi, double yi)
    {
        int ix = FindGridIndex(xGrid, xi);
        int iy = FindGridIndex(yGrid, yi);

        double x0 = xGrid[ix], x1 = xGrid[ix + 1];
        double y0 = yGrid[iy], y1 = yGrid[iy + 1];

        double tx = (xi - x0) / (x1 - x0);
        double ty = (yi - y0) / (y1 - y0);

        double f00 = gridValues[ix, iy];
        double f10 = gridValues[ix + 1, iy];
        double f01 = gridValues[ix, iy + 1];
        double f11 = gridValues[ix + 1, iy + 1];

        return f00 * (1 - tx) * (1 - ty)
             + f10 * tx * (1 - ty)
             + f01 * (1 - tx) * ty
             + f11 * tx * ty;
    }

    /// <summary>
    /// Trilinear interpolation on a regular 3D grid.
    /// </summary>
    /// <param name="xGrid">Sorted x-coordinates.</param>
    /// <param name="yGrid">Sorted y-coordinates.</param>
    /// <param name="zGrid">Sorted z-coordinates.</param>
    /// <param name="gridValues">Grid values [ix, iy, iz] of size Nx × Ny × Nz.</param>
    /// <param name="xi">Query x.</param>
    /// <param name="yi">Query y.</param>
    /// <param name="zi">Query z.</param>
    public static double Trilinear(
        double[] xGrid, double[] yGrid, double[] zGrid,
        double[,,] gridValues,
        double xi, double yi, double zi)
    {
        int ix = FindGridIndex(xGrid, xi);
        int iy = FindGridIndex(yGrid, yi);
        int iz = FindGridIndex(zGrid, zi);

        double tx = (xi - xGrid[ix]) / (xGrid[ix + 1] - xGrid[ix]);
        double ty = (yi - yGrid[iy]) / (yGrid[iy + 1] - yGrid[iy]);
        double tz = (zi - zGrid[iz]) / (zGrid[iz + 1] - zGrid[iz]);

        // Interpolate along z at all 4 corners
        double c00 = gridValues[ix, iy, iz] * (1 - tz) + gridValues[ix, iy, iz + 1] * tz;
        double c10 = gridValues[ix + 1, iy, iz] * (1 - tz) + gridValues[ix + 1, iy, iz + 1] * tz;
        double c01 = gridValues[ix, iy + 1, iz] * (1 - tz) + gridValues[ix, iy + 1, iz + 1] * tz;
        double c11 = gridValues[ix + 1, iy + 1, iz] * (1 - tz) + gridValues[ix + 1, iy + 1, iz + 1] * tz;

        // Interpolate along y
        double c0 = c00 * (1 - ty) + c01 * ty;
        double c1 = c10 * (1 - ty) + c11 * ty;

        // Interpolate along x
        return c0 * (1 - tx) + c1 * tx;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static double Distance(double[] a, double[] b)
    {
        double sum = 0;
        for (int i = 0; i < a.Length; i++)
        {
            double d = a[i] - b[i];
            sum += d * d;
        }
        return Math.Sqrt(sum);
    }

    private double EstimateEpsilon()
    {
        // Use average nearest-neighbor distance as shape parameter
        double totalDist = 0;
        for (int i = 0; i < _n; i++)
        {
            double minDist = double.MaxValue;
            for (int j = 0; j < _n; j++)
            {
                if (i == j) continue;
                double d = Distance(_points[i], _points[j]);
                if (d < minDist) minDist = d;
            }
            totalDist += minDist;
        }
        return totalDist / _n;
    }

    private double[] SolveRBFWeights(RbfKernel kernel, double epsilon)
    {
        // Build Φ matrix
        var phi = new double[_n, _n];
        for (int i = 0; i < _n; i++)
        {
            for (int j = 0; j < _n; j++)
            {
                double r = Distance(_points[i], _points[j]);
                phi[i, j] = KernelFunction(r, epsilon, kernel);
            }
        }

        // Solve Φw = f via Gaussian elimination
        return SolveLinearSystem(phi, _values);
    }

    private static double KernelFunction(double r, double epsilon, RbfKernel kernel)
    {
        return kernel switch
        {
            RbfKernel.Gaussian => Math.Exp(-(r * r) / (epsilon * epsilon)),
            RbfKernel.Multiquadric => Math.Sqrt(1 + (r * r) / (epsilon * epsilon)),
            RbfKernel.InverseMultiquadric => 1.0 / Math.Sqrt(1 + (r * r) / (epsilon * epsilon)),
            RbfKernel.ThinPlateSpline => r > 0 ? r * r * Math.Log(r) : 0,
            RbfKernel.Cubic => r * r * r,
            _ => throw new ArgumentOutOfRangeException(nameof(kernel))
        };
    }

    private static double[] SolveLinearSystem(double[,] A, double[] b)
    {
        int n = b.Length;
        var aug = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = A[i, j];
            aug[i, n] = b[i];
        }

        // Gaussian elimination with partial pivoting
        for (int col = 0; col < n; col++)
        {
            int maxRow = col;
            for (int row = col + 1; row < n; row++)
                if (Math.Abs(aug[row, col]) > Math.Abs(aug[maxRow, col]))
                    maxRow = row;

            if (maxRow != col)
                for (int j = 0; j <= n; j++)
                    (aug[col, j], aug[maxRow, j]) = (aug[maxRow, j], aug[col, j]);

            double pivot = aug[col, col];
            if (Math.Abs(pivot) < 1e-14)
                throw new InvalidOperationException("Singular or near-singular RBF matrix. Try different epsilon or kernel.");

            for (int row = col + 1; row < n; row++)
            {
                double factor = aug[row, col] / pivot;
                for (int j = col; j <= n; j++)
                    aug[row, j] -= factor * aug[col, j];
            }
        }

        // Back substitution
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = aug[i, n];
            for (int j = i + 1; j < n; j++)
                sum -= aug[i, j] * x[j];
            x[i] = sum / aug[i, i];
        }
        return x;
    }

    private static int FindGridIndex(double[] grid, double val)
    {
        if (val <= grid[0]) return 0;
        if (val >= grid[grid.Length - 1]) return grid.Length - 2;

        int lo = 0, hi = grid.Length - 1;
        while (lo < hi - 1)
        {
            int mid = (lo + hi) / 2;
            if (val < grid[mid])
                hi = mid;
            else
                lo = mid;
        }
        return lo;
    }
}

/// <summary>Radial basis function kernel types.</summary>
public enum RbfKernel
{
    /// <summary>$\phi(r) = e^{-(r/\varepsilon)^2}$</summary>
    Gaussian,

    /// <summary>$\phi(r) = \sqrt{1 + (r/\varepsilon)^2}$</summary>
    Multiquadric,

    /// <summary>$\phi(r) = 1/\sqrt{1 + (r/\varepsilon)^2}$</summary>
    InverseMultiquadric,

    /// <summary>$\phi(r) = r^2 \ln(r)$</summary>
    ThinPlateSpline,

    /// <summary>$\phi(r) = r^3$</summary>
    Cubic
}
