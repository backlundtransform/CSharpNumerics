using System;

namespace CSharpNumerics.Numerics.SignalProcessing;

/// <summary>
/// Savitzky-Golay smoothing filter that fits successive sub-sets of adjacent data points
/// with a low-degree polynomial by the method of linear least squares.
/// Preserves peak shape and higher moments better than moving-average filters.
/// </summary>
public class SavitzkyGolayFilter
{
    private readonly int _windowSize;
    private readonly int _polynomialOrder;
    private readonly double[] _coefficients;

    /// <summary>
    /// Creates a Savitzky-Golay filter with specified window size and polynomial order.
    /// </summary>
    /// <param name="windowSize">Number of points in the smoothing window (must be odd and ≥ 3).</param>
    /// <param name="polynomialOrder">Order of the fitting polynomial (must be less than windowSize).</param>
    public SavitzkyGolayFilter(int windowSize, int polynomialOrder)
    {
        if (windowSize < 3 || windowSize % 2 == 0)
            throw new ArgumentException("Window size must be odd and at least 3.", nameof(windowSize));
        if (polynomialOrder < 0 || polynomialOrder >= windowSize)
            throw new ArgumentException("Polynomial order must be non-negative and less than window size.", nameof(polynomialOrder));

        _windowSize = windowSize;
        _polynomialOrder = polynomialOrder;
        _coefficients = ComputeCoefficients(windowSize, polynomialOrder);
    }

    /// <summary>
    /// Gets the convolution coefficients for this filter configuration.
    /// </summary>
    public double[] Coefficients => (double[])_coefficients.Clone();

    /// <summary>
    /// Applies the Savitzky-Golay filter to the input signal.
    /// Boundary points are handled by reducing the window size symmetrically.
    /// </summary>
    /// <param name="signal">Input signal array.</param>
    /// <returns>Filtered signal of the same length.</returns>
    public double[] Apply(double[] signal)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (signal.Length < _windowSize)
            throw new ArgumentException($"Signal length must be at least {_windowSize}.", nameof(signal));

        int n = signal.Length;
        int halfWindow = _windowSize / 2;
        var result = new double[n];

        // Interior points: full convolution with precomputed coefficients
        for (int i = halfWindow; i < n - halfWindow; i++)
        {
            double sum = 0.0;
            for (int j = -halfWindow; j <= halfWindow; j++)
            {
                sum += _coefficients[j + halfWindow] * signal[i + j];
            }
            result[i] = sum;
        }

        // Boundary points: use reduced-size filter or mirror extension
        for (int i = 0; i < halfWindow; i++)
        {
            result[i] = ComputeBoundaryPoint(signal, i);
        }
        for (int i = n - halfWindow; i < n; i++)
        {
            result[i] = ComputeBoundaryPoint(signal, i);
        }

        return result;
    }

    /// <summary>
    /// Applies the filter and also computes the first derivative (smoothed).
    /// </summary>
    /// <param name="signal">Input signal array.</param>
    /// <param name="spacing">Spacing between sample points (e.g., dt).</param>
    /// <returns>Smoothed first derivative of the signal.</returns>
    public double[] ApplyDerivative(double[] signal, double spacing = 1.0)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (signal.Length < _windowSize)
            throw new ArgumentException($"Signal length must be at least {_windowSize}.", nameof(signal));

        int n = signal.Length;
        int halfWindow = _windowSize / 2;
        var derivCoeffs = ComputeDerivativeCoefficients(_windowSize, _polynomialOrder);
        var result = new double[n];

        for (int i = halfWindow; i < n - halfWindow; i++)
        {
            double sum = 0.0;
            for (int j = -halfWindow; j <= halfWindow; j++)
            {
                sum += derivCoeffs[j + halfWindow] * signal[i + j];
            }
            result[i] = sum / spacing;
        }

        // Boundary: forward/backward difference
        for (int i = 0; i < halfWindow && i < n - 1; i++)
        {
            result[i] = (signal[i + 1] - signal[i]) / spacing;
        }
        for (int i = Math.Max(halfWindow, n - halfWindow); i < n; i++)
        {
            result[i] = (signal[i] - signal[i - 1]) / spacing;
        }

        return result;
    }

    private double ComputeBoundaryPoint(double[] signal, int index)
    {
        // Use mirror extension for boundary handling
        int n = signal.Length;
        int halfWindow = _windowSize / 2;
        double sum = 0.0;

        for (int j = -halfWindow; j <= halfWindow; j++)
        {
            int idx = index + j;
            // Mirror reflection at boundaries
            if (idx < 0) idx = -idx;
            if (idx >= n) idx = 2 * (n - 1) - idx;
            idx = Math.Max(0, Math.Min(n - 1, idx));
            sum += _coefficients[j + halfWindow] * signal[idx];
        }

        return sum;
    }

    /// <summary>
    /// Computes the Savitzky-Golay convolution coefficients using the pseudo-inverse method.
    /// The coefficients are the zeroth row of (J^T J)^{-1} J^T where J is the Vandermonde matrix.
    /// </summary>
    private static double[] ComputeCoefficients(int windowSize, int polyOrder)
    {
        int halfWindow = windowSize / 2;
        int m = windowSize;
        int p = polyOrder + 1;

        // Build Vandermonde matrix J[i, k] = i^k where i ranges from -halfWindow to halfWindow
        var J = new double[m, p];
        for (int i = 0; i < m; i++)
        {
            double x = i - halfWindow;
            double xPow = 1.0;
            for (int k = 0; k < p; k++)
            {
                J[i, k] = xPow;
                xPow *= x;
            }
        }

        // Compute J^T * J (p×p)
        var JtJ = new double[p, p];
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < m; k++)
                    sum += J[k, i] * J[k, j];
                JtJ[i, j] = sum;
            }
        }

        // Invert J^T * J using Gauss-Jordan elimination
        var inv = InvertMatrix(JtJ, p);

        // Compute (J^T J)^{-1} J^T — we only need the zeroth row for smoothing
        var coeffs = new double[m];
        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;
            for (int k = 0; k < p; k++)
            {
                sum += inv[0, k] * J[i, k];
            }
            coeffs[i] = sum;
        }

        return coeffs;
    }

    /// <summary>
    /// Computes derivative coefficients (first row of (J^T J)^{-1} J^T).
    /// </summary>
    private static double[] ComputeDerivativeCoefficients(int windowSize, int polyOrder)
    {
        int halfWindow = windowSize / 2;
        int m = windowSize;
        int p = polyOrder + 1;

        if (p < 2) p = 2; // Need at least linear for derivative

        var J = new double[m, p];
        for (int i = 0; i < m; i++)
        {
            double x = i - halfWindow;
            double xPow = 1.0;
            for (int k = 0; k < p; k++)
            {
                J[i, k] = xPow;
                xPow *= x;
            }
        }

        var JtJ = new double[p, p];
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < m; k++)
                    sum += J[k, i] * J[k, j];
                JtJ[i, j] = sum;
            }
        }

        var inv = InvertMatrix(JtJ, p);

        // Row 1 (derivative coefficient row)
        var coeffs = new double[m];
        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;
            for (int k = 0; k < p; k++)
            {
                sum += inv[1, k] * J[i, k];
            }
            coeffs[i] = sum;
        }

        return coeffs;
    }

    private static double[,] InvertMatrix(double[,] matrix, int n)
    {
        var aug = new double[n, 2 * n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = matrix[i, j];
            aug[i, i + n] = 1.0;
        }

        for (int col = 0; col < n; col++)
        {
            // Partial pivoting
            int maxRow = col;
            double maxVal = Math.Abs(aug[col, col]);
            for (int row = col + 1; row < n; row++)
            {
                if (Math.Abs(aug[row, col]) > maxVal)
                {
                    maxVal = Math.Abs(aug[row, col]);
                    maxRow = row;
                }
            }

            if (maxRow != col)
            {
                for (int j = 0; j < 2 * n; j++)
                {
                    double tmp = aug[col, j];
                    aug[col, j] = aug[maxRow, j];
                    aug[maxRow, j] = tmp;
                }
            }

            double pivot = aug[col, col];
            if (Math.Abs(pivot) < 1e-15)
                throw new InvalidOperationException("Matrix is singular.");

            for (int j = 0; j < 2 * n; j++)
                aug[col, j] /= pivot;

            for (int row = 0; row < n; row++)
            {
                if (row == col) continue;
                double factor = aug[row, col];
                for (int j = 0; j < 2 * n; j++)
                    aug[row, j] -= factor * aug[col, j];
            }
        }

        var result = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i, j] = aug[i, j + n];

        return result;
    }
}
