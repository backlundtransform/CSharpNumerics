using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Fitting;
using System;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Time-series detrending methods: polynomial, moving average, median filter,
/// Savitzky–Golay smoothing, and differencing.
/// </summary>
public static class TimeSeriesDetrending
{
    /// <summary>
    /// Polynomial detrending: fits a polynomial of the given degree and subtracts it.
    /// Uses <see cref="LeastSquaresFitter"/> internally.
    /// </summary>
    /// <param name="times">Observation times (independent variable).</param>
    /// <param name="values">Observed values.</param>
    /// <param name="degree">Polynomial degree (1 = linear, 2 = quadratic, …).</param>
    /// <returns>Tuple of (detrended values, trend).</returns>
    public static (VectorN Detrended, VectorN Trend) PolynomialDetrend(
        VectorN times,
        VectorN values,
        int degree = 1)
    {
        if (times.Length != values.Length)
            throw new ArgumentException("times and values must have the same length.");
        if (degree < 0)
            throw new ArgumentException("degree must be >= 0.");

        FittingResult fit = LeastSquaresFitter.Fit(times, values, degree);
        VectorN detrended = values - fit.FittedValues;
        return (detrended, fit.FittedValues);
    }

    /// <summary>
    /// Moving-average detrending with a centred window.
    /// Edge points use a narrower (asymmetric) window.
    /// </summary>
    /// <param name="values">Observed values (assumed equally spaced).</param>
    /// <param name="windowSize">Window size (must be odd; if even it is incremented).</param>
    public static (VectorN Detrended, VectorN Trend) MovingAverageDetrend(
        VectorN values,
        int windowSize)
    {
        if (windowSize < 1)
            throw new ArgumentException("windowSize must be >= 1.");
        if (windowSize % 2 == 0) windowSize++;

        int n = values.Length;
        double[] y = values.Values;
        double[] trend = new double[n];
        int half = windowSize / 2;

        for (int i = 0; i < n; i++)
        {
            int lo = Math.Max(0, i - half);
            int hi = Math.Min(n - 1, i + half);
            double sum = 0;
            for (int j = lo; j <= hi; j++)
                sum += y[j];
            trend[i] = sum / (hi - lo + 1);
        }

        double[] detrended = new double[n];
        for (int i = 0; i < n; i++)
            detrended[i] = y[i] - trend[i];

        return (new VectorN(detrended), new VectorN(trend));
    }

    /// <summary>
    /// Median-filter detrending with a centred window. Robust to outliers.
    /// </summary>
    /// <param name="values">Observed values (assumed equally spaced).</param>
    /// <param name="windowSize">Window size (must be odd; if even it is incremented).</param>
    public static (VectorN Detrended, VectorN Trend) MedianFilterDetrend(
        VectorN values,
        int windowSize)
    {
        if (windowSize < 1)
            throw new ArgumentException("windowSize must be >= 1.");
        if (windowSize % 2 == 0) windowSize++;

        int n = values.Length;
        double[] y = values.Values;
        double[] trend = new double[n];
        int half = windowSize / 2;

        for (int i = 0; i < n; i++)
        {
            int lo = Math.Max(0, i - half);
            int hi = Math.Min(n - 1, i + half);
            int count = hi - lo + 1;
            double[] window = new double[count];
            for (int j = 0; j < count; j++)
                window[j] = y[lo + j];
            Array.Sort(window);
            trend[i] = count % 2 == 1
                ? window[count / 2]
                : (window[count / 2 - 1] + window[count / 2]) / 2.0;
        }

        double[] detrended = new double[n];
        for (int i = 0; i < n; i++)
            detrended[i] = y[i] - trend[i];

        return (new VectorN(detrended), new VectorN(trend));
    }

    /// <summary>
    /// Savitzky–Golay smoothing/detrending. Fits a local polynomial by least squares
    /// in a sliding window and uses the fitted value at each centre point as the trend.
    /// </summary>
    /// <param name="values">Observed values (assumed equally spaced).</param>
    /// <param name="windowSize">Window size (must be odd and > polyOrder).</param>
    /// <param name="polyOrder">Polynomial order for local fits (default 3).</param>
    public static (VectorN Detrended, VectorN Trend) SavitzkyGolayDetrend(
        VectorN values,
        int windowSize,
        int polyOrder = 3)
    {
        if (windowSize < 1)
            throw new ArgumentException("windowSize must be >= 1.");
        if (windowSize % 2 == 0) windowSize++;
        if (polyOrder < 0 || polyOrder >= windowSize)
            throw new ArgumentException("polyOrder must be >= 0 and < windowSize.");

        int n = values.Length;
        double[] y = values.Values;
        int half = windowSize / 2;

        // Precompute convolution coefficients for the zeroth row of (A^T A)^{-1} A^T
        // where A is the Vandermonde of [-half, ..., half].
        double[] coeffs = ComputeSGCoefficients(half, polyOrder);

        double[] trend = new double[n];

        for (int i = 0; i < n; i++)
        {
            if (i >= half && i < n - half)
            {
                // Full window: apply precomputed coefficients
                double val = 0;
                for (int j = -half; j <= half; j++)
                    val += coeffs[j + half] * y[i + j];
                trend[i] = val;
            }
            else
            {
                // Edge: fit local polynomial with available points
                int lo = Math.Max(0, i - half);
                int hi = Math.Min(n - 1, i + half);
                int count = hi - lo + 1;
                int effectiveOrder = Math.Min(polyOrder, count - 1);

                double center = i;
                double val = 0;
                double[] localCoeffs = ComputeSGCoefficientsLocal(lo, hi, i, effectiveOrder);
                for (int j = 0; j < count; j++)
                    val += localCoeffs[j] * y[lo + j];
                trend[i] = val;
            }
        }

        double[] detrended = new double[n];
        for (int i = 0; i < n; i++)
            detrended[i] = y[i] - trend[i];

        return (new VectorN(detrended), new VectorN(trend));
    }

    /// <summary>
    /// Differencing: Δ^d y[i] = y[i] − y[i−1] (applied d times).
    /// Output length is n − order.
    /// </summary>
    /// <param name="values">Observed values.</param>
    /// <param name="order">Differencing order (default 1).</param>
    public static VectorN Difference(VectorN values, int order = 1)
    {
        if (order < 1)
            throw new ArgumentException("order must be >= 1.");
        if (values.Length <= order)
            throw new ArgumentException("values length must be > order.");

        double[] current = (double[])values.Values.Clone();
        int len = current.Length;

        for (int d = 0; d < order; d++)
        {
            double[] next = new double[len - 1];
            for (int i = 0; i < len - 1; i++)
                next[i] = current[i + 1] - current[i];
            current = next;
            len--;
        }

        return new VectorN(current);
    }

    /// <summary>
    /// Compute SG convolution coefficients for a symmetric window of size 2*half+1.
    /// Solves (A^T A)^{-1} A^T for row 0 (smoothing, i.e. zeroth derivative).
    /// </summary>
    private static double[] ComputeSGCoefficients(int half, int polyOrder)
    {
        int winSize = 2 * half + 1;
        int p = polyOrder + 1;

        // Build Vandermonde matrix A (winSize × p)
        double[,] A = new double[winSize, p];
        for (int i = 0; i < winSize; i++)
        {
            double x = i - half;
            double xpow = 1.0;
            for (int j = 0; j < p; j++)
            {
                A[i, j] = xpow;
                xpow *= x;
            }
        }

        // Compute (A^T A)
        double[,] ATA = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double sum = 0;
                for (int k = 0; k < winSize; k++)
                    sum += A[k, i] * A[k, j];
                ATA[i, j] = sum;
            }

        // Invert ATA
        double[,] inv = InvertSmall(ATA, p);

        // Coefficients: row 0 of inv × A^T → smoothing weights
        double[] coeffs = new double[winSize];
        for (int i = 0; i < winSize; i++)
        {
            double sum = 0;
            for (int j = 0; j < p; j++)
                sum += inv[0, j] * A[i, j];
            coeffs[i] = sum;
        }

        return coeffs;
    }

    /// <summary>
    /// Compute SG coefficients for an asymmetric window (edge handling).
    /// Returns weights such that trend[center] = sum(weights[j] * y[lo+j]).
    /// </summary>
    private static double[] ComputeSGCoefficientsLocal(int lo, int hi, int center, int polyOrder)
    {
        int count = hi - lo + 1;
        int p = polyOrder + 1;

        double[,] A = new double[count, p];
        for (int i = 0; i < count; i++)
        {
            double x = (lo + i) - center;
            double xpow = 1.0;
            for (int j = 0; j < p; j++)
            {
                A[i, j] = xpow;
                xpow *= x;
            }
        }

        double[,] ATA = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double sum = 0;
                for (int k = 0; k < count; k++)
                    sum += A[k, i] * A[k, j];
                ATA[i, j] = sum;
            }

        double[,] inv = InvertSmall(ATA, p);

        double[] coeffs = new double[count];
        for (int i = 0; i < count; i++)
        {
            double sum = 0;
            for (int j = 0; j < p; j++)
                sum += inv[0, j] * A[i, j];
            coeffs[i] = sum;
        }

        return coeffs;
    }

    /// <summary>Gauss-Jordan inversion for small matrices.</summary>
    private static double[,] InvertSmall(double[,] A, int n)
    {
        double[,] aug = new double[n, 2 * n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = A[i, j];
            aug[i, i + n] = 1.0;
        }

        for (int col = 0; col < n; col++)
        {
            int pivot = col;
            double maxVal = Math.Abs(aug[col, col]);
            for (int row = col + 1; row < n; row++)
            {
                double v = Math.Abs(aug[row, col]);
                if (v > maxVal) { maxVal = v; pivot = row; }
            }

            if (pivot != col)
            {
                for (int j = 0; j < 2 * n; j++)
                {
                    double tmp = aug[col, j];
                    aug[col, j] = aug[pivot, j];
                    aug[pivot, j] = tmp;
                }
            }

            double diag = aug[col, col];
            if (Math.Abs(diag) < 1e-30) diag = 1e-30;

            for (int j = 0; j < 2 * n; j++)
                aug[col, j] /= diag;

            for (int row = 0; row < n; row++)
            {
                if (row == col) continue;
                double factor = aug[row, col];
                for (int j = 0; j < 2 * n; j++)
                    aug[row, j] -= factor * aug[col, j];
            }
        }

        double[,] result = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i, j] = aug[i, j + n];

        return result;
    }
}
