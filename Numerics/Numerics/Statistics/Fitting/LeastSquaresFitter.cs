using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Ordinary least-squares fitting for polynomials and general design matrices.
/// Solves the normal equations (X^T X) β = X^T y.
/// </summary>
public static class LeastSquaresFitter
{
    /// <summary>
    /// Fits a polynomial of the specified degree to (x, y) data.
    /// Degree 1 yields linear regression: y = a₀ + a₁x.
    /// Coefficients are ordered [a₀, a₁, …, a_d] where d is the degree.
    /// </summary>
    public static FittingResult Fit(VectorN x, VectorN y, int degree = 1)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("Vectors x and y must have the same length.");
        if (x.Length <= degree)
            throw new ArgumentException("Number of data points must exceed the polynomial degree.");
        if (degree < 0)
            throw new ArgumentException("Degree must be non-negative.");

        int n = x.Length;
        int p = degree + 1;
        double[,] X = FittingSolver.BuildPolynomialDesignMatrix(x, degree);
        return FitDesignMatrix(X, y.Values, n, p);
    }

    /// <summary>
    /// Fits a multiple linear regression model using the supplied design matrix.
    /// Include a column of ones in the design matrix if an intercept term is desired.
    /// </summary>
    public static FittingResult Fit(double[,] designMatrix, VectorN y)
    {
        int n = designMatrix.GetLength(0);
        int p = designMatrix.GetLength(1);
        if (n != y.Length)
            throw new ArgumentException("Design matrix row count must match y length.");
        if (n <= p)
            throw new ArgumentException("Number of observations must exceed number of parameters.");

        return FitDesignMatrix(designMatrix, y.Values, n, p);
    }

    private static FittingResult FitDesignMatrix(double[,] X, double[] y, int n, int p)
    {
        double[,] XtX = FittingSolver.MultiplyATA(X, n, p);
        double[] Xty = FittingSolver.MultiplyATb(X, y, n, p);
        double[] beta = FittingSolver.Solve(XtX, Xty);
        double[] fitted = FittingSolver.ComputeFitted(X, beta, n, p);

        double ssRes = 0;
        for (int i = 0; i < n; i++)
        {
            double r = y[i] - fitted[i];
            ssRes += r * r;
        }
        double s2 = (n > p) ? ssRes / (n - p) : 0.0;

        double[,] XtXInv = FittingSolver.Invert(XtX, p);
        double[] se = FittingSolver.ComputeStandardErrors(XtXInv, s2, p);

        return new FittingResult(
            new VectorN(beta), new VectorN(y), new VectorN(fitted), new VectorN(se));
    }
}
