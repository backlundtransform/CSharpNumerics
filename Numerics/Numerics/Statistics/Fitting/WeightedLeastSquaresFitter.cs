using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Weighted least-squares fitting. Each observation is assigned a weight that
/// controls its influence on the fit. Solves (X^T W X) β = X^T W y where W = diag(weights).
/// </summary>
public static class WeightedLeastSquaresFitter
{
    /// <summary>
    /// Fits a weighted polynomial of the specified degree.
    /// </summary>
    public static FittingResult Fit(VectorN x, VectorN y, VectorN weights, int degree = 1)
    {
        if (x.Length != y.Length || x.Length != weights.Length)
            throw new ArgumentException("Vectors x, y, and weights must have the same length.");
        if (x.Length <= degree)
            throw new ArgumentException("Number of data points must exceed the polynomial degree.");
        if (degree < 0)
            throw new ArgumentException("Degree must be non-negative.");

        int n = x.Length;
        int p = degree + 1;
        double[,] X = FittingSolver.BuildPolynomialDesignMatrix(x, degree);
        return FitDesignMatrix(X, y.Values, weights.Values, n, p);
    }

    /// <summary>
    /// Fits a weighted multiple regression model using the supplied design matrix and weights.
    /// </summary>
    public static FittingResult Fit(double[,] designMatrix, VectorN y, VectorN weights)
    {
        int n = designMatrix.GetLength(0);
        int p = designMatrix.GetLength(1);
        if (n != y.Length || n != weights.Length)
            throw new ArgumentException("Design matrix rows, y, and weights must have the same length.");
        if (n <= p)
            throw new ArgumentException("Number of observations must exceed number of parameters.");

        return FitDesignMatrix(designMatrix, y.Values, weights.Values, n, p);
    }

    internal static FittingResult FitDesignMatrix(double[,] X, double[] y, double[] weights, int n, int p)
    {
        double[,] XtWX = FittingSolver.MultiplyATWA(X, weights, n, p);
        double[] XtWy = FittingSolver.MultiplyATWb(X, weights, y, n, p);
        double[] beta = FittingSolver.Solve(XtWX, XtWy);
        double[] fitted = FittingSolver.ComputeFitted(X, beta, n, p);

        double ssRes = 0;
        for (int i = 0; i < n; i++)
        {
            double r = y[i] - fitted[i];
            ssRes += weights[i] * r * r;
        }
        double s2 = (n > p) ? ssRes / (n - p) : 0.0;

        double[,] XtWXInv = FittingSolver.Invert(XtWX, p);
        double[] se = FittingSolver.ComputeStandardErrors(XtWXInv, s2, p);

        return new FittingResult(
            new VectorN(beta), new VectorN(y), new VectorN(fitted), new VectorN(se));
    }
}
