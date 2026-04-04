using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Robust regression using Iteratively Reweighted Least Squares (IRLS).
/// Down-weights outliers according to the chosen weight function so that
/// the fit is resistant to anomalous observations.
/// </summary>
public static class RobustFitter
{
    private const double DefaultHuberC = 1.345;
    private const double DefaultTukeyC = 4.685;
    private const double DefaultCauchyC = 2.385;
    private const double DefaultAndrewsC = 1.339;

    /// <summary>
    /// Robustly fits a polynomial of the specified degree.
    /// </summary>
    /// <param name="x">Independent variable values.</param>
    /// <param name="y">Dependent variable values.</param>
    /// <param name="degree">Polynomial degree (1 = linear).</param>
    /// <param name="weightFunction">Weight function for outlier down-weighting.</param>
    /// <param name="tuningConstant">Tuning constant; 0 selects the default for the chosen weight function.</param>
    /// <param name="maxIterations">Maximum IRLS iterations.</param>
    /// <param name="tolerance">Convergence tolerance on coefficient change.</param>
    public static FittingResult Fit(
        VectorN x,
        VectorN y,
        int degree = 1,
        RobustWeightFunction weightFunction = RobustWeightFunction.Huber,
        double tuningConstant = 0,
        int maxIterations = 50,
        double tolerance = 1e-6)
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
        return FitDesignMatrix(X, y.Values, n, p, weightFunction, tuningConstant, maxIterations, tolerance);
    }

    /// <summary>
    /// Robustly fits a model to the supplied design matrix and response.
    /// </summary>
    public static FittingResult Fit(
        double[,] designMatrix,
        VectorN y,
        RobustWeightFunction weightFunction = RobustWeightFunction.Huber,
        double tuningConstant = 0,
        int maxIterations = 50,
        double tolerance = 1e-6)
    {
        int n = designMatrix.GetLength(0);
        int p = designMatrix.GetLength(1);
        if (n != y.Length)
            throw new ArgumentException("Design matrix row count must match y length.");
        if (n <= p)
            throw new ArgumentException("Number of observations must exceed number of parameters.");

        return FitDesignMatrix(designMatrix, y.Values, n, p, weightFunction, tuningConstant, maxIterations, tolerance);
    }

    private static FittingResult FitDesignMatrix(
        double[,] X, double[] y, int n, int p,
        RobustWeightFunction wf, double c,
        int maxIterations, double tolerance)
    {
        if (c <= 0)
            c = DefaultC(wf);

        // Step 1: initial OLS fit
        double[,] XtX = FittingSolver.MultiplyATA(X, n, p);
        double[] Xty = FittingSolver.MultiplyATb(X, y, n, p);
        double[] beta = FittingSolver.Solve(XtX, Xty);

        double[] weights = new double[n];

        for (int iter = 0; iter < maxIterations; iter++)
        {
            double[] fitted = FittingSolver.ComputeFitted(X, beta, n, p);

            // Compute residuals
            double[] residuals = new double[n];
            for (int i = 0; i < n; i++)
                residuals[i] = y[i] - fitted[i];

            // Scale estimate via MAD
            double mad = FittingSolver.MAD(residuals);
            double scale = mad / 0.6745;
            if (scale < 1e-15)
                scale = 1.0; // All residuals nearly identical — use unit scale

            // Compute weights
            for (int i = 0; i < n; i++)
            {
                double u = residuals[i] / scale;
                weights[i] = ComputeWeight(u, c, wf);
            }

            // WLS step
            double[,] XtWX = FittingSolver.MultiplyATWA(X, weights, n, p);
            double[] XtWy = FittingSolver.MultiplyATWb(X, weights, y, n, p);

            double[] betaNew;
            try { betaNew = FittingSolver.Solve(XtWX, XtWy); }
            catch (InvalidOperationException) { break; }

            // Check convergence
            double maxDelta = 0;
            for (int j = 0; j < p; j++)
                maxDelta = Math.Max(maxDelta, Math.Abs(betaNew[j] - beta[j]));

            beta = betaNew;
            if (maxDelta < tolerance) break;
        }

        double[] finalFitted = FittingSolver.ComputeFitted(X, beta, n, p);

        // Standard errors: use final weighted residual variance
        double ssRes = 0;
        for (int i = 0; i < n; i++)
        {
            double r = y[i] - finalFitted[i];
            ssRes += weights[i] * r * r;
        }
        double s2 = (n > p) ? ssRes / (n - p) : 0.0;

        double[] se;
        try
        {
            double[,] XtWX = FittingSolver.MultiplyATWA(X, weights, n, p);
            double[,] XtWXInv = FittingSolver.Invert(XtWX, p);
            se = FittingSolver.ComputeStandardErrors(XtWXInv, s2, p);
        }
        catch (InvalidOperationException)
        {
            se = new double[p];
        }

        return new FittingResult(
            new VectorN(beta), new VectorN(y), new VectorN(finalFitted), new VectorN(se));
    }

    private static double DefaultC(RobustWeightFunction wf)
    {
        switch (wf)
        {
            case RobustWeightFunction.Huber: return DefaultHuberC;
            case RobustWeightFunction.TukeyBisquare: return DefaultTukeyC;
            case RobustWeightFunction.Cauchy: return DefaultCauchyC;
            case RobustWeightFunction.Andrews: return DefaultAndrewsC;
            default: return DefaultHuberC;
        }
    }

    private static double ComputeWeight(double u, double c, RobustWeightFunction wf)
    {
        double absU = Math.Abs(u);
        switch (wf)
        {
            case RobustWeightFunction.Huber:
                return absU <= c ? 1.0 : c / absU;

            case RobustWeightFunction.TukeyBisquare:
                if (absU > c) return 0.0;
                double t = 1.0 - (u / c) * (u / c);
                return t * t;

            case RobustWeightFunction.Cauchy:
                return 1.0 / (1.0 + (u / c) * (u / c));

            case RobustWeightFunction.Andrews:
                if (absU >= c * Math.PI) return 0.0;
                double v = u / c;
                return Math.Abs(v) < 1e-15 ? 1.0 : Math.Sin(v) / v;

            default:
                return 1.0;
        }
    }
}
