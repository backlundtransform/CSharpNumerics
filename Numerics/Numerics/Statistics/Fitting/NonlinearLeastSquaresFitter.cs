using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Nonlinear least-squares fitting using the Levenberg-Marquardt algorithm.
/// The Jacobian is computed numerically via central finite differences.
/// </summary>
public static class NonlinearLeastSquaresFitter
{
    /// <summary>
    /// Fits a nonlinear model y = f(x; θ) to observed (x, y) data.
    /// </summary>
    /// <param name="model">Model function: f(x, θ) → predicted y.</param>
    /// <param name="x">Independent variable values.</param>
    /// <param name="y">Observed dependent variable values.</param>
    /// <param name="initialParameters">Initial parameter guess.</param>
    /// <param name="maxIterations">Maximum number of LM iterations.</param>
    /// <param name="tolerance">Convergence tolerance on parameter change.</param>
    public static FittingResult Fit(
        Func<double, VectorN, double> model,
        VectorN x,
        VectorN y,
        VectorN initialParameters,
        int maxIterations = 200,
        double tolerance = 1e-10)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("Vectors x and y must have the same length.");
        if (initialParameters.Length == 0)
            throw new ArgumentException("Initial parameters must be provided.");

        int n = x.Length;
        int p = initialParameters.Length;

        double[] xArr = x.Values;
        double[] yArr = y.Values;
        double[] theta = (double[])initialParameters.Values.Clone();
        double lambda = 1e-3;
        const double lambdaUp = 10.0;
        const double lambdaDown = 0.1;

        double cost = ComputeCost(model, xArr, yArr, theta, n);

        for (int iter = 0; iter < maxIterations; iter++)
        {
            double[] residuals = ComputeResiduals(model, xArr, yArr, theta, n);
            double[,] J = ComputeJacobian(model, xArr, theta, n, p);

            double[,] JtJ = FittingSolver.MultiplyATA(J, n, p);
            double[] Jtr = FittingSolver.MultiplyATb(J, residuals, n, p);

            // Marquardt damping: J^T J + λ · diag(J^T J)
            for (int i = 0; i < p; i++)
                JtJ[i, i] += lambda * (JtJ[i, i] > 0 ? JtJ[i, i] : 1.0);

            double[] delta;
            try { delta = FittingSolver.Solve(JtJ, Jtr); }
            catch (InvalidOperationException) { break; }

            double[] thetaNew = new double[p];
            for (int i = 0; i < p; i++)
                thetaNew[i] = theta[i] + delta[i];

            double costNew = ComputeCost(model, xArr, yArr, thetaNew, n);

            if (costNew < cost)
            {
                theta = thetaNew;
                cost = costNew;
                lambda *= lambdaDown;
            }
            else
            {
                lambda *= lambdaUp;
            }

            double deltaMax = 0;
            for (int i = 0; i < p; i++)
                deltaMax = Math.Max(deltaMax, Math.Abs(delta[i]));
            if (deltaMax < tolerance) break;
        }

        double[] fitted = new double[n];
        VectorN thetaVec = new VectorN(theta);
        for (int i = 0; i < n; i++)
            fitted[i] = model(xArr[i], thetaVec);

        // Standard errors via Jacobian at solution
        double[] se = ComputeStandardErrorsAtSolution(model, xArr, yArr, theta, fitted, n, p);

        return new FittingResult(
            thetaVec, y, new VectorN(fitted), new VectorN(se));
    }

    private static double[] ComputeResiduals(
        Func<double, VectorN, double> model, double[] x, double[] y, double[] theta, int n)
    {
        VectorN thetaVec = new VectorN(theta);
        double[] r = new double[n];
        for (int i = 0; i < n; i++)
            r[i] = y[i] - model(x[i], thetaVec);
        return r;
    }

    private static double ComputeCost(
        Func<double, VectorN, double> model, double[] x, double[] y, double[] theta, int n)
    {
        VectorN thetaVec = new VectorN(theta);
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            double r = y[i] - model(x[i], thetaVec);
            sum += r * r;
        }
        return sum;
    }

    private static double[,] ComputeJacobian(
        Func<double, VectorN, double> model, double[] x, double[] theta, int n, int p)
    {
        double[,] J = new double[n, p];
        double[] thetaPerturbed = (double[])theta.Clone();
        const double h = 1e-7;

        for (int j = 0; j < p; j++)
        {
            double original = theta[j];
            double step = Math.Max(h, Math.Abs(original) * h);

            thetaPerturbed[j] = original + step;
            VectorN pPlus = new VectorN(thetaPerturbed);
            double[] fPlus = new double[n];
            for (int i = 0; i < n; i++)
                fPlus[i] = model(x[i], pPlus);

            thetaPerturbed[j] = original - step;
            VectorN pMinus = new VectorN(thetaPerturbed);
            double[] fMinus = new double[n];
            for (int i = 0; i < n; i++)
                fMinus[i] = model(x[i], pMinus);

            thetaPerturbed[j] = original;
            double denom = 2.0 * step;
            for (int i = 0; i < n; i++)
                J[i, j] = (fPlus[i] - fMinus[i]) / denom;
        }

        return J;
    }

    private static double[] ComputeStandardErrorsAtSolution(
        Func<double, VectorN, double> model, double[] x, double[] y,
        double[] theta, double[] fitted, int n, int p)
    {
        try
        {
            double[,] J = ComputeJacobian(model, x, theta, n, p);
            double ssRes = 0;
            for (int i = 0; i < n; i++)
            {
                double r = y[i] - fitted[i];
                ssRes += r * r;
            }
            double s2 = (n > p) ? ssRes / (n - p) : 0.0;

            double[,] JtJ = FittingSolver.MultiplyATA(J, n, p);
            double[,] cov = FittingSolver.Invert(JtJ, p);
            return FittingSolver.ComputeStandardErrors(cov, s2, p);
        }
        catch (InvalidOperationException)
        {
            return new double[p];
        }
    }
}
