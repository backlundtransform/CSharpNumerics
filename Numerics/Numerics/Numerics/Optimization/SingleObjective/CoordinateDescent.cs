using System;

namespace CSharpNumerics.Numerics.Optimization.SingleObjective;

/// <summary>
/// Coordinate descent optimiser for L1/L2 regularised problems.
/// Each iteration cycles through all coordinates and applies the
/// closed-form update with soft-thresholding for the L1 (Lasso) penalty.
/// Used by Lasso, ElasticNet, and similar models.
/// </summary>
public class CoordinateDescent
{
    public int MaxIterations { get; set; }
    public double Tolerance { get; set; }

    public CoordinateDescent(int maxIterations = 1000, double tolerance = 1e-7)
    {
        MaxIterations = maxIterations;
        Tolerance = tolerance;
    }

    /// <summary>
    /// Minimise ½‖Xw − y‖² + l1·‖w‖₁ + ½·l2·‖w‖²
    /// using cyclic coordinate descent with soft-thresholding.
    /// </summary>
    /// <param name="X">Design matrix (n × p), as a 2D array [row, col].</param>
    /// <param name="y">Target vector of length n.</param>
    /// <param name="l1">L1 (Lasso) penalty strength.</param>
    /// <param name="l2">L2 (Ridge) penalty strength.</param>
    /// <param name="skipBiasRegularisation">If true, coordinate 0 is not penalised (intercept).</param>
    /// <returns>Optimised weight vector of length p.</returns>
    public double[] Solve(double[,] X, double[] y, double l1, double l2,
        bool skipBiasRegularisation = false)
    {
        int n = X.GetLength(0);
        int p = X.GetLength(1);
        var w = new double[p];

        for (int iter = 0; iter < MaxIterations; iter++)
        {
            double maxChange = 0;

            for (int j = 0; j < p; j++)
            {
                double rho = 0;
                double norm = 0;

                for (int i = 0; i < n; i++)
                {
                    double pred = 0;
                    for (int k = 0; k < p; k++)
                        if (k != j)
                            pred += X[i, k] * w[k];

                    rho += X[i, j] * (y[i] - pred);
                    norm += X[i, j] * X[i, j];
                }

                double oldW = w[j];

                if (skipBiasRegularisation && j == 0)
                {
                    w[j] = norm > 0 ? rho / norm : 0;
                }
                else
                {
                    double numerator = SoftThreshold(rho / norm, l1);
                    double denominator = 1.0 + l2;
                    w[j] = norm > 0 ? numerator / denominator : 0;
                }

                double change = Math.Abs(w[j] - oldW);
                if (change > maxChange)
                    maxChange = change;
            }

            if (maxChange < Tolerance)
                break;
        }

        return w;
    }

    private static double SoftThreshold(double z, double gamma)
        => Math.Sign(z) * Math.Max(Math.Abs(z) - gamma, 0.0);
}
