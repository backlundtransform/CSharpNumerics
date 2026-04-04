using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Parameter estimation utilities: confidence intervals for regression coefficients,
/// prediction intervals for new observations, maximum likelihood estimation,
/// and method-of-moments estimators for common distributions.
/// </summary>
public static class ParameterEstimation
{
    #region Confidence & Prediction Intervals

    /// <summary>
    /// Computes two-sided confidence intervals for each regression coefficient.
    /// CI = β_i ± t_{α/2, df} · SE_i.
    /// </summary>
    /// <param name="result">A fitting result containing coefficients, standard errors, and degrees of freedom.</param>
    /// <param name="confidenceLevel">Confidence level (e.g., 0.95 for 95 %).</param>
    public static (VectorN Lower, VectorN Upper) ConfidenceIntervals(
        FittingResult result,
        double confidenceLevel = 0.95)
    {
        double alpha = 1.0 - confidenceLevel;
        double tCrit = FittingSolver.TQuantile(1.0 - alpha / 2.0, result.DegreesOfFreedom);

        int p = result.ParameterCount;
        double[] lower = new double[p];
        double[] upper = new double[p];

        for (int i = 0; i < p; i++)
        {
            double margin = tCrit * result.StandardErrors[i];
            lower[i] = result.Coefficients[i] - margin;
            upper[i] = result.Coefficients[i] + margin;
        }

        return (new VectorN(lower), new VectorN(upper));
    }

    /// <summary>
    /// Computes a prediction interval for a new observation at a given design-matrix row.
    /// PI = ŷ ± t_{α/2, df} · s · √(1 + x'^T (X^T X)^{-1} x').
    /// </summary>
    /// <param name="result">A fitting result from a linear or polynomial fit.</param>
    /// <param name="designMatrix">The design matrix used during fitting (n × p).</param>
    /// <param name="xNewRow">The new observation's design row (length p).</param>
    /// <param name="confidenceLevel">Confidence level (e.g., 0.95).</param>
    public static (double Lower, double Upper) PredictionInterval(
        FittingResult result,
        double[,] designMatrix,
        VectorN xNewRow,
        double confidenceLevel = 0.95)
    {
        int n = designMatrix.GetLength(0);
        int p = designMatrix.GetLength(1);

        double[,] XtX = FittingSolver.MultiplyATA(designMatrix, n, p);
        double[,] XtXInv = FittingSolver.Invert(XtX, p);

        // Residual variance (unbiased)
        double ssRes = 0;
        for (int i = 0; i < result.Residuals.Length; i++)
            ssRes += result.Residuals[i] * result.Residuals[i];
        double s2 = result.DegreesOfFreedom > 0 ? ssRes / result.DegreesOfFreedom : 0.0;

        // x'^T (X^T X)^{-1} x'
        double leverage = 0;
        for (int i = 0; i < p; i++)
        {
            double sum = 0;
            for (int j = 0; j < p; j++)
                sum += XtXInv[i, j] * xNewRow[j];
            leverage += xNewRow[i] * sum;
        }

        double yHat = 0;
        for (int i = 0; i < p; i++)
            yHat += result.Coefficients[i] * xNewRow[i];

        double alphaCi = 1.0 - confidenceLevel;
        double tCrit = FittingSolver.TQuantile(1.0 - alphaCi / 2.0, result.DegreesOfFreedom);
        double margin = tCrit * Math.Sqrt(s2 * (1.0 + leverage));

        return (yHat - margin, yHat + margin);
    }

    #endregion

    #region Maximum Likelihood Estimation

    /// <summary>
    /// Finds the parameters θ that minimise a negative log-likelihood function
    /// using gradient descent with backtracking line search.
    /// </summary>
    /// <param name="negativeLogLikelihood">Function mapping parameter vector → −ln L(θ).</param>
    /// <param name="initialParameters">Starting point for the optimisation.</param>
    /// <param name="maxIterations">Maximum number of gradient descent steps.</param>
    /// <param name="tolerance">Convergence tolerance on gradient norm.</param>
    public static VectorN MaximumLikelihoodEstimate(
        Func<VectorN, double> negativeLogLikelihood,
        VectorN initialParameters,
        int maxIterations = 1000,
        double tolerance = 1e-8)
    {
        int p = initialParameters.Length;
        double[] theta = (double[])initialParameters.Values.Clone();
        const double h = 1e-7;

        for (int iter = 0; iter < maxIterations; iter++)
        {
            VectorN thetaVec = new VectorN(theta);
            double f0 = negativeLogLikelihood(thetaVec);

            // Numerical gradient (central differences)
            double[] grad = new double[p];
            for (int j = 0; j < p; j++)
            {
                double orig = theta[j];
                double step = Math.Max(h, Math.Abs(orig) * h);

                theta[j] = orig + step;
                double fPlus = negativeLogLikelihood(new VectorN(theta));
                theta[j] = orig - step;
                double fMinus = negativeLogLikelihood(new VectorN(theta));
                theta[j] = orig;

                grad[j] = (fPlus - fMinus) / (2.0 * step);
            }

            double gradNorm = 0;
            for (int j = 0; j < p; j++)
                gradNorm += grad[j] * grad[j];
            gradNorm = Math.Sqrt(gradNorm);

            if (gradNorm < tolerance) break;

            // Backtracking line search
            double alpha = 1.0;
            const double c1 = 1e-4;
            const double rho = 0.5;

            for (int ls = 0; ls < 40; ls++)
            {
                double[] thetaNew = new double[p];
                for (int j = 0; j < p; j++)
                    thetaNew[j] = theta[j] - alpha * grad[j];

                double fNew = negativeLogLikelihood(new VectorN(thetaNew));
                if (fNew <= f0 - c1 * alpha * gradNorm * gradNorm)
                {
                    theta = thetaNew;
                    break;
                }
                alpha *= rho;
            }
        }

        return new VectorN(theta);
    }

    #endregion

    #region Method of Moments

    /// <summary>
    /// Method-of-moments estimators for a normal distribution.
    /// Returns (μ̂ = sample mean, σ̂² = sample variance).
    /// </summary>
    public static (double Mean, double Variance) MethodOfMomentsNormal(VectorN data)
    {
        int n = data.Length;
        double mean = 0;
        for (int i = 0; i < n; i++) mean += data[i];
        mean /= n;

        double variance = 0;
        for (int i = 0; i < n; i++)
            variance += (data[i] - mean) * (data[i] - mean);
        variance /= (n > 1 ? n - 1 : 1);

        return (mean, variance);
    }

    /// <summary>
    /// Method-of-moments estimator for an exponential distribution.
    /// Returns the estimated rate λ̂ = 1 / x̄.
    /// </summary>
    public static double MethodOfMomentsExponential(VectorN data)
    {
        double mean = 0;
        for (int i = 0; i < data.Length; i++) mean += data[i];
        mean /= data.Length;
        return mean > 0 ? 1.0 / mean : 0.0;
    }

    /// <summary>
    /// Method-of-moments estimator for a Poisson distribution.
    /// Returns the estimated rate λ̂ = x̄.
    /// </summary>
    public static double MethodOfMomentsPoisson(VectorN data)
    {
        double mean = 0;
        for (int i = 0; i < data.Length; i++) mean += data[i];
        return mean / data.Length;
    }

    /// <summary>
    /// Method-of-moments estimators for a gamma distribution.
    /// Returns (α̂ = x̄²/s², β̂ = s²/x̄) where α is shape and β is scale.
    /// </summary>
    public static (double Shape, double Scale) MethodOfMomentsGamma(VectorN data)
    {
        int n = data.Length;
        double mean = 0;
        for (int i = 0; i < n; i++) mean += data[i];
        mean /= n;

        double variance = 0;
        for (int i = 0; i < n; i++)
            variance += (data[i] - mean) * (data[i] - mean);
        variance /= (n > 1 ? n - 1 : 1);

        if (variance <= 0 || mean <= 0) return (1.0, mean);
        double shape = mean * mean / variance;
        double scale = variance / mean;
        return (shape, scale);
    }

    #endregion
}
