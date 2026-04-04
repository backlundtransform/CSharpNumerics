using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Static methods for evaluating the quality of a fitted model.
/// </summary>
public static class GoodnessOfFit
{
    /// <summary>Coefficient of determination: R² = 1 − SS_res / SS_tot.</summary>
    public static double RSquared(VectorN observed, VectorN predicted)
    {
        int n = observed.Length;
        double yMean = 0;
        for (int i = 0; i < n; i++) yMean += observed[i];
        yMean /= n;

        double ssTot = 0, ssRes = 0;
        for (int i = 0; i < n; i++)
        {
            ssTot += (observed[i] - yMean) * (observed[i] - yMean);
            double r = observed[i] - predicted[i];
            ssRes += r * r;
        }

        return ssTot > 0 ? 1.0 - ssRes / ssTot : 0.0;
    }

    /// <summary>Adjusted R² penalising model complexity: 1 − (1−R²)(n−1)/(n−p).</summary>
    public static double AdjustedRSquared(double rSquared, int observationCount, int parameterCount)
    {
        int dof = observationCount - parameterCount;
        return dof > 0
            ? 1.0 - (1.0 - rSquared) * (observationCount - 1.0) / dof
            : rSquared;
    }

    /// <summary>Root mean squared error: √(Σ(y−ŷ)² / n).</summary>
    public static double RMSE(VectorN observed, VectorN predicted)
    {
        int n = observed.Length;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            double r = observed[i] - predicted[i];
            sum += r * r;
        }
        return n > 0 ? Math.Sqrt(sum / n) : 0.0;
    }

    /// <summary>Mean absolute error: Σ|y−ŷ| / n.</summary>
    public static double MAE(VectorN observed, VectorN predicted)
    {
        int n = observed.Length;
        double sum = 0;
        for (int i = 0; i < n; i++)
            sum += Math.Abs(observed[i] - predicted[i]);
        return n > 0 ? sum / n : 0.0;
    }

    /// <summary>Sum of squared residuals: Σ(y−ŷ)².</summary>
    public static double SSE(VectorN observed, VectorN predicted)
    {
        double sum = 0;
        for (int i = 0; i < observed.Length; i++)
        {
            double r = observed[i] - predicted[i];
            sum += r * r;
        }
        return sum;
    }

    /// <summary>Total sum of squares: Σ(y−ȳ)².</summary>
    public static double SST(VectorN observed)
    {
        int n = observed.Length;
        double yMean = 0;
        for (int i = 0; i < n; i++) yMean += observed[i];
        yMean /= n;

        double sum = 0;
        for (int i = 0; i < n; i++)
            sum += (observed[i] - yMean) * (observed[i] - yMean);
        return sum;
    }

    /// <summary>Mean squared error: SS_res / n.</summary>
    public static double MSE(VectorN observed, VectorN predicted)
    {
        int n = observed.Length;
        return n > 0 ? SSE(observed, predicted) / n : 0.0;
    }

    /// <summary>Akaike information criterion: n·ln(SS_res/n) + 2p.</summary>
    public static double AIC(VectorN observed, VectorN predicted, int parameterCount)
    {
        int n = observed.Length;
        double sse = SSE(observed, predicted);
        return n * Math.Log(sse / n) + 2.0 * parameterCount;
    }

    /// <summary>Bayesian information criterion: n·ln(SS_res/n) + p·ln(n).</summary>
    public static double BIC(VectorN observed, VectorN predicted, int parameterCount)
    {
        int n = observed.Length;
        double sse = SSE(observed, predicted);
        return n * Math.Log(sse / n) + parameterCount * Math.Log(n);
    }

    /// <summary>
    /// Performs residual analysis: computes standardized residuals, Durbin-Watson statistic,
    /// mean, standard deviation, skewness, and excess kurtosis of residuals.
    /// </summary>
    public static ResidualAnalysis AnalyzeResiduals(VectorN observed, VectorN predicted)
    {
        int n = observed.Length;
        double[] residuals = new double[n];
        double mean = 0;
        for (int i = 0; i < n; i++)
        {
            residuals[i] = observed[i] - predicted[i];
            mean += residuals[i];
        }
        mean /= n;

        double variance = 0;
        for (int i = 0; i < n; i++)
            variance += (residuals[i] - mean) * (residuals[i] - mean);
        variance /= (n > 1 ? n - 1 : 1);
        double stdDev = Math.Sqrt(variance);

        // Standardized residuals
        double[] std = new double[n];
        for (int i = 0; i < n; i++)
            std[i] = stdDev > 0 ? (residuals[i] - mean) / stdDev : 0.0;

        // Durbin-Watson
        double dwNum = 0, dwDen = 0;
        for (int i = 0; i < n; i++)
        {
            dwDen += residuals[i] * residuals[i];
            if (i > 0)
            {
                double diff = residuals[i] - residuals[i - 1];
                dwNum += diff * diff;
            }
        }
        double dw = dwDen > 0 ? dwNum / dwDen : 2.0;

        // Skewness and excess kurtosis
        double skew = 0, kurt = 0;
        if (stdDev > 0 && n > 2)
        {
            for (int i = 0; i < n; i++)
            {
                double z = (residuals[i] - mean) / stdDev;
                double z2 = z * z;
                skew += z * z2;
                kurt += z2 * z2;
            }
            skew /= n;
            kurt = kurt / n - 3.0;
        }

        return new ResidualAnalysis(new VectorN(std), dw, mean, stdDev, skew, kurt);
    }
}
