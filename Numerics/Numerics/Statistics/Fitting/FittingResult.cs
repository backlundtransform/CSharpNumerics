using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Encapsulates the results of a curve-fitting operation, including coefficients,
/// residuals, goodness-of-fit metrics, and standard errors.
/// </summary>
public class FittingResult
{
    /// <summary>Estimated model coefficients (e.g. [intercept, β₁, β₂, …] for regression).</summary>
    public VectorN Coefficients { get; }

    /// <summary>Residuals: observed[i] − fitted[i].</summary>
    public VectorN Residuals { get; }

    /// <summary>Fitted (predicted) values.</summary>
    public VectorN FittedValues { get; }

    /// <summary>Coefficient of determination (R²).</summary>
    public double RSquared { get; }

    /// <summary>Adjusted R² accounting for the number of predictors.</summary>
    public double AdjustedRSquared { get; }

    /// <summary>Root mean squared error.</summary>
    public double RMSE { get; }

    /// <summary>Mean absolute error.</summary>
    public double MAE { get; }

    /// <summary>Residual degrees of freedom (n − p).</summary>
    public int DegreesOfFreedom { get; }

    /// <summary>Standard errors of the estimated coefficients.</summary>
    public VectorN StandardErrors { get; }

    /// <summary>Number of estimated parameters.</summary>
    public int ParameterCount { get; }

    /// <summary>Number of observations used in the fit.</summary>
    public int ObservationCount { get; }

    public FittingResult(
        VectorN coefficients,
        VectorN observedValues,
        VectorN fittedValues,
        VectorN standardErrors)
    {
        int n = observedValues.Length;
        int p = coefficients.Length;

        Coefficients = coefficients;
        FittedValues = fittedValues;
        StandardErrors = standardErrors;
        ObservationCount = n;
        ParameterCount = p;
        DegreesOfFreedom = n - p;

        double[] residuals = new double[n];
        double ssRes = 0;
        double yMean = 0;
        double absResSum = 0;

        for (int i = 0; i < n; i++)
            yMean += observedValues[i];
        yMean /= n;

        double ssTot = 0;
        for (int i = 0; i < n; i++)
        {
            double r = observedValues[i] - fittedValues[i];
            residuals[i] = r;
            ssRes += r * r;
            ssTot += (observedValues[i] - yMean) * (observedValues[i] - yMean);
            absResSum += Math.Abs(r);
        }

        Residuals = new VectorN(residuals);

        RSquared = ssTot > 0 ? 1.0 - ssRes / ssTot : 0.0;
        int dof = n - p;
        AdjustedRSquared = (dof > 0 && ssTot > 0)
            ? 1.0 - (1.0 - RSquared) * (n - 1.0) / dof
            : RSquared;
        RMSE = n > 0 ? Math.Sqrt(ssRes / n) : 0.0;
        MAE = n > 0 ? absResSum / n : 0.0;
    }
}
