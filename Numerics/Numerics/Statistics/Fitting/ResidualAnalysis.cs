using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Contains the results of a residual analysis, including standardized residuals,
/// autocorrelation (Durbin-Watson), and distributional shape measures.
/// </summary>
public class ResidualAnalysis
{
    /// <summary>Residuals divided by their standard deviation.</summary>
    public VectorN StandardizedResiduals { get; }

    /// <summary>Durbin-Watson statistic for first-order autocorrelation (≈ 2 means no autocorrelation).</summary>
    public double DurbinWatsonStatistic { get; }

    /// <summary>Mean of the residuals (should be near 0 for a well-specified model).</summary>
    public double MeanResidual { get; }

    /// <summary>Standard deviation of the residuals.</summary>
    public double ResidualStandardDeviation { get; }

    /// <summary>Skewness of the residual distribution (0 = symmetric).</summary>
    public double Skewness { get; }

    /// <summary>Excess kurtosis of the residual distribution (0 = normal).</summary>
    public double Kurtosis { get; }

    public ResidualAnalysis(
        VectorN standardizedResiduals,
        double durbinWatsonStatistic,
        double meanResidual,
        double residualStandardDeviation,
        double skewness,
        double kurtosis)
    {
        StandardizedResiduals = standardizedResiduals;
        DurbinWatsonStatistic = durbinWatsonStatistic;
        MeanResidual = meanResidual;
        ResidualStandardDeviation = residualStandardDeviation;
        Skewness = skewness;
        Kurtosis = kurtosis;
    }
}
