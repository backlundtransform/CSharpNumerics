namespace CSharpNumerics.Statistics.Fitting;

/// <summary>
/// Weight function used by the robust fitter (IRLS) to down-weight outlying residuals.
/// </summary>
public enum RobustWeightFunction
{
    /// <summary>Huber weight: w = 1 for |u| ≤ c, w = c/|u| otherwise. Default c = 1.345.</summary>
    Huber,

    /// <summary>Tukey bisquare: w = (1 − (u/c)²)² for |u| ≤ c, 0 otherwise. Default c = 4.685.</summary>
    TukeyBisquare,

    /// <summary>Cauchy/Lorentzian: w = 1 / (1 + (u/c)²). Default c = 2.385.</summary>
    Cauchy,

    /// <summary>Andrews sine wave: w = sin(u/c)/(u/c) for |u| &lt; cπ, 0 otherwise. Default c = 1.339.</summary>
    Andrews
}
