namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Shape to use when fitting peaks.
/// </summary>
public enum PeakShape
{
    /// <summary>Gaussian: A · exp(−(x−μ)² / (2σ²)) + baseline</summary>
    Gaussian,

    /// <summary>Lorentzian (Cauchy): A · γ² / ((x−x₀)² + γ²) + baseline</summary>
    Lorentzian
}
