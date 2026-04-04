namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Result of fitting a single peak (Gaussian or Lorentzian).
/// </summary>
public class PeakResult
{
    /// <summary>Index of the peak centre in the original data.</summary>
    public int PeakIndex { get; }

    /// <summary>Fitted centre position (x-axis).</summary>
    public double Center { get; }

    /// <summary>Fitted amplitude.</summary>
    public double Amplitude { get; }

    /// <summary>Fitted width parameter (σ for Gaussian, γ for Lorentzian).</summary>
    public double Width { get; }

    /// <summary>Fitted baseline offset.</summary>
    public double Baseline { get; }

    /// <summary>Shape that was fitted.</summary>
    public PeakShape Shape { get; }

    /// <summary>R² goodness-of-fit for the peak region.</summary>
    public double RSquared { get; }

    public PeakResult(int peakIndex, double center, double amplitude, double width,
        double baseline, PeakShape shape, double rSquared)
    {
        PeakIndex = peakIndex;
        Center = center;
        Amplitude = amplitude;
        Width = width;
        Baseline = baseline;
        Shape = shape;
        RSquared = rSquared;
    }
}
