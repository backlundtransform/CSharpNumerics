using CSharpNumerics.Engines.Exoplanet.Enums;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public class TransitDetectionConfig
{
    public double MinPeriodDays { get; set; } = 0.5;
    public double MaxPeriodDays { get; set; } = 100.0;
    public double MinTransitDepthPpm { get; set; } = 100.0;
    public double SnrThreshold { get; set; } = 7.0;
    public int MaxPlanets { get; set; } = 3;
    public DetrendingMethod DetrendingMethod { get; set; } = DetrendingMethod.MedianFilter;
    public PeriodSearchMethod PeriodSearchMethod { get; set; } = PeriodSearchMethod.BLS;
    public int DetrendingWindowSize { get; set; } = 101;
    public int DetrendingPolyDegree { get; set; } = 3;
    public double OutlierSigmaThreshold { get; set; } = 5.0;
    public int NumTrialPeriods { get; set; } = 5000;
    public int PhaseBins { get; set; } = 200;
}
