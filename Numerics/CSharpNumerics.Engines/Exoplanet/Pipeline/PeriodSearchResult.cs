using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public class PeriodSearchResult
{
    public double BestPeriod { get; }
    public double BestPower { get; }
    public double BestEpoch { get; }
    public double BestDepth { get; }
    public double BestDurationFraction { get; }
    public double FalseAlarmProbability { get; }
    public VectorN Periods { get; }
    public VectorN Power { get; }

    public PeriodSearchResult(double bestPeriod, double bestPower, double bestEpoch,
        double bestDepth, double bestDurationFraction, double falseAlarmProbability,
        VectorN periods, VectorN power)
    {
        BestPeriod = bestPeriod;
        BestPower = bestPower;
        BestEpoch = bestEpoch;
        BestDepth = bestDepth;
        BestDurationFraction = bestDurationFraction;
        FalseAlarmProbability = falseAlarmProbability;
        Periods = periods;
        Power = power;
    }
}
