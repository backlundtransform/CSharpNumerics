using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Result of a Box-fitting Least Squares (BLS) transit search.
/// </summary>
public class BLSResult
{
    /// <summary>Trial periods at which SR was evaluated.</summary>
    public VectorN Periods { get; }

    /// <summary>Signal Residue (SR) at each trial period.</summary>
    public VectorN Power { get; }

    /// <summary>Index of the best (highest SR) trial period.</summary>
    public int BestIndex { get; }

    /// <summary>Best-fit period.</summary>
    public double BestPeriod => Periods[BestIndex];

    /// <summary>SR at the best period.</summary>
    public double BestPower => Power[BestIndex];

    /// <summary>Best-fit transit epoch (mid-transit time).</summary>
    public double BestEpoch { get; }

    /// <summary>Best-fit transit duration as a fraction of the period.</summary>
    public double BestDurationFraction { get; }

    /// <summary>Best-fit transit depth (flux decrement).</summary>
    public double BestDepth { get; }

    public BLSResult(
        VectorN periods,
        VectorN power,
        int bestIndex,
        double bestEpoch,
        double bestDurationFraction,
        double bestDepth)
    {
        Periods = periods;
        Power = power;
        BestIndex = bestIndex;
        BestEpoch = bestEpoch;
        BestDurationFraction = bestDurationFraction;
        BestDepth = bestDepth;
    }
}
