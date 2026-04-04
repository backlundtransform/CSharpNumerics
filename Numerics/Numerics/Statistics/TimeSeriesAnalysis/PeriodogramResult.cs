using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Result of a periodogram computation (Lomb-Scargle or similar).
/// </summary>
public class PeriodogramResult
{
    /// <summary>Trial frequencies at which power was evaluated.</summary>
    public VectorN Frequencies { get; }

    /// <summary>Spectral power at each trial frequency.</summary>
    public VectorN Power { get; }

    /// <summary>Trial periods (1 / frequency) for convenience.</summary>
    public VectorN Periods { get; }

    /// <summary>Index of the peak power value.</summary>
    public int BestIndex { get; }

    /// <summary>Frequency of the strongest peak.</summary>
    public double BestFrequency => Frequencies[BestIndex];

    /// <summary>Period of the strongest peak.</summary>
    public double BestPeriod => Periods[BestIndex];

    /// <summary>Power at the strongest peak.</summary>
    public double BestPower => Power[BestIndex];

    public PeriodogramResult(VectorN frequencies, VectorN power)
    {
        Frequencies = frequencies;
        Power = power;

        double[] periods = new double[frequencies.Length];
        int bestIdx = 0;
        double bestPow = power[0];

        for (int i = 0; i < frequencies.Length; i++)
        {
            periods[i] = frequencies[i] > 0 ? 1.0 / frequencies[i] : double.PositiveInfinity;
            if (power[i] > bestPow)
            {
                bestPow = power[i];
                bestIdx = i;
            }
        }

        Periods = new VectorN(periods);
        BestIndex = bestIdx;
    }
}
