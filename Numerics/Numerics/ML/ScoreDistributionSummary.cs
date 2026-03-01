using CSharpNumerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML;

/// <summary>
/// Descriptive statistics summary for a collection of scores.
/// Used by experiment results to characterize how scores are distributed
/// across all tested pipeline configurations.
///
/// All properties are computed from the library's
/// <see cref="DescriptiveStatisticsExtensions"/>.
/// </summary>
public class ScoreDistributionSummary
{
    /// <summary>Number of scores (pipeline configurations).</summary>
    public int Count { get; }

    /// <summary>Arithmetic mean score.</summary>
    public double Mean { get; }

    /// <summary>Median score (50th percentile).</summary>
    public double Median { get; }

    /// <summary>Population standard deviation of scores.</summary>
    public double StandardDeviation { get; }

    /// <summary>Minimum score.</summary>
    public double Min { get; }

    /// <summary>Maximum score.</summary>
    public double Max { get; }

    /// <summary>Score range (Max - Min).</summary>
    public double Range { get; }

    /// <summary>25th percentile (Q1).</summary>
    public double P25 { get; }

    /// <summary>75th percentile (Q3).</summary>
    public double P75 { get; }

    /// <summary>Interquartile range (Q3 - Q1).</summary>
    public double InterquartileRange { get; }

    /// <summary>
    /// Fisher's sample skewness.
    /// Positive = right-tailed (many low scores, few high).
    /// Negative = left-tailed (many high scores, few low).
    /// Requires at least 3 scores.
    /// </summary>
    public double Skewness { get; }

    /// <summary>
    /// Sample excess kurtosis.
    /// High = heavy tails (score outliers).
    /// 0 = normal-like distribution.
    /// Requires at least 4 scores.
    /// </summary>
    public double Kurtosis { get; }

    /// <summary>
    /// Confidence interval (lower, upper) for the mean score
    /// at the specified confidence level.
    /// </summary>
    public (double Lower, double Upper) ConfidenceInterval { get; }

    /// <summary>Confidence level used for the interval (default 0.95).</summary>
    public double ConfidenceLevel { get; }

    /// <summary>
    /// Build a summary from a collection of scores.
    /// </summary>
    /// <param name="scores">The score values to summarize.</param>
    /// <param name="confidenceLevel">Confidence level for the interval (default 0.95).</param>
    public ScoreDistributionSummary(IEnumerable<double> scores, double confidenceLevel = 0.95)
    {
        var list = scores.ToList();
        Count = list.Count;
        ConfidenceLevel = confidenceLevel;

        if (Count == 0)
            return;

        Func<double, double> identity = x => x;

        Mean = list.Average();
        Median = list.Median(identity);
        StandardDeviation = list.StandardDeviation(identity);
        Min = list.Min();
        Max = list.Max();
        Range = list.Range(identity);

        if (Count >= 2)
        {
            P25 = list.Percentile(identity, 25);
            P75 = list.Percentile(identity, 75);
            InterquartileRange = P75 - P25;
            ConfidenceInterval = list.ConfidenceIntervals(identity, confidenceLevel);
        }

        if (Count >= 3)
            Skewness = list.Skewness(identity);

        if (Count >= 4)
            Kurtosis = list.Kurtosis(identity);
    }

    public override string ToString()
    {
        return $"N={Count}, Mean={Mean:F4}, Median={Median:F4}, " +
               $"StdDev={StandardDeviation:F4}, IQR={InterquartileRange:F4}, " +
               $"Range=[{Min:F4}, {Max:F4}], " +
               $"Skew={Skewness:F3}, Kurt={Kurtosis:F3}";
    }
}
