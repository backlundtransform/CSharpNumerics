using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Statistics.MonteCarlo
{
    /// <summary>
    /// Stores the aggregated results of a Monte Carlo simulation.
    /// Provides descriptive statistics, percentiles, histograms and
    /// access to the raw samples for further analysis.
    /// </summary>
    public class MonteCarloResult
    {
        private readonly double[] _sortedSamples;

        /// <summary>
        /// All raw sample values in the order they were generated.
        /// </summary>
        public double[] Samples { get; }

        /// <summary>Number of iterations (samples) in the simulation.</summary>
        public int Count => Samples.Length;

        /// <summary>Sample mean.</summary>
        public double Mean { get; }

        /// <summary>Sample variance (using n-1 divisor).</summary>
        public double Variance { get; }

        /// <summary>Sample standard deviation.</summary>
        public double StandardDeviation { get; }

        /// <summary>Minimum observed value.</summary>
        public double Min { get; }

        /// <summary>Maximum observed value.</summary>
        public double Max { get; }

        /// <summary>Median (50th percentile).</summary>
        public double Median { get; }

        /// <summary>
        /// Creates a result from the provided sample array.
        /// All statistics are computed eagerly.
        /// </summary>
        /// <param name="samples">The raw sample data.</param>
        public MonteCarloResult(double[] samples)
        {
            if (samples == null || samples.Length == 0)
                throw new ArgumentException("samples must not be null or empty.");

            Samples = samples;
            _sortedSamples = (double[])samples.Clone();
            Array.Sort(_sortedSamples);

            Mean = samples.Average();
            Min = _sortedSamples[0];
            Max = _sortedSamples[_sortedSamples.Length - 1];
            Median = Percentile(50);

            double sumSq = 0;
            for (int i = 0; i < samples.Length; i++)
            {
                double diff = samples[i] - Mean;
                sumSq += diff * diff;
            }
            Variance = samples.Length > 1 ? sumSq / (samples.Length - 1) : 0.0;
            StandardDeviation = Math.Sqrt(Variance);
        }

        /// <summary>
        /// Returns the value at the given percentile (0–100) using linear interpolation.
        /// </summary>
        /// <param name="p">Percentile (0–100).</param>
        public double Percentile(double p)
        {
            if (p < 0 || p > 100)
                throw new ArgumentOutOfRangeException(nameof(p), "Percentile must be in [0, 100].");

            if (_sortedSamples.Length == 1)
                return _sortedSamples[0];

            double rank = (p / 100.0) * (_sortedSamples.Length - 1);
            int lower = (int)Math.Floor(rank);
            int upper = (int)Math.Ceiling(rank);

            if (lower == upper)
                return _sortedSamples[lower];

            double fraction = rank - lower;
            return _sortedSamples[lower] * (1.0 - fraction) + _sortedSamples[upper] * fraction;
        }

        /// <summary>
        /// Returns a symmetric confidence interval [lower, upper] around the mean
        /// based on the empirical percentiles.
        /// </summary>
        /// <param name="confidenceLevel">Confidence level, e.g. 0.95 for 95%.</param>
        public (double lower, double upper) ConfidenceInterval(double confidenceLevel)
        {
            if (confidenceLevel <= 0 || confidenceLevel >= 1)
                throw new ArgumentOutOfRangeException(nameof(confidenceLevel), "Must be in (0, 1).");

            double alpha = (1.0 - confidenceLevel) / 2.0;
            double lower = Percentile(alpha * 100.0);
            double upper = Percentile((1.0 - alpha) * 100.0);
            return (lower, upper);
        }

        /// <summary>
        /// Returns the fraction of samples that satisfy the predicate.
        /// Useful for computing probabilities, e.g. P(X > threshold).
        /// </summary>
        /// <param name="predicate">A function that returns true for "events of interest".</param>
        public double Probability(Func<double, bool> predicate)
        {
            int hits = 0;
            for (int i = 0; i < Samples.Length; i++)
            {
                if (predicate(Samples[i]))
                    hits++;
            }
            return (double)hits / Samples.Length;
        }

        /// <summary>
        /// Generates a histogram of the samples with the specified number of bins.
        /// Returns an array of (binCenter, count) tuples.
        /// </summary>
        /// <param name="bins">Number of bins (default 20).</param>
        public (double binCenter, int count)[] Histogram(int bins = 20)
        {
            if (bins < 1)
                throw new ArgumentException("bins must be at least 1.");

            double range = Max - Min;
            if (range == 0)
                return new[] { (Min, Count) };

            double binWidth = range / bins;
            var histogram = new (double binCenter, int count)[bins];

            for (int i = 0; i < bins; i++)
                histogram[i] = (Min + binWidth * (i + 0.5), 0);

            for (int i = 0; i < Samples.Length; i++)
            {
                int bin = (int)((Samples[i] - Min) / binWidth);
                if (bin >= bins) bin = bins - 1;
                histogram[bin] = (histogram[bin].binCenter, histogram[bin].count + 1);
            }

            return histogram;
        }

        /// <summary>
        /// Standard error of the mean (SEM = σ / √n).
        /// Indicates the precision of the Monte Carlo estimate.
        /// </summary>
        public double StandardError => StandardDeviation / Math.Sqrt(Count);
    }
}
