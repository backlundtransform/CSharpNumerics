using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// Common interface for probability distributions.
    /// Every distribution can report its parameters, compute PDF/PMF and CDF,
    /// and generate random samples.
    /// </summary>
    public interface IDistribution
    {
        /// <summary>
        /// The theoretical mean (expected value) of the distribution.
        /// </summary>
        double Mean { get; }

        /// <summary>
        /// The theoretical variance of the distribution.
        /// </summary>
        double Variance { get; }

        /// <summary>
        /// The theoretical standard deviation (√Variance).
        /// </summary>
        double StandardDeviation { get; }

        /// <summary>
        /// Evaluates the probability density function (continuous)
        /// or probability mass function (discrete) at x.
        /// </summary>
        /// <param name="x">The point at which to evaluate.</param>
        double Pdf(double x);

        /// <summary>
        /// Evaluates the cumulative distribution function P(X ≤ x).
        /// </summary>
        /// <param name="x">The point at which to evaluate.</param>
        double Cdf(double x);

        /// <summary>
        /// Generates a single random sample from this distribution.
        /// </summary>
        /// <param name="rng">The random generator to use.</param>
        double Sample(RandomGenerator rng);

        /// <summary>
        /// Generates n random samples from this distribution.
        /// </summary>
        /// <param name="rng">The random generator to use.</param>
        /// <param name="count">Number of samples.</param>
        double[] Samples(RandomGenerator rng, int count);
    }
}
