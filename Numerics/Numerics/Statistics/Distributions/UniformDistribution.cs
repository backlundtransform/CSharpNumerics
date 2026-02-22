using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The continuous uniform distribution over [a, b].
    /// PDF = 1/(b-a) on [a,b], 0 elsewhere.
    /// </summary>
    public class UniformDistribution : IDistribution
    {
        /// <summary>Lower bound (inclusive).</summary>
        public double A { get; }

        /// <summary>Upper bound (exclusive for sampling, inclusive for PDF).</summary>
        public double B { get; }

        /// <summary>
        /// Creates a uniform distribution over [a, b].
        /// </summary>
        /// <param name="a">Lower bound.</param>
        /// <param name="b">Upper bound.</param>
        public UniformDistribution(double a = 0.0, double b = 1.0)
        {
            if (a >= b)
                throw new ArgumentException("a must be less than b.");

            A = a;
            B = b;
        }

        /// <inheritdoc />
        public double Mean => (A + B) / 2.0;

        /// <inheritdoc />
        public double Variance => Math.Pow(B - A, 2) / 12.0;

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Variance);

        /// <inheritdoc />
        public double Pdf(double x)
        {
            if (x < A || x > B) return 0.0;
            return 1.0 / (B - A);
        }

        /// <inheritdoc />
        public double Cdf(double x)
        {
            if (x < A) return 0.0;
            if (x >= B) return 1.0;
            return (x - A) / (B - A);
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng) => rng.NextUniform(A, B);

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count) => rng.UniformSamples(count, A, B);
    }
}
