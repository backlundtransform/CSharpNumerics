using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The Bernoulli distribution Bernoulli(p).
    /// P(X=1) = p, P(X=0) = 1-p.
    /// </summary>
    public class BernoulliDistribution : IDistribution
    {
        /// <summary>Probability of success.</summary>
        public double P { get; }

        /// <summary>
        /// Creates a Bernoulli distribution with the given success probability.
        /// </summary>
        /// <param name="p">Probability of success (0 ≤ p ≤ 1).</param>
        public BernoulliDistribution(double p)
        {
            if (p < 0 || p > 1)
                throw new ArgumentException("p must be in [0, 1].");

            P = p;
        }

        /// <inheritdoc />
        public double Mean => P;

        /// <inheritdoc />
        public double Variance => P * (1.0 - P);

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Variance);

        /// <summary>
        /// PMF: returns P for x=1, (1-P) for x=0, 0 otherwise.
        /// </summary>
        public double Pdf(double x)
        {
            if (Math.Abs(x - 1.0) < 1e-9) return P;
            if (Math.Abs(x) < 1e-9) return 1.0 - P;
            return 0.0;
        }

        /// <summary>
        /// CDF: 0 for x &lt; 0, (1-P) for 0 ≤ x &lt; 1, 1 for x ≥ 1.
        /// </summary>
        public double Cdf(double x)
        {
            if (x < 0) return 0.0;
            if (x < 1) return 1.0 - P;
            return 1.0;
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng) => rng.NextBernoulli(P);

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = rng.NextBernoulli(P);
            return samples;
        }
    }
}
