using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The exponential distribution Exp(λ) for modelling time between events.
    /// PDF = λ * exp(-λx) for x ≥ 0.
    /// </summary>
    public class ExponentialDistribution : IDistribution
    {
        /// <summary>Rate parameter (λ).</summary>
        public double Lambda { get; }

        /// <summary>
        /// Creates an exponential distribution with the given rate parameter.
        /// </summary>
        /// <param name="lambda">Rate parameter (λ > 0).</param>
        public ExponentialDistribution(double lambda)
        {
            if (lambda <= 0)
                throw new ArgumentException("lambda must be positive.");

            Lambda = lambda;
        }

        /// <inheritdoc />
        public double Mean => 1.0 / Lambda;

        /// <inheritdoc />
        public double Variance => 1.0 / (Lambda * Lambda);

        /// <inheritdoc />
        public double StandardDeviation => 1.0 / Lambda;

        /// <inheritdoc />
        public double Pdf(double x)
        {
            if (x < 0) return 0.0;
            return Lambda * Math.Exp(-Lambda * x);
        }

        /// <inheritdoc />
        public double Cdf(double x)
        {
            if (x < 0) return 0.0;
            return 1.0 - Math.Exp(-Lambda * x);
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng) => rng.NextExponential(Lambda);

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = rng.NextExponential(Lambda);
            return samples;
        }
    }
}
