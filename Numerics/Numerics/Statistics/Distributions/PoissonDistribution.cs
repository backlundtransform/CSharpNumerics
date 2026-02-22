using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The Poisson distribution Poisson(λ) for modelling counts of events.
    /// PMF = (λ^k * e^(-λ)) / k!.
    /// </summary>
    public class PoissonDistribution : IDistribution
    {
        /// <summary>Expected number of events (λ).</summary>
        public double Lambda { get; }

        /// <summary>
        /// Creates a Poisson distribution with the given mean λ.
        /// </summary>
        /// <param name="lambda">Expected number of events (λ > 0).</param>
        public PoissonDistribution(double lambda)
        {
            if (lambda <= 0)
                throw new ArgumentException("lambda must be positive.");

            Lambda = lambda;
        }

        /// <inheritdoc />
        public double Mean => Lambda;

        /// <inheritdoc />
        public double Variance => Lambda;

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Lambda);

        /// <summary>
        /// Evaluates the probability mass function at integer k.
        /// For non-integer x the value is 0.
        /// </summary>
        public double Pdf(double x)
        {
            int k = (int)Math.Round(x);
            if (k < 0 || Math.Abs(x - k) > 1e-9) return 0.0;

            return Math.Exp(k * Math.Log(Lambda) - Lambda - LogFactorial(k));
        }

        /// <summary>
        /// Evaluates the CDF P(X ≤ x) = Σ_{k=0..⌊x⌋} PMF(k).
        /// </summary>
        public double Cdf(double x)
        {
            if (x < 0) return 0.0;
            int n = (int)Math.Floor(x);
            double sum = 0.0;

            for (int k = 0; k <= n; k++)
                sum += Pdf(k);

            return Math.Min(sum, 1.0);
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng) => rng.NextPoisson(Lambda);

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = rng.NextPoisson(Lambda);
            return samples;
        }

        private static double LogFactorial(int n)
        {
            double result = 0;
            for (int i = 2; i <= n; i++)
                result += Math.Log(i);
            return result;
        }
    }
}
