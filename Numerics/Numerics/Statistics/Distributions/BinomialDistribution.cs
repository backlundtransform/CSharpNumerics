using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The binomial distribution Bin(n, p): the number of successes
    /// in n independent Bernoulli(p) trials.
    /// PMF = C(n,k) * p^k * (1-p)^(n-k).
    /// </summary>
    public class BinomialDistribution : IDistribution
    {
        /// <summary>Number of trials.</summary>
        public int N { get; }

        /// <summary>Probability of success per trial.</summary>
        public double P { get; }

        /// <summary>
        /// Creates a binomial distribution with n trials and success probability p.
        /// </summary>
        /// <param name="n">Number of trials (n ≥ 0).</param>
        /// <param name="p">Success probability (0 ≤ p ≤ 1).</param>
        public BinomialDistribution(int n, double p)
        {
            if (n < 0)
                throw new ArgumentException("n must be non-negative.");
            if (p < 0 || p > 1)
                throw new ArgumentException("p must be in [0, 1].");

            N = n;
            P = p;
        }

        /// <inheritdoc />
        public double Mean => N * P;

        /// <inheritdoc />
        public double Variance => N * P * (1.0 - P);

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Variance);

        /// <summary>
        /// PMF = C(n,k) * p^k * (1-p)^(n-k).
        /// Returns 0 for non-integer or out-of-range x.
        /// </summary>
        public double Pdf(double x)
        {
            int k = (int)Math.Round(x);
            if (k < 0 || k > N || Math.Abs(x - k) > 1e-9)
                return 0.0;

            return Math.Exp(LogBinomialCoefficient(N, k)
                            + k * Math.Log(P)
                            + (N - k) * Math.Log(1.0 - P));
        }

        /// <summary>
        /// CDF P(X ≤ x) = Σ_{k=0..⌊x⌋} PMF(k).
        /// </summary>
        public double Cdf(double x)
        {
            if (x < 0) return 0.0;
            if (x >= N) return 1.0;

            int upper = (int)Math.Floor(x);
            double sum = 0.0;
            for (int k = 0; k <= upper; k++)
                sum += Pdf(k);
            return Math.Min(sum, 1.0);
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng) => rng.NextBinomial(N, P);

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = rng.NextBinomial(N, P);
            return samples;
        }

        private static double LogBinomialCoefficient(int n, int k)
        {
            return LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k);
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
