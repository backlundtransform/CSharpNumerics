using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The F-distribution F(d1, d2) with d1 and d2 degrees of freedom.
    /// Arises as the ratio of two scaled chi-squared distributions.
    /// Used in ANOVA and variance-ratio tests.
    /// </summary>
    public class FDistribution : IDistribution
    {
        /// <summary>Numerator degrees of freedom (d1).</summary>
        public double D1 { get; }

        /// <summary>Denominator degrees of freedom (d2).</summary>
        public double D2 { get; }

        /// <summary>
        /// Creates an F-distribution with the specified degrees of freedom.
        /// </summary>
        /// <param name="d1">Numerator degrees of freedom (d1 > 0).</param>
        /// <param name="d2">Denominator degrees of freedom (d2 > 0).</param>
        public FDistribution(double d1, double d2)
        {
            if (d1 <= 0)
                throw new ArgumentException("d1 must be positive.");
            if (d2 <= 0)
                throw new ArgumentException("d2 must be positive.");
            D1 = d1;
            D2 = d2;
        }

        /// <inheritdoc />
        public double Mean => D2 > 2 ? D2 / (D2 - 2.0) : double.NaN;

        /// <inheritdoc />
        public double Variance => D2 > 4
            ? (2.0 * D2 * D2 * (D1 + D2 - 2.0)) / (D1 * (D2 - 2.0) * (D2 - 2.0) * (D2 - 4.0))
            : double.NaN;

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Variance);

        /// <inheritdoc />
        public double Pdf(double x)
        {
            if (x <= 0) return 0;

            double d1 = D1, d2 = D2;
            double lnNum = (d1 / 2.0) * Math.Log(d1) + (d2 / 2.0) * Math.Log(d2)
                         + ((d1 / 2.0) - 1.0) * Math.Log(x);
            double lnDen = ((d1 + d2) / 2.0) * Math.Log(d1 * x + d2)
                         + LogBeta(d1 / 2.0, d2 / 2.0);
            return Math.Exp(lnNum - lnDen);
        }

        /// <summary>
        /// CDF P(X ≤ x) using the regularized incomplete beta function.
        /// F CDF = I_{d1*x/(d1*x+d2)}(d1/2, d2/2)
        /// </summary>
        public double Cdf(double x)
        {
            if (x <= 0) return 0;
            double bx = D1 * x / (D1 * x + D2);
            return StudentTDistribution.RegularizedIncompleteBeta(D1 / 2.0, D2 / 2.0, bx);
        }

        /// <summary>
        /// Returns the value x such that P(X ≤ x) = p (inverse CDF).
        /// Uses bisection.
        /// </summary>
        /// <param name="p">Probability in (0, 1).</param>
        public double InverseCdf(double p)
        {
            if (p <= 0 || p >= 1)
                throw new ArgumentOutOfRangeException(nameof(p), "p must be in (0, 1).");

            double lo = 0.0;
            double hi = 100.0;
            while (Cdf(hi) < p) hi *= 2;

            for (int i = 0; i < 100; i++)
            {
                double mid = (lo + hi) / 2.0;
                if (Cdf(mid) < p)
                    lo = mid;
                else
                    hi = mid;
            }
            return (lo + hi) / 2.0;
        }

        /// <summary>
        /// Returns the upper-tail p-value: P(X ≥ f) = 1 - CDF(f).
        /// </summary>
        public double UpperTailPValue(double fStatistic)
        {
            return 1.0 - Cdf(fStatistic);
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng)
        {
            // F = (χ²₁/d1) / (χ²₂/d2)
            var chi1 = new ChiSquaredDistribution(D1);
            var chi2 = new ChiSquaredDistribution(D2);
            return (chi1.Sample(rng) / D1) / (chi2.Sample(rng) / D2);
        }

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = Sample(rng);
            return samples;
        }

        private static double LogBeta(double a, double b)
        {
            return StudentTDistribution.LogGamma(a) + StudentTDistribution.LogGamma(b)
                   - StudentTDistribution.LogGamma(a + b);
        }
    }
}
