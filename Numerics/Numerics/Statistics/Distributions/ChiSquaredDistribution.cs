using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The chi-squared distribution χ²(k) with k degrees of freedom.
    /// This is the distribution of the sum of k independent standard normal random variables squared.
    /// PDF = (x^(k/2-1) * e^(-x/2)) / (2^(k/2) * Γ(k/2)), for x > 0.
    /// </summary>
    public class ChiSquaredDistribution : IDistribution
    {
        /// <summary>Degrees of freedom (k).</summary>
        public double DegreesOfFreedom { get; }

        /// <summary>
        /// Creates a chi-squared distribution with the specified degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom (k > 0).</param>
        public ChiSquaredDistribution(double degreesOfFreedom)
        {
            if (degreesOfFreedom <= 0)
                throw new ArgumentException("Degrees of freedom must be positive.");
            DegreesOfFreedom = degreesOfFreedom;
        }

        /// <inheritdoc />
        public double Mean => DegreesOfFreedom;

        /// <inheritdoc />
        public double Variance => 2.0 * DegreesOfFreedom;

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Variance);

        /// <inheritdoc />
        public double Pdf(double x)
        {
            if (x <= 0) return 0;
            double k2 = DegreesOfFreedom / 2.0;
            return Math.Exp((k2 - 1.0) * Math.Log(x) - x / 2.0 - k2 * Math.Log(2.0) - StudentTDistribution.LogGamma(k2));
        }

        /// <summary>
        /// CDF P(X ≤ x) using the regularized lower incomplete gamma function.
        /// P(X ≤ x) = γ(k/2, x/2) / Γ(k/2) = P(k/2, x/2)
        /// </summary>
        public double Cdf(double x)
        {
            if (x <= 0) return 0;
            return LowerRegularizedGamma(DegreesOfFreedom / 2.0, x / 2.0);
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

            double lo = 0.0, hi = DegreesOfFreedom + 10.0 * Math.Sqrt(2.0 * DegreesOfFreedom);
            // Expand upper bound if needed
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
        /// Returns the upper-tail p-value: P(X ≥ x) = 1 - CDF(x).
        /// </summary>
        public double UpperTailPValue(double chiSquaredStatistic)
        {
            return 1.0 - Cdf(chiSquaredStatistic);
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng)
        {
            double sum = 0;
            int k = (int)Math.Round(DegreesOfFreedom);
            for (int i = 0; i < k; i++)
            {
                double z = rng.NextGaussian();
                sum += z * z;
            }
            return sum;
        }

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = Sample(rng);
            return samples;
        }

        /// <summary>
        /// Regularized lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a).
        /// Uses series expansion for x &lt; a+1, continued fraction otherwise.
        /// </summary>
        internal static double LowerRegularizedGamma(double a, double x)
        {
            if (x < 0) return 0;
            if (x == 0) return 0;

            if (x < a + 1.0)
            {
                // Series expansion
                return GammaSeries(a, x);
            }
            else
            {
                // Continued fraction
                return 1.0 - GammaContinuedFraction(a, x);
            }
        }

        private static double GammaSeries(double a, double x)
        {
            const int maxIterations = 200;
            const double epsilon = 1e-14;

            double lnGamma = StudentTDistribution.LogGamma(a);
            double ap = a;
            double sum = 1.0 / a;
            double del = sum;

            for (int n = 1; n <= maxIterations; n++)
            {
                ap += 1.0;
                del *= x / ap;
                sum += del;
                if (Math.Abs(del) < Math.Abs(sum) * epsilon)
                    break;
            }

            return sum * Math.Exp(-x + a * Math.Log(x) - lnGamma);
        }

        private static double GammaContinuedFraction(double a, double x)
        {
            const int maxIterations = 200;
            const double epsilon = 1e-14;

            double lnGamma = StudentTDistribution.LogGamma(a);
            double b = x + 1.0 - a;
            double c = 1.0 / 1e-30;
            double d = 1.0 / b;
            double h = d;

            for (int i = 1; i <= maxIterations; i++)
            {
                double an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (Math.Abs(d) < 1e-30) d = 1e-30;
                c = b + an / c;
                if (Math.Abs(c) < 1e-30) c = 1e-30;
                d = 1.0 / d;
                double del = d * c;
                h *= del;
                if (Math.Abs(del - 1.0) < epsilon)
                    break;
            }

            return Math.Exp(-x + a * Math.Log(x) - lnGamma) * h;
        }
    }
}
