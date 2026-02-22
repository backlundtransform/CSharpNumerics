using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// The normal (Gaussian) distribution N(μ, σ²).
    /// PDF = (1 / (σ√(2π))) * exp(-(x-μ)² / (2σ²)).
    /// CDF is approximated using the Abramowitz–Stegun rational approximation for the error function.
    /// </summary>
    public class NormalDistribution : IDistribution
    {
        /// <summary>Mean (μ) of the distribution.</summary>
        public double Mu { get; }

        /// <summary>Standard deviation (σ) of the distribution.</summary>
        public double Sigma { get; }

        /// <summary>
        /// Creates a normal distribution with the specified mean and standard deviation.
        /// </summary>
        /// <param name="mu">Mean (default 0).</param>
        /// <param name="sigma">Standard deviation (default 1).</param>
        public NormalDistribution(double mu = 0.0, double sigma = 1.0)
        {
            if (sigma <= 0)
                throw new ArgumentException("sigma must be positive.");

            Mu = mu;
            Sigma = sigma;
        }

        /// <inheritdoc />
        public double Mean => Mu;

        /// <inheritdoc />
        public double Variance => Sigma * Sigma;

        /// <inheritdoc />
        public double StandardDeviation => Sigma;

        /// <inheritdoc />
        public double Pdf(double x)
        {
            double z = (x - Mu) / Sigma;
            return Math.Exp(-0.5 * z * z) / (Sigma * Math.Sqrt(2.0 * Math.PI));
        }

        /// <inheritdoc />
        public double Cdf(double x)
        {
            double z = (x - Mu) / (Sigma * Math.Sqrt(2.0));
            return 0.5 * (1.0 + Erf(z));
        }

        /// <inheritdoc />
        public double Sample(RandomGenerator rng) => rng.NextGaussian(Mu, Sigma);

        /// <inheritdoc />
        public double[] Samples(RandomGenerator rng, int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = rng.NextGaussian(Mu, Sigma);
            return samples;
        }

        /// <summary>
        /// Returns the value z such that P(X ≤ z) = p (inverse CDF / quantile function).
        /// Uses the rational approximation by Peter Acklam.
        /// </summary>
        /// <param name="p">Probability in (0, 1).</param>
        public double InverseCdf(double p)
        {
            if (p <= 0 || p >= 1)
                throw new ArgumentOutOfRangeException(nameof(p), "p must be in (0, 1).");

            // Standard normal quantile
            double z = RationalApproxInvCdf(p);
            return Mu + Sigma * z;
        }

        /// <summary>
        /// Abramowitz–Stegun approximation for the error function erf(x).
        /// Maximum error ≈ 1.5 × 10⁻⁷.
        /// </summary>
        private static double Erf(double x)
        {
            double sign = x < 0 ? -1.0 : 1.0;
            x = Math.Abs(x);

            const double a1 = 0.254829592;
            const double a2 = -0.284496736;
            const double a3 = 1.421413741;
            const double a4 = -1.453152027;
            const double a5 = 1.061405429;
            const double p = 0.3275911;

            double t = 1.0 / (1.0 + p * x);
            double y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

            return sign * y;
        }

        /// <summary>
        /// Rational approximation for the inverse CDF of the standard normal distribution.
        /// Peter Acklam's algorithm — relative error < 1.15 × 10⁻⁹.
        /// </summary>
        private static double RationalApproxInvCdf(double p)
        {
            const double pLow = 0.02425;
            const double pHigh = 1.0 - pLow;

            double[] a = { -3.969683028665376e+01, 2.209460984245205e+02,
                           -2.759285104469687e+02, 1.383577518672690e+02,
                           -3.066479806614716e+01, 2.506628277459239e+00 };

            double[] b = { -5.447609879822406e+01, 1.615858368580409e+02,
                           -1.556989798598866e+02, 6.680131188771972e+01,
                           -1.328068155288572e+01 };

            double[] c = { -7.784894002430293e-03, -3.223964580411365e-01,
                           -2.400758277161838e+00, -2.549732539343734e+00,
                            4.374664141464968e+00, 2.938163982698783e+00 };

            double[] d = {  7.784695709041462e-03, 3.224671290700398e-01,
                            2.445134137142996e+00, 3.754408661907416e+00 };

            double q, r;

            if (p < pLow)
            {
                q = Math.Sqrt(-2.0 * Math.Log(p));
                return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                       ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
            }
            else if (p <= pHigh)
            {
                q = p - 0.5;
                r = q * q;
                return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
                       (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
            }
            else
            {
                q = Math.Sqrt(-2.0 * Math.Log(1.0 - p));
                return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
            }
        }
    }
}
