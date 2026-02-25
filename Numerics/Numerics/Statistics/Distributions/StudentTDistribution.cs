using System;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.Distributions
{
    /// <summary>
    /// Student's t-distribution with ν degrees of freedom.
    /// Used for small-sample inference on the mean of a normally distributed population.
    /// PDF = Γ((ν+1)/2) / (√(νπ) Γ(ν/2)) * (1 + x²/ν)^(-(ν+1)/2)
    /// </summary>
    public class StudentTDistribution : IDistribution
    {
        /// <summary>Degrees of freedom (ν).</summary>
        public double DegreesOfFreedom { get; }

        /// <summary>
        /// Creates a Student's t-distribution with the specified degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">Degrees of freedom (ν > 0).</param>
        public StudentTDistribution(double degreesOfFreedom)
        {
            if (degreesOfFreedom <= 0)
                throw new ArgumentException("Degrees of freedom must be positive.");
            DegreesOfFreedom = degreesOfFreedom;
        }

        /// <inheritdoc />
        public double Mean => DegreesOfFreedom > 1 ? 0.0 : double.NaN;

        /// <inheritdoc />
        public double Variance => DegreesOfFreedom > 2
            ? DegreesOfFreedom / (DegreesOfFreedom - 2.0)
            : (DegreesOfFreedom > 1 ? double.PositiveInfinity : double.NaN);

        /// <inheritdoc />
        public double StandardDeviation => Math.Sqrt(Variance);

        /// <inheritdoc />
        public double Pdf(double x)
        {
            double v = DegreesOfFreedom;
            double coeff = Math.Exp(LogGamma((v + 1) / 2.0) - LogGamma(v / 2.0))
                           / Math.Sqrt(v * Math.PI);
            return coeff * Math.Pow(1.0 + x * x / v, -(v + 1) / 2.0);
        }

        /// <summary>
        /// CDF P(X ≤ x) using the regularized incomplete beta function.
        /// </summary>
        public double Cdf(double x)
        {
            double v = DegreesOfFreedom;
            double t2 = x * x;
            double bx = v / (v + t2);

            double ibeta = RegularizedIncompleteBeta(v / 2.0, 0.5, bx);

            if (x >= 0)
                return 1.0 - 0.5 * ibeta;
            else
                return 0.5 * ibeta;
        }

        /// <summary>
        /// Returns the value t such that P(X ≤ t) = p (inverse CDF / quantile function).
        /// Uses bisection on the CDF.
        /// </summary>
        /// <param name="p">Probability in (0, 1).</param>
        public double InverseCdf(double p)
        {
            if (p <= 0 || p >= 1)
                throw new ArgumentOutOfRangeException(nameof(p), "p must be in (0, 1).");

            // Bisection method
            double lo = -1000.0, hi = 1000.0;
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

        /// <inheritdoc />
        public double Sample(RandomGenerator rng)
        {
            // Generate t-distributed sample via ratio of normal / sqrt(chi-squared/v)
            double z = rng.NextGaussian();
            double chi2 = 0;
            int v = (int)Math.Round(DegreesOfFreedom);
            for (int i = 0; i < v; i++)
            {
                double g = rng.NextGaussian();
                chi2 += g * g;
            }
            return z / Math.Sqrt(chi2 / v);
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
        /// Computes the two-tailed p-value for a given t-statistic: P(|T| ≥ |t|).
        /// </summary>
        /// <param name="tStatistic">The observed t value.</param>
        public double TwoTailedPValue(double tStatistic)
        {
            return 2.0 * (1.0 - Cdf(Math.Abs(tStatistic)));
        }

        /// <summary>
        /// Logarithm of the Gamma function using Stirling's approximation (Lanczos).
        /// </summary>
        internal static double LogGamma(double x)
        {
            double[] coef = { 76.18009172947146, -86.50532032941677,
                              24.01409824083091, -1.231739572450155,
                              0.001208650973866179, -0.000005395239384953 };
            double y = x;
            double tmp = x + 5.5;
            tmp -= (x + 0.5) * Math.Log(tmp);
            double ser = 1.000000000190015;
            for (int j = 0; j < 6; j++)
            {
                y += 1.0;
                ser += coef[j] / y;
            }
            return -tmp + Math.Log(2.5066282746310005 * ser / x);
        }

        /// <summary>
        /// Regularized incomplete beta function I_x(a, b) using a continued fraction expansion.
        /// </summary>
        internal static double RegularizedIncompleteBeta(double a, double b, double x)
        {
            if (x < 0 || x > 1) return 0;
            if (x == 0 || x == 1) return x;

            double lnBeta = LogGamma(a) + LogGamma(b) - LogGamma(a + b);
            double front = Math.Exp(Math.Log(x) * a + Math.Log(1.0 - x) * b - lnBeta);

            // Use continued fraction (Lentz's algorithm)
            if (x < (a + 1.0) / (a + b + 2.0))
            {
                return front * BetaContinuedFraction(a, b, x) / a;
            }
            else
            {
                return 1.0 - front * BetaContinuedFraction(b, a, 1.0 - x) / b;
            }
        }

        private static double BetaContinuedFraction(double a, double b, double x)
        {
            const int maxIterations = 200;
            const double epsilon = 1e-14;

            double qab = a + b;
            double qap = a + 1.0;
            double qam = a - 1.0;
            double c = 1.0;
            double d = 1.0 - qab * x / qap;
            if (Math.Abs(d) < 1e-30) d = 1e-30;
            d = 1.0 / d;
            double h = d;

            for (int m = 1; m <= maxIterations; m++)
            {
                int m2 = 2 * m;

                // Even step
                double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
                d = 1.0 + aa * d;
                if (Math.Abs(d) < 1e-30) d = 1e-30;
                c = 1.0 + aa / c;
                if (Math.Abs(c) < 1e-30) c = 1e-30;
                d = 1.0 / d;
                h *= d * c;

                // Odd step
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
                d = 1.0 + aa * d;
                if (Math.Abs(d) < 1e-30) d = 1e-30;
                c = 1.0 + aa / c;
                if (Math.Abs(c) < 1e-30) c = 1e-30;
                d = 1.0 / d;
                double del = d * c;
                h *= del;

                if (Math.Abs(del - 1.0) < epsilon)
                    break;
            }

            return h;
        }
    }
}
