using System;

namespace CSharpNumerics.Statistics.Random
{
    /// <summary>
    /// Provides advanced random sampling methods built on top of <see cref="System.Random"/>.
    /// Supports uniform, Gaussian (Box-Muller), exponential, and other sampling strategies.
    /// Thread-safe when each thread uses its own instance.
    /// </summary>
    public class RandomGenerator
    {
        private readonly System.Random _rng;
        private double? _spareGaussian;

        /// <summary>
        /// Initializes a new <see cref="RandomGenerator"/> with a random seed.
        /// </summary>
        public RandomGenerator()
        {
            _rng = new System.Random();
        }

        /// <summary>
        /// Initializes a new <see cref="RandomGenerator"/> with the specified seed for reproducible results.
        /// </summary>
        /// <param name="seed">The seed value.</param>
        public RandomGenerator(int seed)
        {
            _rng = new System.Random(seed);
        }

        /// <summary>
        /// Returns a random double in [0, 1).
        /// </summary>
        public double NextDouble() => _rng.NextDouble();

        /// <summary>
        /// Returns a random integer in [0, maxExclusive).
        /// </summary>
        /// <param name="maxExclusive">The exclusive upper bound.</param>
        public int NextInt(int maxExclusive) => _rng.Next(maxExclusive);

        /// <summary>
        /// Returns a random integer in [minInclusive, maxExclusive).
        /// </summary>
        /// <param name="minInclusive">The inclusive lower bound.</param>
        /// <param name="maxExclusive">The exclusive upper bound.</param>
        public int NextInt(int minInclusive, int maxExclusive) => _rng.Next(minInclusive, maxExclusive);

        /// <summary>
        /// Returns a uniformly distributed random double in [min, max).
        /// </summary>
        /// <param name="min">Lower bound (inclusive).</param>
        /// <param name="max">Upper bound (exclusive).</param>
        public double NextUniform(double min, double max)
        {
            if (min >= max)
                throw new ArgumentException("min must be less than max.");

            return min + _rng.NextDouble() * (max - min);
        }

        /// <summary>
        /// Returns a normally distributed random value using the Box-Muller transform.
        /// Uses the spare value from the previous call when available.
        /// </summary>
        /// <param name="mean">Mean of the distribution (default 0).</param>
        /// <param name="standardDeviation">Standard deviation (default 1).</param>
        public double NextGaussian(double mean = 0.0, double standardDeviation = 1.0)
        {
            if (standardDeviation < 0)
                throw new ArgumentException("standardDeviation must be non-negative.");

            if (_spareGaussian.HasValue)
            {
                var spare = _spareGaussian.Value;
                _spareGaussian = null;
                return mean + standardDeviation * spare;
            }

            double u, v, s;
            do
            {
                u = 2.0 * _rng.NextDouble() - 1.0;
                v = 2.0 * _rng.NextDouble() - 1.0;
                s = u * u + v * v;
            } while (s >= 1.0 || s == 0.0);

            var factor = Math.Sqrt(-2.0 * Math.Log(s) / s);
            _spareGaussian = v * factor;
            return mean + standardDeviation * (u * factor);
        }

        /// <summary>
        /// Returns an exponentially distributed random value with the given rate parameter λ.
        /// Uses inverse transform sampling: X = -ln(U) / λ.
        /// </summary>
        /// <param name="lambda">Rate parameter (λ > 0).</param>
        public double NextExponential(double lambda)
        {
            if (lambda <= 0)
                throw new ArgumentException("lambda must be positive.");

            return -Math.Log(1.0 - _rng.NextDouble()) / lambda;
        }

        /// <summary>
        /// Returns a Poisson-distributed random integer with the given mean λ.
        /// Uses Knuth's algorithm for small λ and the rejection method for large λ.
        /// </summary>
        /// <param name="lambda">Expected number of events (λ > 0).</param>
        public int NextPoisson(double lambda)
        {
            if (lambda <= 0)
                throw new ArgumentException("lambda must be positive.");

            if (lambda < 30)
            {
                // Knuth's algorithm
                double l = Math.Exp(-lambda);
                int k = 0;
                double p = 1.0;

                do
                {
                    k++;
                    p *= _rng.NextDouble();
                } while (p > l);

                return k - 1;
            }
            else
            {
                // Rejection method (normal approximation with correction)
                double c = 0.767 - 3.36 / lambda;
                double beta = Math.PI / Math.Sqrt(3.0 * lambda);
                double alpha = beta * lambda;
                double k = Math.Log(c) - lambda - Math.Log(beta);

                while (true)
                {
                    double u = _rng.NextDouble();
                    double x = (alpha - Math.Log((1 - u) / u)) / beta;
                    int n = (int)Math.Floor(x + 0.5);
                    if (n < 0) continue;

                    double v = _rng.NextDouble();
                    double y = alpha - beta * n;
                    double lhs = y + Math.Log(v / Math.Pow(1.0 + Math.Exp(y), 2));
                    double rhs = k + n * Math.Log(lambda) - LogFactorial(n);

                    if (lhs <= rhs)
                        return n;
                }
            }
        }

        /// <summary>
        /// Returns a Bernoulli trial result: 1 with probability p, 0 otherwise.
        /// </summary>
        /// <param name="p">Probability of success (0 ≤ p ≤ 1).</param>
        public int NextBernoulli(double p)
        {
            if (p < 0 || p > 1)
                throw new ArgumentException("p must be in [0, 1].");

            return _rng.NextDouble() < p ? 1 : 0;
        }

        /// <summary>
        /// Returns a binomially distributed random integer: the number of successes
        /// in n independent Bernoulli trials each with success probability p.
        /// </summary>
        /// <param name="n">Number of trials.</param>
        /// <param name="p">Probability of success per trial.</param>
        public int NextBinomial(int n, double p)
        {
            if (n < 0)
                throw new ArgumentException("n must be non-negative.");
            if (p < 0 || p > 1)
                throw new ArgumentException("p must be in [0, 1].");

            int successes = 0;
            for (int i = 0; i < n; i++)
            {
                if (_rng.NextDouble() < p)
                    successes++;
            }
            return successes;
        }

        /// <summary>
        /// Generates an array of n samples from a standard normal distribution N(0,1).
        /// </summary>
        /// <param name="count">Number of samples.</param>
        public double[] GaussianSamples(int count)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = NextGaussian();
            return samples;
        }

        /// <summary>
        /// Generates an array of n uniform samples in [min, max).
        /// </summary>
        /// <param name="count">Number of samples.</param>
        /// <param name="min">Lower bound.</param>
        /// <param name="max">Upper bound.</param>
        public double[] UniformSamples(int count, double min = 0.0, double max = 1.0)
        {
            var samples = new double[count];
            for (int i = 0; i < count; i++)
                samples[i] = NextUniform(min, max);
            return samples;
        }

        /// <summary>
        /// Shuffles the input array in-place using the Fisher-Yates algorithm.
        /// </summary>
        /// <typeparam name="T">Element type.</typeparam>
        /// <param name="array">The array to shuffle.</param>
        public void Shuffle<T>(T[] array)
        {
            for (int i = array.Length - 1; i > 0; i--)
            {
                int j = _rng.Next(i + 1);
                var tmp = array[i];
                array[i] = array[j];
                array[j] = tmp;
            }
        }

        /// <summary>
        /// Selects k elements from the array uniformly at random without replacement.
        /// </summary>
        /// <typeparam name="T">Element type.</typeparam>
        /// <param name="source">The source array.</param>
        /// <param name="k">Number of elements to select.</param>
        public T[] Sample<T>(T[] source, int k)
        {
            if (k > source.Length)
                throw new ArgumentException("k cannot exceed the source length.");

            var copy = (T[])source.Clone();
            Shuffle(copy);

            var result = new T[k];
            Array.Copy(copy, result, k);
            return result;
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
