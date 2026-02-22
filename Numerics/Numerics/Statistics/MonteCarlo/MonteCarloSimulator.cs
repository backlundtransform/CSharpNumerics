using System;
using System.Collections.Generic;
using CSharpNumerics.Statistics.Distributions;
using CSharpNumerics.Statistics.Random;

namespace CSharpNumerics.Statistics.MonteCarlo
{
    /// <summary>
    /// A general-purpose Monte Carlo simulation engine.
    /// Runs a user-defined model (or a plain function) many times with
    /// stochastic inputs and collects the results into a <see cref="MonteCarloResult"/>.
    /// 
    /// Supports three usage patterns:
    /// <list type="number">
    ///   <item>Custom model via <see cref="IMonteCarloModel"/>.</item>
    ///   <item>Inline function that receives a <see cref="RandomGenerator"/>.</item>
    ///   <item>Integration estimation — computes E[f(X)] where X ~ distribution.</item>
    /// </list>
    /// 
    /// <example>
    /// <code>
    /// // Estimate π by sampling random points in the unit square
    /// var sim = new MonteCarloSimulator(seed: 42);
    /// var result = sim.Run(rng =>
    /// {
    ///     double x = rng.NextUniform(0, 1);
    ///     double y = rng.NextUniform(0, 1);
    ///     return x * x + y * y &lt;= 1.0 ? 1.0 : 0.0;
    /// }, iterations: 1_000_000);
    /// double piEstimate = result.Mean * 4.0;
    /// </code>
    /// </example>
    /// </summary>
    public class MonteCarloSimulator
    {
        private readonly RandomGenerator _rng;

        /// <summary>
        /// Creates a simulator with a random seed.
        /// </summary>
        public MonteCarloSimulator()
        {
            _rng = new RandomGenerator();
        }

        /// <summary>
        /// Creates a simulator with the specified seed for reproducible results.
        /// </summary>
        /// <param name="seed">Random seed.</param>
        public MonteCarloSimulator(int seed)
        {
            _rng = new RandomGenerator(seed);
        }

        /// <summary>
        /// Creates a simulator using an existing <see cref="RandomGenerator"/> instance.
        /// </summary>
        /// <param name="rng">The random generator.</param>
        public MonteCarloSimulator(RandomGenerator rng)
        {
            _rng = rng ?? throw new ArgumentNullException(nameof(rng));
        }

        /// <summary>
        /// Runs a simulation using an <see cref="IMonteCarloModel"/>.
        /// Each iteration calls <see cref="IMonteCarloModel.Evaluate"/>.
        /// </summary>
        /// <param name="model">The model to evaluate.</param>
        /// <param name="iterations">Number of iterations.</param>
        /// <returns>Aggregated simulation results.</returns>
        public MonteCarloResult Run(IMonteCarloModel model, int iterations)
        {
            if (model == null) throw new ArgumentNullException(nameof(model));
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            var samples = new double[iterations];
            for (int i = 0; i < iterations; i++)
                samples[i] = model.Evaluate(_rng);

            return new MonteCarloResult(samples);
        }

        /// <summary>
        /// Runs a simulation using an inline function.
        /// The function receives the <see cref="RandomGenerator"/> and returns a scalar outcome.
        /// </summary>
        /// <param name="trialFunc">Function to execute per iteration.</param>
        /// <param name="iterations">Number of iterations.</param>
        /// <returns>Aggregated simulation results.</returns>
        public MonteCarloResult Run(Func<RandomGenerator, double> trialFunc, int iterations)
        {
            if (trialFunc == null) throw new ArgumentNullException(nameof(trialFunc));
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            var samples = new double[iterations];
            for (int i = 0; i < iterations; i++)
                samples[i] = trialFunc(_rng);

            return new MonteCarloResult(samples);
        }

        /// <summary>
        /// Monte Carlo integration: estimates E[f(X)] where X is drawn from the given distribution.
        /// Useful for computing integrals: ∫ f(x) p(x) dx ≈ (1/N) Σ f(xᵢ) where xᵢ ~ p(x).
        /// </summary>
        /// <param name="function">The function f(x) to integrate.</param>
        /// <param name="distribution">The distribution to sample X from.</param>
        /// <param name="iterations">Number of samples.</param>
        /// <returns>Aggregated results where Mean ≈ E[f(X)].</returns>
        public MonteCarloResult Integrate(Func<double, double> function, IDistribution distribution, int iterations)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            var samples = new double[iterations];
            for (int i = 0; i < iterations; i++)
            {
                double x = distribution.Sample(_rng);
                samples[i] = function(x);
            }

            return new MonteCarloResult(samples);
        }

        /// <summary>
        /// Monte Carlo integration over a uniform interval [a, b]:
        /// ∫ₐᵇ f(x) dx ≈ (b-a) * (1/N) Σ f(xᵢ), where xᵢ ~ U(a,b).
        /// The result Mean already contains the (b-a) scaling factor.
        /// </summary>
        /// <param name="function">The integrand f(x).</param>
        /// <param name="a">Lower integration bound.</param>
        /// <param name="b">Upper integration bound.</param>
        /// <param name="iterations">Number of samples.</param>
        /// <returns>Result where Mean ≈ ∫ₐᵇ f(x) dx.</returns>
        public MonteCarloResult IntegrateUniform(Func<double, double> function, double a, double b, int iterations)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (a >= b) throw new ArgumentException("a must be less than b.");
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            double scale = b - a;
            var samples = new double[iterations];
            for (int i = 0; i < iterations; i++)
            {
                double x = _rng.NextUniform(a, b);
                samples[i] = scale * function(x);
            }

            return new MonteCarloResult(samples);
        }

        /// <summary>
        /// Runs a multi-dimensional simulation where each trial produces a vector of outcomes.
        /// Returns one <see cref="MonteCarloResult"/> per output dimension.
        /// </summary>
        /// <param name="trialFunc">Function that returns an array of outputs per trial.</param>
        /// <param name="iterations">Number of iterations.</param>
        /// <returns>Array of results, one per output dimension.</returns>
        public MonteCarloResult[] RunMultivariate(Func<RandomGenerator, double[]> trialFunc, int iterations)
        {
            if (trialFunc == null) throw new ArgumentNullException(nameof(trialFunc));
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            // Run first iteration to determine dimensionality
            var first = trialFunc(_rng);
            int dims = first.Length;
            var allSamples = new double[dims][];
            for (int d = 0; d < dims; d++)
            {
                allSamples[d] = new double[iterations];
                allSamples[d][0] = first[d];
            }

            for (int i = 1; i < iterations; i++)
            {
                var result = trialFunc(_rng);
                for (int d = 0; d < dims; d++)
                    allSamples[d][i] = result[d];
            }

            var results = new MonteCarloResult[dims];
            for (int d = 0; d < dims; d++)
                results[d] = new MonteCarloResult(allSamples[d]);

            return results;
        }

        /// <summary>
        /// Runs convergence analysis: executes the trial function and records
        /// the running mean at each step. Useful for verifying that the simulation
        /// has converged.
        /// </summary>
        /// <param name="trialFunc">Function to execute per iteration.</param>
        /// <param name="iterations">Number of iterations.</param>
        /// <returns>Array of running means (index i = mean of first i+1 samples).</returns>
        public double[] ConvergenceCurve(Func<RandomGenerator, double> trialFunc, int iterations)
        {
            if (trialFunc == null) throw new ArgumentNullException(nameof(trialFunc));
            if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

            var runningMeans = new double[iterations];
            double sum = 0;

            for (int i = 0; i < iterations; i++)
            {
                sum += trialFunc(_rng);
                runningMeans[i] = sum / (i + 1);
            }

            return runningMeans;
        }
    }
}
