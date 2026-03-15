using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Materials.Nuclear.DecayChains
{
    /// <summary>
    /// Models a sequential radioactive decay chain where a parent isotope
    /// decays through one or more daughter products.
    /// Uses the Bateman equations to compute mass/activity at any time t.
    /// <para>
    /// Example: Cs-137 → Ba-137m → Ba-137 (stable).
    /// </para>
    /// </summary>
    public class DecayChain
    {
        private readonly List<ChainStep> _steps = new List<ChainStep>();

        /// <summary>The parent isotope at the head of the chain.</summary>
        public Isotope Parent { get; }

        /// <summary>Number of decay steps in the chain.</summary>
        public int StepCount => _steps.Count;

        /// <summary>
        /// Creates a decay chain starting from the given parent isotope.
        /// </summary>
        public DecayChain(Isotope parent)
        {
            Parent = parent;
        }

        /// <summary>
        /// Adds a decay step: the last isotope in the chain decays to
        /// <paramref name="daughter"/> with the given branching ratio.
        /// </summary>
        /// <param name="daughter">The daughter isotope.</param>
        /// <param name="branchingRatio">Fraction of decays producing this daughter (0..1). Default 1.0.</param>
        public DecayChain AddStep(Isotope daughter, double branchingRatio = 1.0)
        {
            if (branchingRatio < 0 || branchingRatio > 1)
                throw new ArgumentOutOfRangeException(nameof(branchingRatio), "Must be between 0 and 1.");

            var parentIso = _steps.Count == 0 ? Parent : _steps[_steps.Count - 1].Daughter;
            _steps.Add(new ChainStep(parentIso, daughter, branchingRatio));
            return this;
        }

        /// <summary>
        /// Returns all isotopes in the chain (parent + all daughters).
        /// </summary>
        public Isotope[] AllIsotopes()
        {
            var result = new Isotope[_steps.Count + 1];
            result[0] = Parent;
            for (int i = 0; i < _steps.Count; i++)
                result[i + 1] = _steps[i].Daughter;
            return result;
        }

        /// <summary>
        /// Evolves the decay chain from the given initial masses over <paramref name="timeSeconds"/>.
        /// Returns the mass (kg) of each isotope at time t using the Bateman equations.
        /// </summary>
        /// <param name="initialMasses">
        /// Initial mass (kg) for each isotope in order: [parent, daughter1, daughter2, ...].
        /// Length must equal <see cref="StepCount"/> + 1.
        /// </param>
        /// <param name="timeSeconds">Elapsed time (seconds).</param>
        /// <returns>Array of masses at time t, same order as <see cref="AllIsotopes"/>.</returns>
        public double[] Evolve(double[] initialMasses, double timeSeconds)
        {
            int n = _steps.Count + 1;
            if (initialMasses == null) throw new ArgumentNullException(nameof(initialMasses));
            if (initialMasses.Length != n)
                throw new ArgumentException($"Expected {n} initial masses, got {initialMasses.Length}.");

            var isotopes = AllIsotopes();
            var lambdas = new double[n];
            var branchingRatios = new double[n]; // branching ratio for step i-1 → i
            for (int i = 0; i < n; i++)
            {
                lambdas[i] = isotopes[i].Lambda;
                branchingRatios[i] = i > 0 ? _steps[i - 1].BranchingRatio : 1.0;
            }

            var result = new double[n];

            // For each isotope j, the Bateman solution for contribution from
            // ancestor i (where i ≤ j) through the chain i → i+1 → ... → j:
            //
            //   N_j(t) = Σ_{i=0..j} initialN_i × BatemanCoeff(i→j, t)
            //
            // BatemanCoeff(i→j, t) = [Π_{k=i..j-1} λ_k × BR_{k+1}] × Σ_{m=i..j} exp(-λ_m t) / Π_{p=i..j, p≠m}(λ_p - λ_m)

            for (int j = 0; j < n; j++)
            {
                double nj = 0;

                for (int i = 0; i <= j; i++)
                {
                    if (initialMasses[i] <= 0) continue;

                    double coeff = BatemanCoefficient(lambdas, branchingRatios, i, j, timeSeconds);
                    nj += initialMasses[i] * coeff;
                }

                result[j] = Math.Max(0, nj);
            }

            return result;
        }

        /// <summary>
        /// Computes the activity (Bq) of each isotope at time t.
        /// </summary>
        public double[] EvolveActivity(double[] initialMasses, double timeSeconds)
        {
            var masses = Evolve(initialMasses, timeSeconds);
            var isotopes = AllIsotopes();
            var activities = new double[masses.Length];
            for (int i = 0; i < masses.Length; i++)
                activities[i] = Decay.Decay.Activity(masses[i], isotopes[i]);
            return activities;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Bateman equation solver
        // ═══════════════════════════════════════════════════════════════

        private static double BatemanCoefficient(
            double[] lambdas, double[] branchingRatios,
            int from, int to, double t)
        {
            if (from == to)
            {
                // Simple exponential decay of the isotope itself
                if (lambdas[from] == 0) return 1.0; // stable
                return Math.Exp(-lambdas[from] * t);
            }

            // Check if any intermediate isotope is stable (blocks the chain)
            for (int k = from; k < to; k++)
            {
                if (lambdas[k] == 0) return 0; // stable isotope blocks chain
            }

            // Product of λ_k × BR_{k+1} for k = from to to-1
            double product = 1.0;
            for (int k = from; k < to; k++)
                product *= lambdas[k] * branchingRatios[k + 1];

            // Summation over m = from..to
            double sum = 0;
            for (int m = from; m <= to; m++)
            {
                if (lambdas[m] == 0)
                {
                    // Stable endpoint: exp(0) = 1
                    double denom = 1.0;
                    for (int p = from; p <= to; p++)
                    {
                        if (p == m) continue;
                        denom *= (lambdas[p] - lambdas[m]);
                    }
                    if (Math.Abs(denom) < 1e-30) continue;
                    sum += 1.0 / denom;
                }
                else
                {
                    double denom = 1.0;
                    bool degenerate = false;
                    for (int p = from; p <= to; p++)
                    {
                        if (p == m) continue;
                        double diff = lambdas[p] - lambdas[m];
                        if (Math.Abs(diff) < 1e-30 * Math.Max(Math.Abs(lambdas[p]), Math.Abs(lambdas[m])))
                        {
                            degenerate = true;
                            break;
                        }
                        denom *= diff;
                    }
                    if (degenerate) continue;
                    sum += Math.Exp(-lambdas[m] * t) / denom;
                }
            }

            return product * sum;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Common chains
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Creates the Cs-137 decay chain: Cs-137 → Ba-137m (94.7%) → Ba-137 (stable).
        /// </summary>
        public static DecayChain Cs137Chain() =>
            new DecayChain(Isotope.Cs137)
                .AddStep(Isotope.Ba137m, 0.947)
                .AddStep(Isotope.Ba137);

        /// <summary>
        /// Creates the I-131 decay chain: I-131 → Xe-131 (stable).
        /// </summary>
        public static DecayChain I131Chain() =>
            new DecayChain(Isotope.I131)
                .AddStep(Isotope.Xe131);

        /// <summary>Helper struct for a decay step.</summary>
        private readonly struct ChainStep
        {
            public readonly Isotope Parent;
            public readonly Isotope Daughter;
            public readonly double BranchingRatio;

            public ChainStep(Isotope parent, Isotope daughter, double branchingRatio)
            {
                Parent = parent;
                Daughter = daughter;
                BranchingRatio = branchingRatio;
            }
        }
    }
}
