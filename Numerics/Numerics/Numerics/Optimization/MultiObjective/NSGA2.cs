using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Numerics.Optimization.MultiObjective;

/// <summary>
/// NSGA-II (Non-dominated Sorting Genetic Algorithm II) by Deb et al. (2002).
/// A population-based multi-objective optimiser that maintains diversity
/// via crowding distance and converges towards the Pareto front.
/// All objectives are minimised.
/// </summary>
public class NSGA2
{
    public int PopulationSize { get; set; }
    public int Generations { get; set; }
    public double CrossoverRate { get; set; }
    public double MutationRate { get; set; }
    public double MutationScale { get; set; }
    public int Seed { get; set; }

    private readonly Func<double[], double[]> _evaluate;
    private readonly int _numVariables;
    private readonly double[] _lowerBounds;
    private readonly double[] _upperBounds;

    /// <summary>
    /// Create a new NSGA-II optimiser.
    /// </summary>
    /// <param name="evaluate">Function mapping decision variables → objective values (all minimised).</param>
    /// <param name="numVariables">Number of decision variables.</param>
    /// <param name="lowerBounds">Lower bounds per variable.</param>
    /// <param name="upperBounds">Upper bounds per variable.</param>
    /// <param name="populationSize">Population size (should be even).</param>
    /// <param name="generations">Number of generations to evolve.</param>
    /// <param name="crossoverRate">Probability of crossover per pair.</param>
    /// <param name="mutationRate">Probability of mutation per variable.</param>
    /// <param name="mutationScale">Scale of Gaussian mutation relative to variable range.</param>
    /// <param name="seed">Random seed for reproducibility.</param>
    public NSGA2(
        Func<double[], double[]> evaluate,
        int numVariables,
        double[] lowerBounds,
        double[] upperBounds,
        int populationSize = 100,
        int generations = 200,
        double crossoverRate = 0.9,
        double mutationRate = 0.1,
        double mutationScale = 0.1,
        int seed = 42)
    {
        _evaluate = evaluate;
        _numVariables = numVariables;
        _lowerBounds = (double[])lowerBounds.Clone();
        _upperBounds = (double[])upperBounds.Clone();
        PopulationSize = populationSize % 2 == 0 ? populationSize : populationSize + 1;
        Generations = generations;
        CrossoverRate = crossoverRate;
        MutationRate = mutationRate;
        MutationScale = mutationScale;
        Seed = seed;
    }

    /// <summary>
    /// Run the NSGA-II algorithm and return the final Pareto front.
    /// </summary>
    public NSGA2Result Run()
    {
        var rnd = new Random(Seed);

        // Initialise random population
        var population = new List<ParetoSolution>(PopulationSize);
        for (int i = 0; i < PopulationSize; i++)
        {
            var vars = RandomIndividual(rnd);
            population.Add(new ParetoSolution(vars, _evaluate(vars)));
        }

        for (int gen = 0; gen < Generations; gen++)
        {
            // Create offspring via crossover + mutation
            var offspring = new List<ParetoSolution>(PopulationSize);
            while (offspring.Count < PopulationSize)
            {
                var parent1 = TournamentSelect(population, rnd);
                var parent2 = TournamentSelect(population, rnd);

                double[] child1, child2;
                if (rnd.NextDouble() < CrossoverRate)
                    SBXCrossover(parent1.Variables, parent2.Variables, out child1, out child2, rnd);
                else
                {
                    child1 = (double[])parent1.Variables.Clone();
                    child2 = (double[])parent2.Variables.Clone();
                }

                PolynomialMutation(child1, rnd);
                PolynomialMutation(child2, rnd);

                Clamp(child1);
                Clamp(child2);

                offspring.Add(new ParetoSolution(child1, _evaluate(child1)));
                offspring.Add(new ParetoSolution(child2, _evaluate(child2)));
            }

            // Combine parent + offspring
            var combined = new List<ParetoSolution>(population.Count + offspring.Count);
            combined.AddRange(population);
            combined.AddRange(offspring);

            // Non-dominated sorting
            var fronts = ParetoFront.NonDominatedSort(combined);

            // Select next generation
            population = new List<ParetoSolution>(PopulationSize);
            foreach (var front in fronts)
            {
                ParetoFront.ComputeCrowdingDistance(front);

                if (population.Count + front.Count <= PopulationSize)
                {
                    population.AddRange(front);
                }
                else
                {
                    int remaining = PopulationSize - population.Count;
                    var sorted = front.OrderByDescending(s => s.CrowdingDistance).Take(remaining);
                    population.AddRange(sorted);
                    break;
                }
            }
        }

        var finalFronts = ParetoFront.NonDominatedSort(population);
        var paretoOptimal = finalFronts.Count > 0 ? finalFronts[0] : new List<ParetoSolution>();
        ParetoFront.ComputeCrowdingDistance(paretoOptimal);

        return new NSGA2Result(population, paretoOptimal);
    }

    // ── Operators ───────────────────────────────────────────────

    private double[] RandomIndividual(Random rnd)
    {
        var vars = new double[_numVariables];
        for (int i = 0; i < _numVariables; i++)
            vars[i] = _lowerBounds[i] + rnd.NextDouble() * (_upperBounds[i] - _lowerBounds[i]);
        return vars;
    }

    private ParetoSolution TournamentSelect(List<ParetoSolution> pop, Random rnd)
    {
        var a = pop[rnd.Next(pop.Count)];
        var b = pop[rnd.Next(pop.Count)];

        if (a.Rank < b.Rank) return a;
        if (b.Rank < a.Rank) return b;
        return a.CrowdingDistance >= b.CrowdingDistance ? a : b;
    }

    /// <summary>Simulated binary crossover (SBX).</summary>
    private void SBXCrossover(double[] p1, double[] p2,
        out double[] c1, out double[] c2, Random rnd, double eta = 20)
    {
        c1 = new double[_numVariables];
        c2 = new double[_numVariables];

        for (int i = 0; i < _numVariables; i++)
        {
            if (rnd.NextDouble() < 0.5 && Math.Abs(p1[i] - p2[i]) > 1e-14)
            {
                double u = rnd.NextDouble();
                double beta = u <= 0.5
                    ? Math.Pow(2.0 * u, 1.0 / (eta + 1))
                    : Math.Pow(1.0 / (2.0 * (1 - u)), 1.0 / (eta + 1));

                c1[i] = 0.5 * ((1 + beta) * p1[i] + (1 - beta) * p2[i]);
                c2[i] = 0.5 * ((1 - beta) * p1[i] + (1 + beta) * p2[i]);
            }
            else
            {
                c1[i] = p1[i];
                c2[i] = p2[i];
            }
        }
    }

    /// <summary>Polynomial mutation.</summary>
    private void PolynomialMutation(double[] individual, Random rnd, double eta = 20)
    {
        for (int i = 0; i < _numVariables; i++)
        {
            if (rnd.NextDouble() < MutationRate)
            {
                double range = _upperBounds[i] - _lowerBounds[i];
                double u = rnd.NextDouble();
                double delta = u < 0.5
                    ? Math.Pow(2.0 * u, 1.0 / (eta + 1)) - 1.0
                    : 1.0 - Math.Pow(2.0 * (1 - u), 1.0 / (eta + 1));

                individual[i] += delta * range * MutationScale;
            }
        }
    }

    private void Clamp(double[] vars)
    {
        for (int i = 0; i < _numVariables; i++)
            vars[i] = Math.Max(_lowerBounds[i], Math.Min(_upperBounds[i], vars[i]));
    }
}

/// <summary>
/// Result of an NSGA-II run.
/// </summary>
public class NSGA2Result
{
    public IReadOnlyList<ParetoSolution> FinalPopulation { get; }
    public IReadOnlyList<ParetoSolution> ParetoFront { get; }

    public NSGA2Result(List<ParetoSolution> population, List<ParetoSolution> paretoFront)
    {
        FinalPopulation = population;
        ParetoFront = paretoFront;
    }

    /// <summary>Number of Pareto-optimal solutions.</summary>
    public int FrontSize => ParetoFront.Count;
}
