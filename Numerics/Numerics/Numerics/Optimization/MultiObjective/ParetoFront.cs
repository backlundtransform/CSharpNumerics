using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Numerics.Optimization.MultiObjective;

/// <summary>
/// Represents a single solution in a multi-objective optimisation.
/// Holds the decision variables and the evaluated objective values.
/// </summary>
public class ParetoSolution
{
    public double[] Variables { get; }
    public double[] Objectives { get; }
    public int Rank { get; internal set; }
    public double CrowdingDistance { get; internal set; }

    public ParetoSolution(double[] variables, double[] objectives)
    {
        Variables = (double[])variables.Clone();
        Objectives = (double[])objectives.Clone();
    }

    /// <summary>
    /// Returns true if this solution dominates <paramref name="other"/>
    /// (all objectives ≤ and at least one strictly &lt;).
    /// </summary>
    public bool Dominates(ParetoSolution other)
    {
        bool atLeastOneStrictlyBetter = false;
        for (int i = 0; i < Objectives.Length; i++)
        {
            if (Objectives[i] > other.Objectives[i]) return false;
            if (Objectives[i] < other.Objectives[i]) atLeastOneStrictlyBetter = true;
        }
        return atLeastOneStrictlyBetter;
    }
}

/// <summary>
/// Computes the Pareto front from a set of solutions using non-dominated sorting.
/// </summary>
public static class ParetoFront
{
    /// <summary>
    /// Partition <paramref name="solutions"/> into non-dominated fronts.
    /// Front 0 is the Pareto-optimal front.
    /// </summary>
    public static List<List<ParetoSolution>> NonDominatedSort(IReadOnlyList<ParetoSolution> solutions)
    {
        int n = solutions.Count;
        var dominationCount = new int[n];
        var dominated = new List<int>[n];

        for (int i = 0; i < n; i++)
            dominated[i] = new List<int>();

        var firstFront = new List<int>();

        for (int p = 0; p < n; p++)
        {
            for (int q = 0; q < n; q++)
            {
                if (p == q) continue;

                if (solutions[p].Dominates(solutions[q]))
                    dominated[p].Add(q);
                else if (solutions[q].Dominates(solutions[p]))
                    dominationCount[p]++;
            }

            if (dominationCount[p] == 0)
            {
                solutions[p].Rank = 0;
                firstFront.Add(p);
            }
        }

        var fronts = new List<List<ParetoSolution>>();
        var currentFrontIndices = firstFront;
        int rank = 0;

        while (currentFrontIndices.Count > 0)
        {
            var front = new List<ParetoSolution>();
            var nextFrontIndices = new List<int>();

            foreach (int p in currentFrontIndices)
            {
                solutions[p].Rank = rank;
                front.Add(solutions[p]);

                foreach (int q in dominated[p])
                {
                    dominationCount[q]--;
                    if (dominationCount[q] == 0)
                        nextFrontIndices.Add(q);
                }
            }

            fronts.Add(front);
            currentFrontIndices = nextFrontIndices;
            rank++;
        }

        return fronts;
    }

    /// <summary>
    /// Returns only the Pareto-optimal (rank-0) solutions.
    /// </summary>
    public static List<ParetoSolution> GetFront(IReadOnlyList<ParetoSolution> solutions)
    {
        var fronts = NonDominatedSort(solutions);
        return fronts.Count > 0 ? fronts[0] : new List<ParetoSolution>();
    }

    /// <summary>
    /// Compute crowding distance for a list of solutions (assumed to share the same rank).
    /// Higher crowding distance = more isolated = more desirable for diversity.
    /// </summary>
    public static void ComputeCrowdingDistance(IList<ParetoSolution> front)
    {
        int n = front.Count;
        if (n == 0) return;

        foreach (var s in front)
            s.CrowdingDistance = 0;

        int numObjectives = front[0].Objectives.Length;

        for (int m = 0; m < numObjectives; m++)
        {
            var sorted = front.OrderBy(s => s.Objectives[m]).ToList();

            sorted[0].CrowdingDistance = double.PositiveInfinity;
            sorted[n - 1].CrowdingDistance = double.PositiveInfinity;

            double range = sorted[n - 1].Objectives[m] - sorted[0].Objectives[m];
            if (range < 1e-15) continue;

            for (int i = 1; i < n - 1; i++)
            {
                sorted[i].CrowdingDistance +=
                    (sorted[i + 1].Objectives[m] - sorted[i - 1].Objectives[m]) / range;
            }
        }
    }
}
