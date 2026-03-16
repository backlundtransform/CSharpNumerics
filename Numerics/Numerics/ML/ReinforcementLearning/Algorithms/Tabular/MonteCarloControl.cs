using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;

/// <summary>
/// First-visit Monte Carlo control with exploring starts.
/// Updates Q-values using full episodic returns (no bootstrapping).
/// Q(s,a) ← Q(s,a) + α [G − Q(s,a)]   (incremental mean)
/// </summary>
public class MonteCarloControl : TabularAgent
{
    public override string Name => "MonteCarloControl";

    public MonteCarloControl(int numStates, int numActions, Func<VectorN, int> stateMapper)
        : base(numStates, numActions, stateMapper) { }

    /// <summary>
    /// Single-step Train is a no-op for MC — all learning happens in EndEpisode.
    /// </summary>
    public override void Train(Transition transition) { }

    public override void EndEpisode(Episode episode)
    {
        if (episode.Length == 0) return;

        double[] returns = episode.DiscountedReturns(Gamma);
        var visited = new HashSet<(int s, int a)>();

        for (int t = 0; t < episode.Length; t++)
        {
            int s = S(episode.Transitions[t].State);
            int a = episode.Transitions[t].Action;

            // First-visit: only update the first time (s,a) appears
            if (visited.Contains((s, a)))
                continue;
            visited.Add((s, a));

            Q[s, a] += LearningRate * (returns[t] - Q[s, a]);
        }
    }

    public override IAgent Clone()
    {
        var clone = new MonteCarloControl(NumStates, NumActions, S)
        {
            LearningRate = LearningRate,
            Gamma = Gamma
        };
        Array.Copy(Q, clone.Q, Q.Length);
        return clone;
    }
}
