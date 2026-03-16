using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;

/// <summary>
/// Off-policy TD(0) Q-Learning.
/// Update: Q(s,a) ← Q(s,a) + α [r + γ max_a' Q(s',a') − Q(s,a)]
/// </summary>
public class QLearning : TabularAgent
{
    public override string Name => "QLearning";

    public QLearning(int numStates, int numActions, Func<VectorN, int> stateMapper)
        : base(numStates, numActions, stateMapper) { }

    public override void Train(Transition t)
    {
        int s = S(t.State);
        int a = t.Action;
        int sNext = S(t.NextState);

        double maxQ = double.NegativeInfinity;
        for (int i = 0; i < NumActions; i++)
            if (Q[sNext, i] > maxQ) maxQ = Q[sNext, i];

        double target = t.Done ? t.Reward : t.Reward + Gamma * maxQ;
        Q[s, a] += LearningRate * (target - Q[s, a]);
    }

    public override IAgent Clone()
    {
        var clone = new QLearning(NumStates, NumActions, S)
        {
            LearningRate = LearningRate,
            Gamma = Gamma
        };
        Array.Copy(Q, clone.Q, Q.Length);
        return clone;
    }
}
