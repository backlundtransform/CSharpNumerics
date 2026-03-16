using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;

/// <summary>
/// On-policy TD(0) SARSA.
/// Update: Q(s,a) ← Q(s,a) + α [r + γ Q(s',a') − Q(s,a)]
/// where a' is the action actually taken in s' (on-policy).
///
/// Since the training loop provides (s,a,r,s',done) without a',
/// SARSA stores the previous transition and updates on the next step when a' is known.
/// </summary>
public class SARSA : TabularAgent
{
    public override string Name => "SARSA";

    private Transition _pending;
    private int _pendingNextAction;
    private bool _hasPending;

    public SARSA(int numStates, int numActions, Func<VectorN, int> stateMapper)
        : base(numStates, numActions, stateMapper) { }

    public override void Train(Transition t)
    {
        if (_hasPending)
        {
            // Now we know a' = t.Action (the action taken in the state that was s')
            int s = S(_pending.State);
            int a = _pending.Action;
            double target = _pending.Done
                ? _pending.Reward
                : _pending.Reward + Gamma * Q[S(_pending.NextState), t.Action];
            Q[s, a] += LearningRate * (target - Q[s, a]);
        }

        _pending = t;
        _hasPending = true;

        // If episode is done, do final update with Q(s',a') = 0
        if (t.Done)
        {
            int s = S(t.State);
            int a = t.Action;
            Q[s, a] += LearningRate * (t.Reward - Q[s, a]);
            _hasPending = false;
        }
    }

    public override void EndEpisode(Episode episode)
    {
        _hasPending = false;
        _pending = null;
    }

    public override IAgent Clone()
    {
        var clone = new SARSA(NumStates, NumActions, S)
        {
            LearningRate = LearningRate,
            Gamma = Gamma
        };
        Array.Copy(Q, clone.Q, Q.Length);
        return clone;
    }
}
