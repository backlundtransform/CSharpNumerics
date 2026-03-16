using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.ReinforcementLearning.Core;

/// <summary>
/// A single (s, a, r, s', done) experience tuple.
/// Supports both discrete actions (int) and continuous actions (VectorN).
/// </summary>
public class Transition
{
    public VectorN State { get; }
    public int Action { get; }
    public VectorN ContinuousAction { get; }
    public bool IsContinuous { get; }
    public double Reward { get; }
    public VectorN NextState { get; }
    public bool Done { get; }

    /// <summary>Create a transition with a discrete action.</summary>
    public Transition(VectorN state, int action, double reward, VectorN nextState, bool done)
    {
        State = state;
        Action = action;
        IsContinuous = false;
        Reward = reward;
        NextState = nextState;
        Done = done;
    }

    /// <summary>Create a transition with a continuous action.</summary>
    public Transition(VectorN state, VectorN continuousAction, double reward, VectorN nextState, bool done)
    {
        State = state;
        ContinuousAction = continuousAction;
        IsContinuous = true;
        Action = -1;
        Reward = reward;
        NextState = nextState;
        Done = done;
    }
}
