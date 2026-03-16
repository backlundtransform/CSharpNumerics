using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Interfaces;

/// <summary>
/// Contract for an RL environment.
/// Discrete-action environments implement Step(int), continuous ones Step(VectorN).
/// </summary>
public interface IEnvironment
{
    /// <summary>Dimensionality of the observation/state vector.</summary>
    int ObservationSize { get; }

    /// <summary>Number of discrete actions, or dimensionality of continuous action space.</summary>
    int ActionSize { get; }

    /// <summary>True if the action space is discrete.</summary>
    bool IsDiscrete { get; }

    /// <summary>Reset to initial state. Returns (state, info).</summary>
    (VectorN state, Dictionary<string, object> info) Reset(int? seed = null);

    /// <summary>Take a discrete action. Returns (nextState, reward, done, info).</summary>
    (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(int action);

    /// <summary>Take a continuous action. Returns (nextState, reward, done, info).</summary>
    (VectorN nextState, double reward, bool done, Dictionary<string, object> info) Step(VectorN action);
}
