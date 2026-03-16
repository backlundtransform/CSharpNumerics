using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.ReinforcementLearning.Interfaces;

/// <summary>
/// Action-selection policy for exploration/exploitation.
/// </summary>
public interface IPolicy
{
    /// <summary>Select a discrete action from Q-values.</summary>
    int SelectAction(VectorN qValues);

    /// <summary>Select a continuous action from mean and std.</summary>
    VectorN SelectAction(VectorN mean, VectorN std);

    /// <summary>Per-episode decay (e.g. ε decay).</summary>
    void Decay();

    IPolicy Clone();
}
