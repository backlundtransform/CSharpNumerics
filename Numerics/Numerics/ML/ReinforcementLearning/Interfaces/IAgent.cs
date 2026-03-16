using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Interfaces;

/// <summary>
/// Contract for an RL agent that can select actions and learn from experience.
/// </summary>
public interface IAgent
{
    string Name { get; }

    /// <summary>Select a discrete action given the current state.</summary>
    int SelectAction(VectorN state);

    /// <summary>Select a continuous action given the current state.</summary>
    VectorN SelectContinuousAction(VectorN state);

    /// <summary>Online single-step update from one transition.</summary>
    void Train(Transition transition);

    /// <summary>Batch update from multiple transitions.</summary>
    void TrainBatch(List<Transition> batch);

    /// <summary>Post-episode hook (e.g. Monte Carlo return computation).</summary>
    void EndEpisode(Episode episode);

    IAgent Clone();
    Dictionary<string, object> GetHyperParameters();
    void SetHyperParameters(Dictionary<string, object> parameters);
}
