using CSharpNumerics.ML.ReinforcementLearning.Core;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Interfaces;

/// <summary>
/// Experience replay buffer for off-policy algorithms.
/// </summary>
public interface IReplayBuffer
{
    void Add(Transition transition);
    List<Transition> Sample(int batchSize);
    int Count { get; }
    int Capacity { get; }
}
