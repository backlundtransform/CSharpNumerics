using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Buffers;

/// <summary>
/// Circular experience replay buffer with uniform random sampling.
/// </summary>
public class ReplayBuffer : IReplayBuffer
{
    private readonly Transition[] _buffer;
    private int _index;
    private int _count;
    private readonly Random _rng;

    public int Count => _count;
    public int Capacity { get; }

    public ReplayBuffer(int capacity, int? seed = null)
    {
        Capacity = capacity;
        _buffer = new Transition[capacity];
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    public void Add(Transition transition)
    {
        _buffer[_index] = transition;
        _index = (_index + 1) % Capacity;
        if (_count < Capacity) _count++;
    }

    public List<Transition> Sample(int batchSize)
    {
        if (batchSize > _count)
            throw new InvalidOperationException(
                $"Cannot sample {batchSize} from buffer with {_count} entries.");

        var batch = new List<Transition>(batchSize);
        for (int i = 0; i < batchSize; i++)
            batch.Add(_buffer[_rng.Next(_count)]);
        return batch;
    }
}
