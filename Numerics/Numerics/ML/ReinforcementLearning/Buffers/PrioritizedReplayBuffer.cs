using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Buffers;

/// <summary>
/// Prioritized experience replay buffer (Schaul et al. 2016).
/// Samples transitions with probability proportional to their TD error priority.
/// Uses a simple array-based proportional priority scheme.
///
/// Priority: p_i = |δ_i| + ε   where δ is the TD error.
/// Probability: P(i) = p_i^α / Σ_j p_j^α
/// </summary>
public class PrioritizedReplayBuffer : IReplayBuffer
{
    private readonly Transition[] _buffer;
    private readonly double[] _priorities;
    private int _index;
    private int _count;
    private readonly Random _rng;

    public double Alpha { get; set; } = 0.6;
    public double Epsilon { get; set; } = 1e-5;
    public double DefaultPriority { get; set; } = 1.0;

    public int Count => _count;
    public int Capacity { get; }

    public PrioritizedReplayBuffer(int capacity, int? seed = null)
    {
        Capacity = capacity;
        _buffer = new Transition[capacity];
        _priorities = new double[capacity];
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    public void Add(Transition transition)
    {
        _buffer[_index] = transition;
        _priorities[_index] = DefaultPriority; // new transitions get max priority
        _index = (_index + 1) % Capacity;
        if (_count < Capacity) _count++;
    }

    /// <summary>Update the priority of a specific buffer index.</summary>
    public void UpdatePriority(int bufferIndex, double tdError)
    {
        _priorities[bufferIndex] = Math.Pow(Math.Abs(tdError) + Epsilon, Alpha);
    }

    public List<Transition> Sample(int batchSize)
    {
        if (batchSize > _count)
            throw new InvalidOperationException(
                $"Cannot sample {batchSize} from buffer with {_count} entries.");

        // Compute sampling probabilities
        double totalPriority = 0;
        for (int i = 0; i < _count; i++)
            totalPriority += _priorities[i];

        var batch = new List<Transition>(batchSize);
        for (int b = 0; b < batchSize; b++)
        {
            double rand = _rng.NextDouble() * totalPriority;
            double cumulative = 0;
            for (int i = 0; i < _count; i++)
            {
                cumulative += _priorities[i];
                if (cumulative >= rand)
                {
                    batch.Add(_buffer[i]);
                    break;
                }
            }
            // Edge case: floating point rounding — add last element
            if (batch.Count <= b)
                batch.Add(_buffer[_count - 1]);
        }

        return batch;
    }

    /// <summary>
    /// Sample with indices (so caller can update priorities after computing TD errors).
    /// Returns (transitions, bufferIndices).
    /// </summary>
    public (List<Transition> transitions, List<int> indices) SampleWithIndices(int batchSize)
    {
        if (batchSize > _count)
            throw new InvalidOperationException(
                $"Cannot sample {batchSize} from buffer with {_count} entries.");

        double totalPriority = 0;
        for (int i = 0; i < _count; i++)
            totalPriority += _priorities[i];

        var transitions = new List<Transition>(batchSize);
        var indices = new List<int>(batchSize);

        for (int b = 0; b < batchSize; b++)
        {
            double rand = _rng.NextDouble() * totalPriority;
            double cumulative = 0;
            int selected = _count - 1;
            for (int i = 0; i < _count; i++)
            {
                cumulative += _priorities[i];
                if (cumulative >= rand)
                {
                    selected = i;
                    break;
                }
            }
            transitions.Add(_buffer[selected]);
            indices.Add(selected);
        }

        return (transitions, indices);
    }
}
