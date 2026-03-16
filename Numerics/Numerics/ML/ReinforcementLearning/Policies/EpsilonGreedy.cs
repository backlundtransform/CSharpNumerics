using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.ReinforcementLearning.Policies;

/// <summary>
/// ε-greedy policy: with probability ε choose a random action, otherwise the greedy (argmax) action.
/// Supports exponential decay: ε ← max(εMin, ε × εDecay) after each call to Decay().
/// </summary>
public class EpsilonGreedy : IPolicy
{
    public double Epsilon { get; set; } = 1.0;
    public double EpsilonMin { get; set; } = 0.01;
    public double EpsilonDecay { get; set; } = 0.995;

    private readonly Random _rng;

    public EpsilonGreedy(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    private EpsilonGreedy(double epsilon, double epsilonMin, double epsilonDecay, Random rng)
    {
        Epsilon = epsilon;
        EpsilonMin = epsilonMin;
        EpsilonDecay = epsilonDecay;
        _rng = rng;
    }

    public int SelectAction(VectorN qValues)
    {
        if (_rng.NextDouble() < Epsilon)
            return _rng.Next(qValues.Length);

        // Greedy: argmax
        int best = 0;
        double bestVal = qValues[0];
        for (int i = 1; i < qValues.Length; i++)
        {
            if (qValues[i] > bestVal)
            {
                bestVal = qValues[i];
                best = i;
            }
        }
        return best;
    }

    public VectorN SelectAction(VectorN mean, VectorN std)
    {
        // Not applicable for ε-greedy; pass through mean
        return mean;
    }

    public void Decay()
    {
        Epsilon = Math.Max(EpsilonMin, Epsilon * EpsilonDecay);
    }

    public IPolicy Clone()
    {
        return new EpsilonGreedy(Epsilon, EpsilonMin, EpsilonDecay, new Random());
    }
}
