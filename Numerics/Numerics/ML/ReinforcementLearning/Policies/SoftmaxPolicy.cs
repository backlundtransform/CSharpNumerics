using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.ReinforcementLearning.Policies;

/// <summary>
/// Boltzmann (softmax) exploration policy.
/// Converts Q-values to action probabilities:
///   P(a) = exp(Q(a)/τ) / Σ exp(Q(a')/τ)
///
/// Temperature τ controls exploration:
///   - High τ → near-uniform (exploratory)
///   - Low τ → near-greedy (exploitative)
///
/// Supports exponential decay: τ ← max(τMin, τ × τDecay) after each Decay().
/// </summary>
public class SoftmaxPolicy : IPolicy
{
    public double Temperature { get; set; } = 1.0;
    public double TemperatureMin { get; set; } = 0.01;
    public double TemperatureDecay { get; set; } = 0.995;

    private readonly Random _rng;

    public SoftmaxPolicy(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    private SoftmaxPolicy(double temperature, double temperatureMin, double temperatureDecay, Random rng)
    {
        Temperature = temperature;
        TemperatureMin = temperatureMin;
        TemperatureDecay = temperatureDecay;
        _rng = rng;
    }

    public int SelectAction(VectorN qValues)
    {
        // Compute softmax probabilities with temperature scaling
        double maxQ = qValues[0];
        for (int i = 1; i < qValues.Length; i++)
            if (qValues[i] > maxQ) maxQ = qValues[i];

        var probs = new double[qValues.Length];
        double sum = 0;
        for (int i = 0; i < qValues.Length; i++)
        {
            probs[i] = Math.Exp((qValues[i] - maxQ) / Temperature);
            sum += probs[i];
        }
        for (int i = 0; i < probs.Length; i++)
            probs[i] /= sum;

        // Sample from the distribution
        double r = _rng.NextDouble();
        double cumulative = 0;
        for (int i = 0; i < probs.Length; i++)
        {
            cumulative += probs[i];
            if (r < cumulative) return i;
        }
        return probs.Length - 1;
    }

    public VectorN SelectAction(VectorN mean, VectorN std)
    {
        // Not applicable for Boltzmann; pass through mean
        return mean;
    }

    public void Decay()
    {
        Temperature = Math.Max(TemperatureMin, Temperature * TemperatureDecay);
    }

    public IPolicy Clone()
    {
        return new SoftmaxPolicy(Temperature, TemperatureMin, TemperatureDecay, new Random());
    }
}
