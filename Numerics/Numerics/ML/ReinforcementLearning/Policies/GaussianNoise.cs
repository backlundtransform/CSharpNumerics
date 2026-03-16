using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.ReinforcementLearning.Policies;

/// <summary>
/// Additive Gaussian noise exploration for continuous action spaces.
/// Adds i.i.d. N(0, σ²) noise to each action dimension.
/// Supports decay: σ ← max(σMin, σ × σDecay) after each Decay().
/// </summary>
public class GaussianNoise : IPolicy
{
    public double Sigma { get; set; } = 0.1;
    public double SigmaMin { get; set; } = 0.01;
    public double SigmaDecay { get; set; } = 0.999;

    private readonly Random _rng;

    public GaussianNoise(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    private GaussianNoise(double sigma, double sigmaMin, double sigmaDecay, Random rng)
    {
        Sigma = sigma;
        SigmaMin = sigmaMin;
        SigmaDecay = sigmaDecay;
        _rng = rng;
    }

    /// <summary>Not applicable for continuous noise — returns argmax of input.</summary>
    public int SelectAction(VectorN qValues)
    {
        int best = 0;
        for (int i = 1; i < qValues.Length; i++)
            if (qValues[i] > qValues[best]) best = i;
        return best;
    }

    /// <summary>
    /// Add Gaussian noise to the mean action.
    /// std parameter is ignored — uses Sigma instead.
    /// </summary>
    public VectorN SelectAction(VectorN mean, VectorN std)
    {
        var noisy = new double[mean.Length];
        for (int i = 0; i < mean.Length; i++)
            noisy[i] = mean[i] + Sigma * SampleGaussian();
        return new VectorN(noisy);
    }

    public void Decay()
    {
        Sigma = Math.Max(SigmaMin, Sigma * SigmaDecay);
    }

    public IPolicy Clone()
    {
        return new GaussianNoise(Sigma, SigmaMin, SigmaDecay, new Random());
    }

    private double SampleGaussian()
    {
        // Box-Muller transform
        double u1 = 1.0 - _rng.NextDouble();
        double u2 = _rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
    }
}
