using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.ML.ReinforcementLearning.Policies;

/// <summary>
/// Ornstein-Uhlenbeck noise process for temporally correlated exploration.
/// Widely used with DDPG for smooth continuous control exploration.
///
/// dx_t = θ(μ - x_t) dt + σ dW_t
///
/// Parameters:
///   θ (Theta) — mean reversion speed (higher = decorrelates faster)
///   μ (Mu) — long-run mean (typically 0)
///   σ (Sigma) — noise volatility
///
/// The process tracks an internal state per action dimension.
/// Call Reset() at the start of each episode.
/// </summary>
public class OrnsteinUhlenbeck : IPolicy
{
    public double Theta { get; set; } = 0.15;
    public double Mu { get; set; } = 0.0;
    public double Sigma { get; set; } = 0.2;
    public double SigmaMin { get; set; } = 0.01;
    public double SigmaDecay { get; set; } = 1.0;
    public double Dt { get; set; } = 1.0;

    private double[] _state;
    private readonly Random _rng;
    private int _actionDim;

    public OrnsteinUhlenbeck(int? seed = null)
    {
        _rng = seed.HasValue ? new Random(seed.Value) : new Random();
    }

    private OrnsteinUhlenbeck(double theta, double mu, double sigma, double sigmaMin,
        double sigmaDecay, double dt, Random rng)
    {
        Theta = theta;
        Mu = mu;
        Sigma = sigma;
        SigmaMin = sigmaMin;
        SigmaDecay = sigmaDecay;
        Dt = dt;
        _rng = rng;
    }

    /// <summary>Reset the OU process state. Call at the start of each episode.</summary>
    public void Reset(int actionDim)
    {
        _actionDim = actionDim;
        _state = new double[actionDim];
        for (int i = 0; i < actionDim; i++)
            _state[i] = Mu;
    }

    /// <summary>Not applicable for OU — returns argmax of input.</summary>
    public int SelectAction(VectorN qValues)
    {
        int best = 0;
        for (int i = 1; i < qValues.Length; i++)
            if (qValues[i] > qValues[best]) best = i;
        return best;
    }

    /// <summary>
    /// Add OU noise to the mean action.
    /// Evolves the internal noise state and adds it to mean.
    /// </summary>
    public VectorN SelectAction(VectorN mean, VectorN std)
    {
        EnsureState(mean.Length);

        var noisy = new double[mean.Length];
        for (int i = 0; i < mean.Length; i++)
        {
            double dx = Theta * (Mu - _state[i]) * Dt + Sigma * Math.Sqrt(Dt) * SampleGaussian();
            _state[i] += dx;
            noisy[i] = mean[i] + _state[i];
        }
        return new VectorN(noisy);
    }

    public void Decay()
    {
        Sigma = Math.Max(SigmaMin, Sigma * SigmaDecay);
    }

    public IPolicy Clone()
    {
        return new OrnsteinUhlenbeck(Theta, Mu, Sigma, SigmaMin, SigmaDecay, Dt, new Random());
    }

    private void EnsureState(int actionDim)
    {
        if (_state == null || _state.Length != actionDim)
            Reset(actionDim);
    }

    private double SampleGaussian()
    {
        double u1 = 1.0 - _rng.NextDouble();
        double u2 = _rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
    }
}
