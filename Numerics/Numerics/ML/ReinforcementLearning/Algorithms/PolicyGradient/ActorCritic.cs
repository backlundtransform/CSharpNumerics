using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;

/// <summary>
/// Advantage Actor-Critic (A2C).
///
/// Two-network architecture:
///   - Actor: π(a|s) — softmax policy network
///   - Critic: V(s) — state-value baseline network
///
/// Per-step TD advantage: A(s,a) = r + γ V(s') - V(s)
/// Actor gradient: ∇J ≈ Σ_t A_t ∇ log π(a_t|s_t)
/// Critic loss: MSE between V(s) and TD target (r + γ V(s'))
/// </summary>
public class ActorCritic : IAgent
{
    public string Name => "A2C";

    // ── Hyperparameters ─────────────────────────────────────
    public int[] ActorHiddenLayers { get; set; } = new[] { 32 };
    public int[] CriticHiddenLayers { get; set; } = new[] { 32 };
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double ActorLearningRate { get; set; } = 0.001;
    public double CriticLearningRate { get; set; } = 0.002;
    public double Gamma { get; set; } = 0.99;
    public double EntropyCoefficient { get; set; } = 0.01;

    // ── Internal state ──────────────────────────────────────
    private NN _actorNet;
    private NN _criticNet;
    private int _actionSize;
    private bool _initialized;
    private Random _rng;

    public void Initialize(int inputSize, int actionSize, int? seed = null)
    {
        _actionSize = actionSize;
        int s = seed ?? 123;
        _actorNet = new NN(inputSize, ActorHiddenLayers, actionSize, Activation, NN.OutputMode.Softmax, s);
        _criticNet = new NN(inputSize, CriticHiddenLayers, 1, Activation, NN.OutputMode.Linear, s + 1);
        _rng = new Random(s);
        _initialized = true;
    }

    // ── Action selection ────────────────────────────────────

    public int SelectAction(VectorN state)
    {
        EnsureInitialized();
        var probs = _actorNet.Forward(state);
        return SampleCategorical(probs);
    }

    /// <summary>Get action probabilities for a state.</summary>
    public VectorN GetActionProbabilities(VectorN state)
    {
        EnsureInitialized();
        return _actorNet.Forward(state);
    }

    /// <summary>Get the estimated state value V(s).</summary>
    public double GetValue(VectorN state)
    {
        EnsureInitialized();
        return _criticNet.Forward(state)[0];
    }

    public VectorN SelectContinuousAction(VectorN state) =>
        new VectorN(new double[] { SelectAction(state) });

    // ── Training (per-step TD actor-critic) ─────────────────

    /// <summary>
    /// Online single-step A2C update:
    ///   1. Compute TD advantage: A = r + γ V(s') - V(s)
    ///   2. Update critic toward TD target
    ///   3. Update actor with policy gradient scaled by advantage
    /// </summary>
    public void Train(Transition transition)
    {
        EnsureInitialized();

        // Critic forward for current and next state
        var vCurrent = _criticNet.Forward(transition.State, out var criticActivations);
        double v = vCurrent[0];
        double vNext = transition.Done ? 0.0 : _criticNet.Forward(transition.NextState)[0];

        double tdTarget = transition.Reward + Gamma * vNext;
        double advantage = tdTarget - v;

        // ── Update critic ───────────────────────────────────
        var criticTarget = new VectorN(new[] { tdTarget });
        _criticNet.ComputeGradients(criticTarget, criticActivations, out var cDw, out var cDb);
        _criticNet.ApplyGradients(cDw, cDb, CriticLearningRate, 1);

        // ── Update actor ────────────────────────────────────
        var probs = _actorNet.Forward(transition.State, out var actorActivations);

        var actorTarget = new double[_actionSize];
        for (int a = 0; a < _actionSize; a++)
        {
            double indicator = (a == transition.Action) ? 1.0 : 0.0;
            // Same as REINFORCE: push output toward action proportional to advantage
            actorTarget[a] = probs[a] - advantage * (indicator - probs[a]);

            // Entropy bonus: encourage exploration by pushing toward uniform
            // ∂H/∂π = -(1 + log π), gradient direction: push prob toward 1/|A|
            if (EntropyCoefficient > 0 && probs[a] > 1e-8)
            {
                double entropyGrad = -(1.0 + Math.Log(probs[a]));
                actorTarget[a] -= EntropyCoefficient * entropyGrad * probs[a] * (1.0 - probs[a]);
            }
        }

        _actorNet.ComputeGradients(new VectorN(actorTarget), actorActivations, out var aDw, out var aDb);
        _actorNet.ApplyGradients(aDw, aDb, ActorLearningRate, 1);
    }

    public void TrainBatch(List<Transition> batch)
    {
        foreach (var t in batch) Train(t);
    }

    public void EndEpisode(Episode episode) { }

    // ── Initialization ──────────────────────────────────────

    private void EnsureInitialized()
    {
        if (!_initialized)
            throw new InvalidOperationException(
                "ActorCritic must be initialized before use. Call Initialize(inputSize, actionSize) " +
                "or use RLExperiment which initializes automatically.");
    }

    // ── Clone / Hyperparameters ─────────────────────────────

    public IAgent Clone() => new ActorCritic
    {
        ActorHiddenLayers = (int[])ActorHiddenLayers.Clone(),
        CriticHiddenLayers = (int[])CriticHiddenLayers.Clone(),
        Activation = Activation,
        ActorLearningRate = ActorLearningRate,
        CriticLearningRate = CriticLearningRate,
        Gamma = Gamma,
        EntropyCoefficient = EntropyCoefficient
    };

    public Dictionary<string, object> GetHyperParameters() => new()
    {
        ["ActorHiddenLayers"] = ActorHiddenLayers,
        ["CriticHiddenLayers"] = CriticHiddenLayers,
        ["Activation"] = Activation,
        ["ActorLearningRate"] = ActorLearningRate,
        ["CriticLearningRate"] = CriticLearningRate,
        ["Gamma"] = Gamma,
        ["EntropyCoefficient"] = EntropyCoefficient
    };

    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("ActorHiddenLayers", out var ah)) ActorHiddenLayers = (int[])ah;
        if (p.TryGetValue("CriticHiddenLayers", out var ch)) CriticHiddenLayers = (int[])ch;
        if (p.TryGetValue("Activation", out var a)) Activation = (ActivationType)a;
        if (p.TryGetValue("ActorLearningRate", out var alr)) ActorLearningRate = Convert.ToDouble(alr);
        if (p.TryGetValue("CriticLearningRate", out var clr)) CriticLearningRate = Convert.ToDouble(clr);
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("EntropyCoefficient", out var ec)) EntropyCoefficient = Convert.ToDouble(ec);
    }

    // ── Helpers ─────────────────────────────────────────────

    private int SampleCategorical(VectorN probs)
    {
        double r = _rng.NextDouble();
        double cumulative = 0;
        for (int i = 0; i < probs.Length; i++)
        {
            cumulative += probs[i];
            if (r < cumulative) return i;
        }
        return probs.Length - 1;
    }
}
