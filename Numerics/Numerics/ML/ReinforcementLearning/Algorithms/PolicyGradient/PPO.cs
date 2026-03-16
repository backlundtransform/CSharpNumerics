using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;

/// <summary>
/// Proximal Policy Optimization (Schulman et al. 2017) — clipped surrogate objective.
///
/// Collects a full episode, then performs multiple epochs of mini-batch updates:
///   L_clip = E[min(r_t(θ) A_t, clip(r_t, 1-ε, 1+ε) A_t)]
///
/// where r_t(θ) = π_new(a_t|s_t) / π_old(a_t|s_t) is the probability ratio.
///
/// Uses GAE (Generalized Advantage Estimation) for advantage computation:
///   A_t = Σ_{l=0}^{T-t} (γλ)^l δ_{t+l}
///   δ_t = r_t + γ V(s_{t+1}) - V(s_t)
/// </summary>
public class PPO : IAgent
{
    public string Name => "PPO";

    // ── Hyperparameters ─────────────────────────────────────
    public int[] ActorHiddenLayers { get; set; } = new[] { 64, 64 };
    public int[] CriticHiddenLayers { get; set; } = new[] { 64, 64 };
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double ActorLearningRate { get; set; } = 0.0003;
    public double CriticLearningRate { get; set; } = 0.001;
    public double Gamma { get; set; } = 0.99;
    public double Lambda { get; set; } = 0.95;
    public double ClipEpsilon { get; set; } = 0.2;
    public int UpdateEpochs { get; set; } = 4;
    public int MiniBatchSize { get; set; } = 64;
    public double EntropyCoefficient { get; set; } = 0.01;

    // ── Internal state ──────────────────────────────────────
    private NN _actorNet;
    private NN _criticNet;
    private int _actionSize;
    private bool _initialized;
    private Random _rng;

    // Episode trajectory
    private List<VectorN> _states;
    private List<int> _actions;
    private List<double> _rewards;
    private List<VectorN> _nextStates;
    private List<bool> _dones;
    private List<double> _oldLogProbs;

    public void Initialize(int inputSize, int actionSize, int? seed = null)
    {
        _actionSize = actionSize;
        int s = seed ?? 123;
        _actorNet = new NN(inputSize, ActorHiddenLayers, actionSize, Activation, NN.OutputMode.Softmax, s);
        _criticNet = new NN(inputSize, CriticHiddenLayers, 1, Activation, NN.OutputMode.Linear, s + 1);
        _rng = new Random(s);
        ClearTrajectory();
        _initialized = true;
    }

    // ── Action selection ────────────────────────────────────

    public int SelectAction(VectorN state)
    {
        EnsureInitialized();
        var probs = _actorNet.Forward(state);
        int action = SampleCategorical(probs);

        // Store trajectory data
        double logProb = Math.Log(Math.Max(probs[action], 1e-10));
        _states.Add(state);
        _actions.Add(action);
        _oldLogProbs.Add(logProb);

        return action;
    }

    /// <summary>Get action probabilities for a state.</summary>
    public VectorN GetActionProbabilities(VectorN state)
    {
        EnsureInitialized();
        return _actorNet.Forward(state);
    }

    /// <summary>Get the critic's value estimate V(s) for a state.</summary>
    public double GetValue(VectorN state)
    {
        EnsureInitialized();
        return _criticNet.Forward(state)[0];
    }

    public VectorN SelectContinuousAction(VectorN state) =>
        new VectorN(new double[] { SelectAction(state) });

    // ── Training ────────────────────────────────────────────

    /// <summary>Store transition data — actual optimization happens in EndEpisode.</summary>
    public void Train(Transition transition)
    {
        EnsureInitialized();
        _rewards.Add(transition.Reward);
        _nextStates.Add(transition.NextState);
        _dones.Add(transition.Done);
    }

    public void TrainBatch(List<Transition> batch)
    {
        foreach (var t in batch) Train(t);
    }

    /// <summary>
    /// Perform PPO update at end of episode:
    ///   1. Compute GAE advantages
    ///   2. Multiple epochs of mini-batch clipped surrogate updates
    /// </summary>
    public void EndEpisode(Episode episode)
    {
        EnsureInitialized();
        int T = _states.Count;
        if (T == 0) return;

        // ── Compute GAE advantages ──────────────────────────
        var advantages = new double[T];
        var returns = new double[T];

        double lastGae = 0;
        for (int t = T - 1; t >= 0; t--)
        {
            double vCurrent = _criticNet.Forward(_states[t])[0];
            double vNext = _dones[t] ? 0.0 : _criticNet.Forward(_nextStates[t])[0];
            double delta = _rewards[t] + Gamma * vNext - vCurrent;
            double doneMultiplier = _dones[t] ? 0.0 : 1.0;
            lastGae = delta + Gamma * Lambda * doneMultiplier * lastGae;
            advantages[t] = lastGae;
            returns[t] = advantages[t] + vCurrent;
        }

        // Normalize advantages
        double mean = advantages.Average();
        double std = Math.Sqrt(advantages.Select(a => (a - mean) * (a - mean)).Average() + 1e-8);
        for (int t = 0; t < T; t++)
            advantages[t] = (advantages[t] - mean) / std;

        // ── Mini-batch PPO updates ──────────────────────────
        var indices = Enumerable.Range(0, T).ToArray();

        for (int epoch = 0; epoch < UpdateEpochs; epoch++)
        {
            Shuffle(indices);

            for (int start = 0; start < T; start += MiniBatchSize)
            {
                int end = Math.Min(start + MiniBatchSize, T);
                UpdateMiniBatch(indices, start, end, advantages, returns);
            }
        }

        ClearTrajectory();
    }

    private void UpdateMiniBatch(int[] indices, int start, int end, double[] advantages, double[] returns)
    {
        var actorWeightGrads = _actorNet.InitWeightGrads();
        var actorBiasGrads = _actorNet.InitBiasGrads();
        var criticWeightGrads = _criticNet.InitWeightGrads();
        var criticBiasGrads = _criticNet.InitBiasGrads();

        int batchSize = end - start;

        for (int b = start; b < end; b++)
        {
            int t = indices[b];

            // ── Actor update (clipped surrogate) ────────────
            var probs = _actorNet.Forward(_states[t], out var actorActivations);
            double newLogProb = Math.Log(Math.Max(probs[_actions[t]], 1e-10));
            double ratio = Math.Exp(newLogProb - _oldLogProbs[t]);

            double adv = advantages[t];
            double surr1 = ratio * adv;
            double surr2 = Math.Clamp(ratio, 1.0 - ClipEpsilon, 1.0 + ClipEpsilon) * adv;
            double ppoObjective = Math.Min(surr1, surr2);

            // Gradient direction: push action prob by clipped advantage
            var actorTarget = new double[_actionSize];
            for (int a = 0; a < _actionSize; a++)
            {
                double indicator = (a == _actions[t]) ? 1.0 : 0.0;
                // Use the effective advantage (clipped)
                double effectiveAdv = (Math.Abs(ratio - 1.0) <= ClipEpsilon || surr1 <= surr2)
                    ? adv : 0.0; // clip kills gradient when ratio strays too far

                actorTarget[a] = probs[a] - effectiveAdv * (indicator - probs[a]);

                // Entropy bonus
                if (EntropyCoefficient > 0 && probs[a] > 1e-8)
                {
                    double entropyGrad = -(1.0 + Math.Log(probs[a]));
                    actorTarget[a] -= EntropyCoefficient * entropyGrad * probs[a] * (1.0 - probs[a]);
                }
            }

            _actorNet.ComputeGradients(new VectorN(actorTarget), actorActivations, out var aDw, out var aDb);
            _actorNet.AccumulateGradients(actorWeightGrads, actorBiasGrads, aDw, aDb);

            // ── Critic update (MSE toward return) ───────────
            var vPred = _criticNet.Forward(_states[t], out var criticActivations);
            var criticTarget = new VectorN(new[] { returns[t] });
            _criticNet.ComputeGradients(criticTarget, criticActivations, out var cDw, out var cDb);
            _criticNet.AccumulateGradients(criticWeightGrads, criticBiasGrads, cDw, cDb);
        }

        _actorNet.ApplyGradients(actorWeightGrads, actorBiasGrads, ActorLearningRate, batchSize);
        _criticNet.ApplyGradients(criticWeightGrads, criticBiasGrads, CriticLearningRate, batchSize);
    }

    // ── Initialization ──────────────────────────────────────

    private void ClearTrajectory()
    {
        _states = new List<VectorN>();
        _actions = new List<int>();
        _rewards = new List<double>();
        _nextStates = new List<VectorN>();
        _dones = new List<bool>();
        _oldLogProbs = new List<double>();
    }

    private void EnsureInitialized()
    {
        if (!_initialized)
            throw new InvalidOperationException(
                "PPO must be initialized before use. Call Initialize(inputSize, actionSize) " +
                "or use RLExperiment which initializes automatically.");
    }

    // ── Clone / Hyperparameters ─────────────────────────────

    public IAgent Clone() => new PPO
    {
        ActorHiddenLayers = (int[])ActorHiddenLayers.Clone(),
        CriticHiddenLayers = (int[])CriticHiddenLayers.Clone(),
        Activation = Activation,
        ActorLearningRate = ActorLearningRate,
        CriticLearningRate = CriticLearningRate,
        Gamma = Gamma,
        Lambda = Lambda,
        ClipEpsilon = ClipEpsilon,
        UpdateEpochs = UpdateEpochs,
        MiniBatchSize = MiniBatchSize,
        EntropyCoefficient = EntropyCoefficient
    };

    public Dictionary<string, object> GetHyperParameters() => new()
    {
        ["ActorHiddenLayers"] = ActorHiddenLayers,
        ["CriticHiddenLayers"] = CriticHiddenLayers,
        ["ActorLearningRate"] = ActorLearningRate,
        ["CriticLearningRate"] = CriticLearningRate,
        ["Gamma"] = Gamma,
        ["Lambda"] = Lambda,
        ["ClipEpsilon"] = ClipEpsilon,
        ["UpdateEpochs"] = UpdateEpochs,
        ["MiniBatchSize"] = MiniBatchSize,
        ["EntropyCoefficient"] = EntropyCoefficient
    };

    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("ActorHiddenLayers", out var ah)) ActorHiddenLayers = (int[])ah;
        if (p.TryGetValue("CriticHiddenLayers", out var ch)) CriticHiddenLayers = (int[])ch;
        if (p.TryGetValue("ActorLearningRate", out var alr)) ActorLearningRate = Convert.ToDouble(alr);
        if (p.TryGetValue("CriticLearningRate", out var clr)) CriticLearningRate = Convert.ToDouble(clr);
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("Lambda", out var l)) Lambda = Convert.ToDouble(l);
        if (p.TryGetValue("ClipEpsilon", out var ce)) ClipEpsilon = Convert.ToDouble(ce);
        if (p.TryGetValue("UpdateEpochs", out var ue)) UpdateEpochs = Convert.ToInt32(ue);
        if (p.TryGetValue("MiniBatchSize", out var mb)) MiniBatchSize = Convert.ToInt32(mb);
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

    private void Shuffle(int[] array)
    {
        for (int i = array.Length - 1; i > 0; i--)
        {
            int j = _rng.Next(i + 1);
            (array[i], array[j]) = (array[j], array[i]);
        }
    }
}
