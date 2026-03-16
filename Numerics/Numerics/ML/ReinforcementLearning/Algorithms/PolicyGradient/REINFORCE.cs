using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;

/// <summary>
/// REINFORCE (Williams, 1992) — vanilla Monte Carlo policy gradient.
///
/// The policy network outputs a softmax probability distribution π(a|s).
/// After each episode, the gradient is:
///   ∇J ≈ Σ_t (G_t - b) ∇ log π(a_t|s_t)
///
/// where G_t is the discounted return from step t and b is a baseline
/// (average return so far) for variance reduction.
/// </summary>
public class REINFORCE : IAgent
{
    public string Name => "REINFORCE";

    // ── Hyperparameters ─────────────────────────────────────
    public int[] HiddenLayers { get; set; } = new[] { 32 };
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double LearningRate { get; set; } = 0.01;
    public double Gamma { get; set; } = 0.99;
    public bool UseBaseline { get; set; } = true;

    // ── Internal state ──────────────────────────────────────
    private NN _policyNet;
    private int _actionSize;
    private bool _initialized;
    private Random _rng;
    private double _baselineReturn;
    private int _episodeCount;

    // Stored log-probs and states for the current episode
    private List<VectorN> _savedStates;
    private List<int> _savedActions;

    public void Initialize(int inputSize, int actionSize, int? seed = null)
    {
        _actionSize = actionSize;
        int s = seed ?? 123;
        _policyNet = new NN(inputSize, HiddenLayers, actionSize, Activation, NN.OutputMode.Softmax, s);
        _rng = new Random(s);
        _savedStates = new List<VectorN>();
        _savedActions = new List<int>();
        _baselineReturn = 0;
        _episodeCount = 0;
        _initialized = true;
    }

    // ── Action selection ────────────────────────────────────

    public int SelectAction(VectorN state)
    {
        EnsureInitialized();
        var probs = _policyNet.Forward(state);
        int action = SampleCategorical(probs);

        // Store for gradient computation in EndEpisode
        _savedStates.Add(state);
        _savedActions.Add(action);

        return action;
    }

    /// <summary>Get action probabilities for a state.</summary>
    public VectorN GetActionProbabilities(VectorN state)
    {
        EnsureInitialized();
        return _policyNet.Forward(state);
    }

    public VectorN SelectContinuousAction(VectorN state) =>
        new VectorN(new double[] { SelectAction(state) });

    // ── Training ────────────────────────────────────────────

    /// <summary>No-op for REINFORCE — all learning happens in EndEpisode.</summary>
    public void Train(Transition transition) { }

    public void TrainBatch(List<Transition> batch) { }

    /// <summary>
    /// Compute Monte Carlo returns and apply REINFORCE gradient update.
    /// </summary>
    public void EndEpisode(Episode episode)
    {
        EnsureInitialized();
        if (episode.Transitions.Count == 0) return;

        var returns = episode.DiscountedReturns(Gamma);

        // Update baseline (running average of episode return)
        _episodeCount++;
        double totalReturn = episode.TotalReturn;
        _baselineReturn += (totalReturn - _baselineReturn) / _episodeCount;

        // Accumulate policy gradients
        var weightGrads = _policyNet.InitWeightGrads();
        var biasGrads = _policyNet.InitBiasGrads();

        for (int t = 0; t < _savedStates.Count && t < returns.Length; t++)
        {
            double G = returns[t];
            double advantage = UseBaseline ? G - _baselineReturn : G;

            var probs = _policyNet.Forward(_savedStates[t], out var activations);

            // Policy gradient: ∇ log π(a|s) * advantage
            // For softmax output: target = π(a|s) + lr * advantage * (1{a} - π(a|s))
            // Equivalently, construct a gradient that pushes probability toward/away from action a
            var gradOutput = new double[_actionSize];
            for (int a = 0; a < _actionSize; a++)
            {
                double indicator = (a == _savedActions[t]) ? 1.0 : 0.0;
                // ∇ log π(a_t|s_t) for softmax = (1{a=a_t} - π(a|s_t))
                // We want to move in direction: advantage * ∇ log π
                // For the MSE-based backward pass: target - output = desired direction
                // So target = output - advantage * (indicator - prob)
                gradOutput[a] = probs[a] - advantage * (indicator - probs[a]);
            }
            var target = new VectorN(gradOutput);

            _policyNet.ComputeGradients(target, activations, out var dW, out var dB);
            _policyNet.AccumulateGradients(weightGrads, biasGrads, dW, dB);
        }

        _policyNet.ApplyGradients(weightGrads, biasGrads, LearningRate, Math.Max(1, _savedStates.Count));

        // Clear episode storage
        _savedStates.Clear();
        _savedActions.Clear();
    }

    // ── Initialization ──────────────────────────────────────

    private void EnsureInitialized()
    {
        if (!_initialized)
            throw new InvalidOperationException(
                "REINFORCE must be initialized before use. Call Initialize(inputSize, actionSize) " +
                "or use RLExperiment which initializes automatically.");
    }

    // ── Clone / Hyperparameters ─────────────────────────────

    public IAgent Clone() => new REINFORCE
    {
        HiddenLayers = (int[])HiddenLayers.Clone(),
        Activation = Activation,
        LearningRate = LearningRate,
        Gamma = Gamma,
        UseBaseline = UseBaseline
    };

    public Dictionary<string, object> GetHyperParameters() => new()
    {
        ["HiddenLayers"] = HiddenLayers,
        ["Activation"] = Activation,
        ["LearningRate"] = LearningRate,
        ["Gamma"] = Gamma,
        ["UseBaseline"] = UseBaseline
    };

    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("HiddenLayers", out var h)) HiddenLayers = (int[])h;
        if (p.TryGetValue("Activation", out var a)) Activation = (ActivationType)a;
        if (p.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("UseBaseline", out var b)) UseBaseline = Convert.ToBoolean(b);
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
