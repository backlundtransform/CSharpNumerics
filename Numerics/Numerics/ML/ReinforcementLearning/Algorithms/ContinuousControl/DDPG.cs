using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.ReinforcementLearning.Buffers;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.ContinuousControl;

/// <summary>
/// Deep Deterministic Policy Gradient (Lillicrap et al. 2016).
///
/// An off-policy actor-critic algorithm for continuous action spaces:
///   - Actor μ(s|θ^μ): deterministic policy, outputs continuous action directly
///   - Critic Q(s,a|θ^Q): evaluates state-action pairs
///   - Target networks μ' and Q' with Polyak averaging: θ' ← τθ + (1-τ)θ'
///   - Experience replay for sample decorrelation
///
/// Actor gradient: ∇J ≈ E[∇_a Q(s,a)|_{a=μ(s)} · ∇_θ μ(s)]
/// Critic loss: L = E[(Q(s,a) - (r + γ Q'(s', μ'(s'))))²]
/// </summary>
public class DDPG : IAgent
{
    public string Name => "DDPG";

    // ── Hyperparameters ─────────────────────────────────────
    public int[] ActorHiddenLayers { get; set; } = new[] { 64, 64 };
    public int[] CriticHiddenLayers { get; set; } = new[] { 64, 64 };
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double ActorLearningRate { get; set; } = 1e-4;
    public double CriticLearningRate { get; set; } = 1e-3;
    public double Gamma { get; set; } = 0.99;
    public double Tau { get; set; } = 0.005;
    public int BatchSize { get; set; } = 64;
    public int MinBufferSize { get; set; } = 128;
    public double ActionScale { get; set; } = 1.0;

    // ── Internal state ──────────────────────────────────────
    private NN _actorNet;
    private NN _criticNet;
    private NN _targetActorNet;
    private NN _targetCriticNet;
    private IReplayBuffer _buffer;
    private int _actionSize;
    private int _inputSize;
    private bool _initialized;

    public void SetReplayBuffer(IReplayBuffer buffer) => _buffer = buffer;

    // ── Initialization ──────────────────────────────────────

    public void Initialize(int inputSize, int actionSize, int? seed = null)
    {
        _inputSize = inputSize;
        _actionSize = actionSize;
        int s = seed ?? 123;

        // Actor: state → action (linear output, scaled by tanh * ActionScale)
        _actorNet = new NN(inputSize, ActorHiddenLayers, actionSize, Activation, NN.OutputMode.Linear, s);
        _targetActorNet = _actorNet.Clone();

        // Critic: (state, action) concatenated → Q-value (1 output)
        _criticNet = new NN(inputSize + actionSize, CriticHiddenLayers, 1, Activation, NN.OutputMode.Linear, s + 1);
        _targetCriticNet = _criticNet.Clone();

        _buffer ??= new ReplayBuffer(100000);
        _initialized = true;
    }

    // ── Action selection ────────────────────────────────────

    /// <summary>Get the deterministic action from the actor (scaled by tanh).</summary>
    public VectorN SelectContinuousAction(VectorN state)
    {
        EnsureInitialized();
        var raw = _actorNet.Forward(state);
        return ScaleAction(raw);
    }

    /// <summary>Get raw actor output (before tanh scaling) for gradient computation.</summary>
    public VectorN GetRawAction(VectorN state)
    {
        EnsureInitialized();
        return _actorNet.Forward(state);
    }

    public int SelectAction(VectorN state)
    {
        // For compatibility — round continuous action to nearest int
        var action = SelectContinuousAction(state);
        return (int)Math.Round(action[0]);
    }

    // ── Training ────────────────────────────────────────────

    public void Train(Transition transition)
    {
        EnsureInitialized();
        _buffer.Add(transition);

        if (_buffer.Count < MinBufferSize) return;

        var batch = _buffer.Sample(BatchSize);
        TrainOnBatch(batch);

        // Soft update target networks
        _targetActorNet.SoftUpdate(_actorNet, Tau);
        _targetCriticNet.SoftUpdate(_criticNet, Tau);
    }

    public void TrainBatch(List<Transition> batch)
    {
        foreach (var t in batch) Train(t);
    }

    public void EndEpisode(Episode episode) { }

    private void TrainOnBatch(List<Transition> batch)
    {
        // ── Update Critic ───────────────────────────────────
        var criticWeightGrads = _criticNet.InitWeightGrads();
        var criticBiasGrads = _criticNet.InitBiasGrads();

        foreach (var t in batch)
        {
            var sa = ConcatStateAction(t.State, t.ContinuousAction);
            var qPred = _criticNet.Forward(sa, out var criticActivations);

            double target;
            if (t.Done)
            {
                target = t.Reward;
            }
            else
            {
                var nextAction = ScaleAction(_targetActorNet.Forward(t.NextState));
                var saNext = ConcatStateAction(t.NextState, nextAction);
                var qNext = _targetCriticNet.Forward(saNext);
                target = t.Reward + Gamma * qNext[0];
            }

            var criticTarget = new VectorN(new[] { target });
            _criticNet.ComputeGradients(criticTarget, criticActivations, out var cDw, out var cDb);
            _criticNet.AccumulateGradients(criticWeightGrads, criticBiasGrads, cDw, cDb);
        }

        _criticNet.ApplyGradients(criticWeightGrads, criticBiasGrads, CriticLearningRate, batch.Count);

        // ── Update Actor ────────────────────────────────────
        // Actor gradient: move action in direction that increases Q
        // ∇_θ J = E[∇_a Q(s,a)|_{a=μ(s)} · ∇_θ μ(s)]
        var actorWeightGrads = _actorNet.InitWeightGrads();
        var actorBiasGrads = _actorNet.InitBiasGrads();

        foreach (var t in batch)
        {
            var rawAction = _actorNet.Forward(t.State, out var actorActivations);
            var scaledAction = ScaleAction(rawAction);
            var sa = ConcatStateAction(t.State, scaledAction);
            var qValue = _criticNet.Forward(sa);

            // Compute ∇_a Q(s,a) via finite differences
            var actionGrad = new double[_actionSize];
            double eps = 1e-4;
            for (int i = 0; i < _actionSize; i++)
            {
                var actionPlus = new double[scaledAction.Length];
                var actionMinus = new double[scaledAction.Length];
                for (int j = 0; j < scaledAction.Length; j++)
                {
                    actionPlus[j] = scaledAction[j];
                    actionMinus[j] = scaledAction[j];
                }
                actionPlus[i] += eps;
                actionMinus[i] -= eps;

                var qPlus = _criticNet.Forward(ConcatStateAction(t.State, new VectorN(actionPlus)));
                var qMinus = _criticNet.Forward(ConcatStateAction(t.State, new VectorN(actionMinus)));
                actionGrad[i] = (qPlus[0] - qMinus[0]) / (2 * eps);
            }

            // Chain rule through tanh: ∂scaled/∂raw = ActionScale * (1 - tanh²(raw))
            // Policy gradient target: push raw action in direction of ∇_a Q · ∂scaled/∂raw
            var actorTarget = new double[_actionSize];
            for (int i = 0; i < _actionSize; i++)
            {
                double tanhVal = Math.Tanh(rawAction[i]);
                double tanhDeriv = ActionScale * (1.0 - tanhVal * tanhVal);
                double grad = actionGrad[i] * tanhDeriv;
                // Target = rawAction + grad (ascent direction — we want to maximize Q)
                actorTarget[i] = rawAction[i] - grad; // ComputeGradients minimizes (target - output)
            }

            _actorNet.ComputeGradients(new VectorN(actorTarget), actorActivations, out var aDw, out var aDb);
            _actorNet.AccumulateGradients(actorWeightGrads, actorBiasGrads, aDw, aDb);
        }

        _actorNet.ApplyGradients(actorWeightGrads, actorBiasGrads, ActorLearningRate, batch.Count);
    }

    // ── Clone / Hyperparameters ─────────────────────────────

    public IAgent Clone() => new DDPG
    {
        ActorHiddenLayers = (int[])ActorHiddenLayers.Clone(),
        CriticHiddenLayers = (int[])CriticHiddenLayers.Clone(),
        Activation = Activation,
        ActorLearningRate = ActorLearningRate,
        CriticLearningRate = CriticLearningRate,
        Gamma = Gamma,
        Tau = Tau,
        BatchSize = BatchSize,
        MinBufferSize = MinBufferSize,
        ActionScale = ActionScale
    };

    public Dictionary<string, object> GetHyperParameters() => new()
    {
        ["ActorHiddenLayers"] = ActorHiddenLayers,
        ["CriticHiddenLayers"] = CriticHiddenLayers,
        ["ActorLearningRate"] = ActorLearningRate,
        ["CriticLearningRate"] = CriticLearningRate,
        ["Gamma"] = Gamma,
        ["Tau"] = Tau,
        ["BatchSize"] = BatchSize,
        ["ActionScale"] = ActionScale
    };

    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("ActorHiddenLayers", out var ah)) ActorHiddenLayers = (int[])ah;
        if (p.TryGetValue("CriticHiddenLayers", out var ch)) CriticHiddenLayers = (int[])ch;
        if (p.TryGetValue("ActorLearningRate", out var alr)) ActorLearningRate = Convert.ToDouble(alr);
        if (p.TryGetValue("CriticLearningRate", out var clr)) CriticLearningRate = Convert.ToDouble(clr);
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("Tau", out var tau)) Tau = Convert.ToDouble(tau);
        if (p.TryGetValue("BatchSize", out var bs)) BatchSize = Convert.ToInt32(bs);
        if (p.TryGetValue("ActionScale", out var asc)) ActionScale = Convert.ToDouble(asc);
    }

    // ── Helpers ─────────────────────────────────────────────

    private VectorN ScaleAction(VectorN raw)
    {
        var scaled = new double[raw.Length];
        for (int i = 0; i < raw.Length; i++)
            scaled[i] = ActionScale * Math.Tanh(raw[i]);
        return new VectorN(scaled);
    }

    private VectorN ConcatStateAction(VectorN state, VectorN action)
    {
        var concat = new double[state.Length + action.Length];
        Array.Copy(state.Values, 0, concat, 0, state.Length);
        Array.Copy(action.Values, 0, concat, state.Length, action.Length);
        return new VectorN(concat);
    }

    private void EnsureInitialized()
    {
        if (!_initialized)
            throw new InvalidOperationException(
                "DDPG must be initialized before use. Call Initialize(inputSize, actionSize) " +
                "or use RLExperiment which initializes automatically.");
    }
}
