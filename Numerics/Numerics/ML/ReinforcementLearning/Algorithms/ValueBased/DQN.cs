using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.ReinforcementLearning.Buffers;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.ValueBased;

/// <summary>
/// Deep Q-Network (Mnih et al. 2015).
/// Uses a neural network to approximate Q(s,a), with:
///   - Experience replay buffer for decorrelation
///   - Target network for stable TD targets (hard copy every N steps)
///   - Huber loss for robust gradient computation
/// </summary>
public class DQN : IAgent
{
    public string Name => "DQN";

    // ── Hyperparameters ─────────────────────────────────────
    public int[] HiddenLayers { get; set; } = new[] { 64, 64 };
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double LearningRate { get; set; } = 0.001;
    public double Gamma { get; set; } = 0.99;
    public int TargetUpdateFrequency { get; set; } = 100;
    public int BatchSize { get; set; } = 32;
    public int MinBufferSize { get; set; } = 64;

    // ── Internal state ──────────────────────────────────────
    protected NN _qNetwork;
    protected NN _targetNetwork;
    private IReplayBuffer _buffer;
    private int _stepCount;
    private int _inputSize;
    private int _actionSize;
    private bool _initialized;

    /// <summary>Inject an external replay buffer (set by RLExperimentBuilder).</summary>
    public void SetReplayBuffer(IReplayBuffer buffer) => _buffer = buffer;

    // ── Action selection ────────────────────────────────────

    public int SelectAction(VectorN state)
    {
        EnsureInitialized(state.Length);
        var qValues = _qNetwork.Forward(state);
        return ArgMax(qValues);
    }

    /// <summary>Get Q-values for a state (used by policy for ε-greedy).</summary>
    public VectorN GetQValues(VectorN state)
    {
        EnsureInitialized(state.Length);
        return _qNetwork.Forward(state);
    }

    public VectorN SelectContinuousAction(VectorN state) =>
        new VectorN(new double[] { SelectAction(state) });

    // ── Training ────────────────────────────────────────────

    public void Train(Transition transition)
    {
        EnsureInitialized(transition.State.Length);

        _buffer.Add(transition);
        _stepCount++;

        if (_buffer.Count < MinBufferSize) return;

        var batch = _buffer.Sample(BatchSize);
        TrainOnBatch(batch);

        if (_stepCount % TargetUpdateFrequency == 0)
            _targetNetwork.CopyWeightsFrom(_qNetwork);
    }

    public void TrainBatch(List<Transition> batch)
    {
        foreach (var t in batch) Train(t);
    }

    public void EndEpisode(Episode episode) { }

    /// <summary>Core training step on a mini-batch. Virtual so DoubleDQN can override target computation.</summary>
    protected virtual void TrainOnBatch(List<Transition> batch)
    {
        var weightGrads = _qNetwork.InitWeightGrads();
        var biasGrads = _qNetwork.InitBiasGrads();

        foreach (var t in batch)
        {
            var qValues = _qNetwork.Forward(t.State, out var activations);

            // TD target: r + γ max_a' Q_target(s', a')
            double target;
            if (t.Done)
            {
                target = t.Reward;
            }
            else
            {
                var qNext = _targetNetwork.Forward(t.NextState);
                target = t.Reward + Gamma * Max(qNext);
            }

            // Compute gradient: only for the chosen action
            var targetVec = new VectorN(qValues.Values);
            double error = target - qValues[t.Action];
            // Huber loss gradient: clip to [-1, 1]
            double clippedError = Math.Max(-1.0, Math.Min(1.0, error));
            var targetValues = (double[])qValues.Values.Clone();
            targetValues[t.Action] = qValues[t.Action] + clippedError;
            targetVec = new VectorN(targetValues);

            _qNetwork.ComputeGradients(targetVec, activations, out var dw, out var db);
            _qNetwork.AccumulateGradients(weightGrads, biasGrads, dw, db);
        }

        _qNetwork.ApplyGradients(weightGrads, biasGrads, LearningRate, batch.Count);
    }

    // ── Initialization ──────────────────────────────────────

    /// <summary>
    /// Initialize networks and buffer. Called lazily on first use when input size is known.
    /// Can also be called explicitly with SetActionSize to pre-configure.
    /// </summary>
    public void Initialize(int inputSize, int actionSize, int? seed = null)
    {
        _inputSize = inputSize;
        _actionSize = actionSize;
        int s = seed ?? 123;
        _qNetwork = new NN(inputSize, HiddenLayers, actionSize, Activation, NN.OutputMode.Linear, s);
        _targetNetwork = new NN(inputSize, HiddenLayers, actionSize, Activation, NN.OutputMode.Linear, s);
        _buffer ??= new ReplayBuffer(10000);
        _initialized = true;
    }

    private void EnsureInitialized(int inputSize)
    {
        if (!_initialized)
            throw new InvalidOperationException(
                "DQN must be initialized before use. Call Initialize(inputSize, actionSize) " +
                "or use RLExperiment which initializes automatically.");
    }

    // ── Clone ───────────────────────────────────────────────

    public virtual IAgent Clone()
    {
        return new DQN
        {
            HiddenLayers = (int[])HiddenLayers.Clone(),
            Activation = Activation,
            LearningRate = LearningRate,
            Gamma = Gamma,
            TargetUpdateFrequency = TargetUpdateFrequency,
            BatchSize = BatchSize,
            MinBufferSize = MinBufferSize
        };
    }

    public Dictionary<string, object> GetHyperParameters() => new()
    {
        ["HiddenLayers"] = HiddenLayers,
        ["Activation"] = Activation,
        ["LearningRate"] = LearningRate,
        ["Gamma"] = Gamma,
        ["TargetUpdateFrequency"] = TargetUpdateFrequency,
        ["BatchSize"] = BatchSize
    };

    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("HiddenLayers", out var h)) HiddenLayers = (int[])h;
        if (p.TryGetValue("Activation", out var a)) Activation = (ActivationType)a;
        if (p.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("TargetUpdateFrequency", out var tu)) TargetUpdateFrequency = Convert.ToInt32(tu);
        if (p.TryGetValue("BatchSize", out var bs)) BatchSize = Convert.ToInt32(bs);
    }

    // ── Helpers ─────────────────────────────────────────────

    protected static int ArgMax(VectorN v)
    {
        int best = 0;
        double bestVal = v[0];
        for (int i = 1; i < v.Length; i++)
            if (v[i] > bestVal) { bestVal = v[i]; best = i; }
        return best;
    }

    protected static double Max(VectorN v)
    {
        double max = v[0];
        for (int i = 1; i < v.Length; i++)
            if (v[i] > max) max = v[i];
        return max;
    }
}
