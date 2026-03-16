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
/// Dueling DQN (Wang et al. 2016).
/// Decomposes Q(s,a) = V(s) + A(s,a) - mean(A(s,·)) using two separate streams
/// that share a common feature backbone.
///
/// Architecture: input → shared hidden → V-stream (1 output) + A-stream (|A| outputs)
/// Q(s,a) = V(s) + A(s,a) - (1/|A|) Σ_a' A(s,a')
/// </summary>
public class DuelingDQN : IAgent
{
    public string Name => "DuelingDQN";

    public int[] SharedLayers { get; set; } = new[] { 64 };
    public int[] ValueLayers { get; set; } = new[] { 32 };
    public int[] AdvantageLayers { get; set; } = new[] { 32 };
    public ActivationType Activation { get; set; } = ActivationType.ReLU;
    public double LearningRate { get; set; } = 0.001;
    public double Gamma { get; set; } = 0.99;
    public int TargetUpdateFrequency { get; set; } = 100;
    public int BatchSize { get; set; } = 32;
    public int MinBufferSize { get; set; } = 64;

    private NN _sharedNet;
    private NN _valueNet;
    private NN _advantageNet;
    private NN _targetSharedNet;
    private NN _targetValueNet;
    private NN _targetAdvantageNet;

    private IReplayBuffer _buffer;
    private int _stepCount;
    private int _actionSize;
    private bool _initialized;

    public void SetReplayBuffer(IReplayBuffer buffer) => _buffer = buffer;

    public void Initialize(int inputSize, int actionSize, int? seed = null)
    {
        _actionSize = actionSize;
        int s = seed ?? 123;

        int sharedOut = SharedLayers.Length > 0 ? SharedLayers[^1] : inputSize;

        _sharedNet = new NN(inputSize, SharedLayers.Length > 1 ? SharedLayers[..^1] : Array.Empty<int>(),
            sharedOut, Activation, NN.OutputMode.Linear, s);
        _valueNet = new NN(sharedOut, ValueLayers, 1, Activation, NN.OutputMode.Linear, s + 1);
        _advantageNet = new NN(sharedOut, AdvantageLayers, actionSize, Activation, NN.OutputMode.Linear, s + 2);

        _targetSharedNet = _sharedNet.Clone();
        _targetValueNet = _valueNet.Clone();
        _targetAdvantageNet = _advantageNet.Clone();

        _buffer ??= new ReplayBuffer(10000);
        _initialized = true;
    }

    private VectorN ComputeQ(NN shared, NN value, NN advantage, VectorN state)
    {
        var features = shared.Forward(state);
        var v = value.Forward(features);     // [V(s)]
        var a = advantage.Forward(features); // [A(s,0), ..., A(s,|A|-1)]

        // Mean advantage
        double meanA = 0;
        for (int i = 0; i < a.Length; i++) meanA += a[i];
        meanA /= a.Length;

        // Q(s,a) = V(s) + A(s,a) - mean(A)
        var q = new double[a.Length];
        for (int i = 0; i < a.Length; i++)
            q[i] = v[0] + a[i] - meanA;
        return new VectorN(q);
    }

    public VectorN GetQValues(VectorN state)
    {
        EnsureInitialized();
        return ComputeQ(_sharedNet, _valueNet, _advantageNet, state);
    }

    public int SelectAction(VectorN state)
    {
        var q = GetQValues(state);
        return ArgMax(q);
    }

    public VectorN SelectContinuousAction(VectorN state) =>
        new VectorN(new double[] { SelectAction(state) });

    public void Train(Transition transition)
    {
        EnsureInitialized();
        _buffer.Add(transition);
        _stepCount++;

        if (_buffer.Count < MinBufferSize) return;

        var batch = _buffer.Sample(BatchSize);

        // For each transition, compute TD error and backprop through advantage net
        // (simplified: we train the advantage stream towards the TD target)
        foreach (var t in batch)
        {
            var features = _sharedNet.Forward(t.State);
            var advValues = _advantageNet.Forward(features, out var advActs);

            double target;
            if (t.Done)
            {
                target = t.Reward;
            }
            else
            {
                var qNextTarget = ComputeQ(_targetSharedNet, _targetValueNet, _targetAdvantageNet, t.NextState);
                target = t.Reward + Gamma * MaxVal(qNextTarget);
            }

            var qCurrent = ComputeQ(_sharedNet, _valueNet, _advantageNet, t.State);
            double error = target - qCurrent[t.Action];
            double clipped = Math.Max(-1.0, Math.Min(1.0, error));

            // Update advantage stream for the taken action
            var advTarget = new VectorN(advValues.Values);
            var targetArr = (double[])advValues.Values.Clone();
            targetArr[t.Action] = advValues[t.Action] + clipped;
            advTarget = new VectorN(targetArr);

            _advantageNet.ComputeGradients(advTarget, advActs, out var dw, out var db);
            _advantageNet.ApplyGradients(dw, db, LearningRate, 1);

            // Update value stream
            var vOut = _valueNet.Forward(features, out var vActs);
            var vTarget = new VectorN(new[] { vOut[0] + clipped });
            _valueNet.ComputeGradients(vTarget, vActs, out var vdw, out var vdb);
            _valueNet.ApplyGradients(vdw, vdb, LearningRate, 1);
        }

        if (_stepCount % TargetUpdateFrequency == 0)
        {
            _targetSharedNet.CopyWeightsFrom(_sharedNet);
            _targetValueNet.CopyWeightsFrom(_valueNet);
            _targetAdvantageNet.CopyWeightsFrom(_advantageNet);
        }
    }

    public void TrainBatch(List<Transition> batch) { foreach (var t in batch) Train(t); }
    public void EndEpisode(Episode episode) { }

    public IAgent Clone() => new DuelingDQN
    {
        SharedLayers = (int[])SharedLayers.Clone(),
        ValueLayers = (int[])ValueLayers.Clone(),
        AdvantageLayers = (int[])AdvantageLayers.Clone(),
        Activation = Activation,
        LearningRate = LearningRate,
        Gamma = Gamma,
        TargetUpdateFrequency = TargetUpdateFrequency,
        BatchSize = BatchSize,
        MinBufferSize = MinBufferSize
    };

    public Dictionary<string, object> GetHyperParameters() => new()
    {
        ["SharedLayers"] = SharedLayers,
        ["ValueLayers"] = ValueLayers,
        ["AdvantageLayers"] = AdvantageLayers,
        ["LearningRate"] = LearningRate,
        ["Gamma"] = Gamma,
        ["TargetUpdateFrequency"] = TargetUpdateFrequency,
        ["BatchSize"] = BatchSize
    };

    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("SharedLayers", out var sl)) SharedLayers = (int[])sl;
        if (p.TryGetValue("ValueLayers", out var vl)) ValueLayers = (int[])vl;
        if (p.TryGetValue("AdvantageLayers", out var al)) AdvantageLayers = (int[])al;
        if (p.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("TargetUpdateFrequency", out var tu)) TargetUpdateFrequency = Convert.ToInt32(tu);
        if (p.TryGetValue("BatchSize", out var bs)) BatchSize = Convert.ToInt32(bs);
    }

    private void EnsureInitialized()
    {
        if (!_initialized)
            throw new InvalidOperationException(
                "DuelingDQN must be initialized before use. Call Initialize(inputSize, actionSize).");
    }

    private static int ArgMax(VectorN v)
    {
        int best = 0;
        for (int i = 1; i < v.Length; i++)
            if (v[i] > v[best]) best = i;
        return best;
    }

    private static double MaxVal(VectorN v)
    {
        double max = v[0];
        for (int i = 1; i < v.Length; i++)
            if (v[i] > max) max = v[i];
        return max;
    }
}
