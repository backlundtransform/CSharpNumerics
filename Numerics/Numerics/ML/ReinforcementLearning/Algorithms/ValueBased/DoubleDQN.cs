using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.ValueBased;

/// <summary>
/// Double DQN (van Hasselt et al. 2016).
/// Decouples action selection from action evaluation to reduce overestimation bias:
///   a* = argmax_a Q_online(s', a)
///   target = r + γ Q_target(s', a*)
/// </summary>
public class DoubleDQN : DQN
{
    public new string Name => "DoubleDQN";

    protected override void TrainOnBatch(List<Transition> batch)
    {
        var weightGrads = _qNetwork.InitWeightGrads();
        var biasGrads = _qNetwork.InitBiasGrads();

        foreach (var t in batch)
        {
            var qValues = _qNetwork.Forward(t.State, out var activations);

            double target;
            if (t.Done)
            {
                target = t.Reward;
            }
            else
            {
                // Double DQN: select action with online net, evaluate with target net
                var qNextOnline = _qNetwork.Forward(t.NextState);
                int bestAction = ArgMax(qNextOnline);
                var qNextTarget = _targetNetwork.Forward(t.NextState);
                target = t.Reward + Gamma * qNextTarget[bestAction];
            }

            double error = target - qValues[t.Action];
            double clippedError = Math.Max(-1.0, Math.Min(1.0, error));
            var targetValues = (double[])qValues.Values.Clone();
            targetValues[t.Action] = qValues[t.Action] + clippedError;
            var targetVec = new VectorN(targetValues);

            _qNetwork.ComputeGradients(targetVec, activations, out var dw, out var db);
            _qNetwork.AccumulateGradients(weightGrads, biasGrads, dw, db);
        }

        _qNetwork.ApplyGradients(weightGrads, biasGrads, LearningRate, batch.Count);
    }

    public override IAgent Clone()
    {
        return new DoubleDQN
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
}
