using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.ReinforcementLearning.Algorithms.Tabular;

/// <summary>
/// Base class for tabular RL agents that maintain a Q-table stored as a Matrix (states × actions).
/// Subclasses implement the specific update rule.
/// </summary>
public abstract class TabularAgent : IAgent
{
    public double LearningRate { get; set; } = 0.1;
    public double Gamma { get; set; } = 0.99;
    public int NumStates { get; }
    public int NumActions { get; }
    public abstract string Name { get; }

    /// <summary>Q(s,a) table — rows are states, columns are actions.</summary>
    protected double[,] Q;

    private readonly Func<VectorN, int> _stateMapper;

    /// <param name="numStates">Total number of discrete states.</param>
    /// <param name="numActions">Total number of discrete actions.</param>
    /// <param name="stateMapper">Maps a VectorN state to a flat integer index.</param>
    protected TabularAgent(int numStates, int numActions, Func<VectorN, int> stateMapper)
    {
        NumStates = numStates;
        NumActions = numActions;
        _stateMapper = stateMapper;
        Q = new double[numStates, numActions];
    }

    /// <summary>Get the Q-values for a given state as a VectorN.</summary>
    public VectorN GetQValues(VectorN state)
    {
        int s = _stateMapper(state);
        var values = new double[NumActions];
        for (int a = 0; a < NumActions; a++)
            values[a] = Q[s, a];
        return new VectorN(values);
    }

    /// <summary>Map a state VectorN to its integer index.</summary>
    protected int S(VectorN state) => _stateMapper(state);

    public int SelectAction(VectorN state) =>
        ArgMax(GetQValues(state));

    public VectorN SelectContinuousAction(VectorN state) =>
        new VectorN(new double[] { SelectAction(state) });

    public abstract void Train(Transition transition);

    public virtual void TrainBatch(List<Transition> batch)
    {
        foreach (var t in batch) Train(t);
    }

    public virtual void EndEpisode(Episode episode) { }

    public abstract IAgent Clone();

    public Dictionary<string, object> GetHyperParameters() =>
        new()
        {
            ["LearningRate"] = LearningRate,
            ["Gamma"] = Gamma
        };

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (parameters.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
    }

    /// <summary>Get the full Q-table as a Matrix (for visualisation).</summary>
    public Matrix GetQTable()
    {
        return new Matrix(Q);
    }

    private static int ArgMax(VectorN v)
    {
        int best = 0;
        double bestVal = v[0];
        for (int i = 1; i < v.Length; i++)
        {
            if (v[i] > bestVal)
            {
                bestVal = v[i];
                best = i;
            }
        }
        return best;
    }
}
