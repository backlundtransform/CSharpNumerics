using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.ML.ReinforcementLearning.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Experiment;

/// <summary>
/// Evaluates a trained agent on N greedy (no exploration) episodes
/// and computes summary statistics.
/// </summary>
public class EpisodeEvaluator
{
    private readonly IEnvironment _env;
    private readonly int _maxStepsPerEpisode;

    public EpisodeEvaluator(IEnvironment env, int maxStepsPerEpisode = 500)
    {
        _env = env;
        _maxStepsPerEpisode = maxStepsPerEpisode;
    }

    /// <summary>
    /// Run N greedy evaluation episodes and return summary statistics.
    /// </summary>
    public EvalResult Evaluate(IAgent agent, int numEpisodes, int? seed = null)
    {
        var returns = new List<double>();

        for (int i = 0; i < numEpisodes; i++)
        {
            int? episodeSeed = seed.HasValue ? seed.Value + i : null;
            double totalReturn = RunGreedyEpisode(agent, episodeSeed);
            returns.Add(totalReturn);
        }

        return new EvalResult(returns);
    }

    private double RunGreedyEpisode(IAgent agent, int? seed)
    {
        var (state, _) = _env.Reset(seed);
        double totalReturn = 0;
        bool continuous = !_env.IsDiscrete;

        for (int step = 0; step < _maxStepsPerEpisode; step++)
        {
            VectorN nextState;
            double reward;
            bool done;

            if (continuous)
            {
                var action = agent.SelectContinuousAction(state);
                (nextState, reward, done, _) = _env.Step(action);
            }
            else
            {
                int action = agent.SelectAction(state);
                (nextState, reward, done, _) = _env.Step(action);
            }

            totalReturn += reward;
            state = nextState;
            if (done) break;
        }

        return totalReturn;
    }
}

/// <summary>
/// Summary statistics from episode evaluation.
/// </summary>
public class EvalResult
{
    public List<double> Returns { get; }
    public int NumEpisodes => Returns.Count;
    public double MeanReturn { get; }
    public double StdDev { get; }
    public double MinReturn { get; }
    public double MaxReturn { get; }
    public double MedianReturn { get; }

    public EvalResult(List<double> returns)
    {
        Returns = returns;
        if (returns.Count == 0) return;

        MeanReturn = returns.Average();
        MinReturn = returns.Min();
        MaxReturn = returns.Max();

        double variance = returns.Select(r => (r - MeanReturn) * (r - MeanReturn)).Average();
        StdDev = Math.Sqrt(variance);

        var sorted = returns.OrderBy(r => r).ToList();
        MedianReturn = sorted.Count % 2 == 0
            ? (sorted[sorted.Count / 2 - 1] + sorted[sorted.Count / 2]) / 2.0
            : sorted[sorted.Count / 2];
    }

    /// <summary>
    /// Compute a confidence interval using the normal approximation.
    /// Returns (lower, upper) bounds.
    /// </summary>
    public (double lower, double upper) ConfidenceInterval(double confidence = 0.95)
    {
        if (Returns.Count <= 1) return (MeanReturn, MeanReturn);

        double z = confidence switch
        {
            >= 0.99 => 2.576,
            >= 0.95 => 1.96,
            >= 0.90 => 1.645,
            _ => 1.96
        };

        double se = StdDev / Math.Sqrt(Returns.Count);
        return (MeanReturn - z * se, MeanReturn + z * se);
    }
}
