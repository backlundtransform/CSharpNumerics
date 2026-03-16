using CSharpNumerics.Statistics.Data;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Core;

/// <summary>
/// Diagnostic output from an RL training run.
/// </summary>
public class TrainingResult
{
    /// <summary>Total return per episode.</summary>
    public List<double> EpisodeReturns { get; } = new();

    /// <summary>Steps per episode.</summary>
    public List<int> EpisodeSteps { get; } = new();

    /// <summary>Loss values recorded during training (optional, per-step or per-episode).</summary>
    public List<double> Losses { get; } = new();

    /// <summary>Exploration parameter over time (e.g. ε).</summary>
    public List<double> ExplorationValues { get; } = new();

    /// <summary>Greedy evaluation returns recorded at periodic intervals.</summary>
    public List<double> EvalReturns { get; } = new();

    /// <summary>Episode return as (episode, return) series for plotting.</summary>
    public List<Serie> ReturnCurve =>
        EpisodeReturns.Select((r, i) => new Serie { Index = i, Value = r }).ToList();

    /// <summary>Training loss as (step, loss) series for plotting.</summary>
    public List<Serie> LossCurve =>
        Losses.Select((l, i) => new Serie { Index = i, Value = l }).ToList();

    /// <summary>Exploration schedule as (episode, value) series for plotting.</summary>
    public List<Serie> ExplorationCurve =>
        ExplorationValues.Select((v, i) => new Serie { Index = i, Value = v }).ToList();

    /// <summary>Average return over all episodes.</summary>
    public double AverageReturn => EpisodeReturns.Count > 0 ? EpisodeReturns.Average() : 0;

    /// <summary>Average return over the last N episodes.</summary>
    public double AverageReturnLastN(int lastN)
    {
        if (EpisodeReturns.Count == 0) return 0;
        return EpisodeReturns.Skip(System.Math.Max(0, EpisodeReturns.Count - lastN)).Average();
    }

    /// <summary>Best single-episode return.</summary>
    public double BestReturn => EpisodeReturns.Count > 0 ? EpisodeReturns.Max() : 0;

    /// <summary>Total number of episodes completed.</summary>
    public int TotalEpisodes => EpisodeReturns.Count;
}
