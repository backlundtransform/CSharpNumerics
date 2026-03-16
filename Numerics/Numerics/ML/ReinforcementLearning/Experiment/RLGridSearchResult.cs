using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Experiment;

/// <summary>
/// Result of a grid search over RL agent configurations.
/// Rankings ordered by descending evaluation return (best first).
/// </summary>
public class RLGridSearchResult
{
    public List<RLGridSearchEntry> Rankings { get; set; } = new();

    /// <summary>Total wall-clock time for the entire grid search.</summary>
    public TimeSpan TotalDuration { get; set; }

    // ── Convenience accessors ───────────────────────────────────

    /// <summary>Best ranked entry (highest eval return).</summary>
    public RLGridSearchEntry Best => Rankings.FirstOrDefault();

    /// <summary>Best average evaluation return.</summary>
    public double BestScore => Best?.AverageEvalReturn ?? double.NaN;

    /// <summary>Name of the best agent type.</summary>
    public string BestAgentName => Best?.AgentName;
}

/// <summary>
/// A single entry in the grid search results — one trained agent configuration.
/// </summary>
public class RLGridSearchEntry
{
    /// <summary>Agent type name (e.g. "DQN", "PPO").</summary>
    public string AgentName { get; set; }

    /// <summary>Hyperparameters for this configuration.</summary>
    public Dictionary<string, object> Parameters { get; set; } = new();

    /// <summary>Human-readable description.</summary>
    public string Description { get; set; }

    /// <summary>Average return over evaluation episodes.</summary>
    public double AverageEvalReturn { get; set; }

    /// <summary>Detailed evaluation result.</summary>
    public EvalResult EvalResult { get; set; }

    /// <summary>Full training result with curves.</summary>
    public RLExperimentResult TrainingResult { get; set; }

    /// <summary>Wall-clock time for training this configuration.</summary>
    public TimeSpan Duration { get; set; }
}

/// <summary>
/// Result of Monte Carlo evaluation — multiple independent training runs.
/// </summary>
public class RLMonteCarloResult
{
    /// <summary>Results from each independent training run.</summary>
    public List<RLExperimentResult> Runs { get; set; } = new();

    /// <summary>Evaluation results from each run (greedy post-training eval).</summary>
    public List<EvalResult> EvalResults { get; set; } = new();

    /// <summary>Mean of the mean evaluation returns across all runs.</summary>
    public double MeanReturn => EvalResults.Count > 0
        ? EvalResults.Average(e => e.MeanReturn) : double.NaN;

    /// <summary>Standard deviation of mean returns across runs.</summary>
    public double StdDev
    {
        get
        {
            if (EvalResults.Count <= 1) return 0;
            double mean = MeanReturn;
            double variance = EvalResults.Select(e => (e.MeanReturn - mean) * (e.MeanReturn - mean)).Average();
            return Math.Sqrt(variance);
        }
    }

    /// <summary>Number of independent runs.</summary>
    public int NumRuns => Runs.Count;

    /// <summary>
    /// Confidence interval on the mean return across runs.
    /// </summary>
    public (double lower, double upper) ConfidenceInterval(double confidence = 0.95)
    {
        if (EvalResults.Count <= 1) return (MeanReturn, MeanReturn);

        double z = confidence switch
        {
            >= 0.99 => 2.576,
            >= 0.95 => 1.96,
            >= 0.90 => 1.645,
            _ => 1.96
        };

        double se = StdDev / Math.Sqrt(EvalResults.Count);
        return (MeanReturn - z * se, MeanReturn + z * se);
    }
}
