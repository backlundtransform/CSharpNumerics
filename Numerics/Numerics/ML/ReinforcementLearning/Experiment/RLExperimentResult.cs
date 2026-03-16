using CSharpNumerics.ML.ReinforcementLearning.Core;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.ReinforcementLearning.Experiment;

/// <summary>
/// Results container for an RL experiment.
/// </summary>
public class RLExperimentResult
{
    /// <summary>The training diagnostics (return curves, loss, etc.).</summary>
    public TrainingResult Training { get; internal set; }

    /// <summary>Name of the agent that was trained.</summary>
    public string AgentName { get; internal set; }

    /// <summary>Hyperparameters used.</summary>
    public Dictionary<string, object> Parameters { get; internal set; }

    /// <summary>Total wall-clock training time.</summary>
    public TimeSpan Duration { get; internal set; }

    // ── Convenience accessors (delegate to Training) ─────────

    /// <summary>Average return over all training episodes.</summary>
    public double AverageReturn => Training.AverageReturn;

    /// <summary>Average return over the last N episodes.</summary>
    public double AverageReturnLastN(int lastN) => Training.AverageReturnLastN(lastN);

    /// <summary>Best single-episode return.</summary>
    public double BestReturn => Training.BestReturn;

    /// <summary>Total episodes trained.</summary>
    public int TotalEpisodes => Training.TotalEpisodes;

    /// <summary>Episode return curve as List&lt;Serie&gt;.</summary>
    public List<Serie> ReturnCurve => Training.ReturnCurve;

    /// <summary>Training loss curve as List&lt;Serie&gt;.</summary>
    public List<Serie> LossCurve => Training.LossCurve;

    /// <summary>Exploration parameter curve as List&lt;Serie&gt;.</summary>
    public List<Serie> ExplorationCurve => Training.ExplorationCurve;
}
