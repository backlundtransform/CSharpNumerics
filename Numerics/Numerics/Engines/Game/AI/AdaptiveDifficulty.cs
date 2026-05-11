using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Game.AI;

/// <summary>
/// Monitors player performance and adjusts AI aggressiveness.
///
/// Tracks metrics like kill/death ratio, survival time, and damage dealt.
/// When the player is struggling, the AI's effective skill is reduced
/// (e.g., slower reaction times, weaker pursuit). When the player dominates,
/// the AI becomes more challenging.
///
/// The difficulty level is a continuous value in [0, 1] where 0 = easiest, 1 = hardest.
/// This value can be used to scale AI policy parameters (exploration noise,
/// decision frequency, accuracy, etc.).
/// </summary>
public class AdaptiveDifficulty
{
    private readonly Queue<double> _recentScores;
    private readonly int _windowSize;
    private double _difficulty;
    private double _targetPerformance;

    /// <summary>Current difficulty level [0, 1]. 0 = easiest, 1 = hardest.</summary>
    public double Difficulty => _difficulty;

    /// <summary>Target player performance ratio (0–1). Default 0.5 = balanced.</summary>
    public double TargetPerformance
    {
        get => _targetPerformance;
        set => _targetPerformance = Math.Clamp(value, 0.05, 0.95);
    }

    /// <summary>How fast difficulty adjusts. Higher = more reactive.</summary>
    public double AdaptationRate { get; set; } = 0.1;

    /// <summary>Minimum difficulty floor.</summary>
    public double MinDifficulty { get; set; } = 0.1;

    /// <summary>Maximum difficulty ceiling.</summary>
    public double MaxDifficulty { get; set; } = 1.0;

    /// <summary>Number of recent scores in the evaluation window.</summary>
    public int WindowSize => _windowSize;

    /// <summary>
    /// Average performance score over the recent window.
    /// </summary>
    public double AveragePerformance =>
        _recentScores.Count > 0 ? _recentScores.Average() : _targetPerformance;

    /// <summary>
    /// Creates an adaptive difficulty controller.
    /// </summary>
    /// <param name="initialDifficulty">Starting difficulty [0, 1].</param>
    /// <param name="targetPerformance">Desired player performance ratio (0–1). 0.5 = balanced.</param>
    /// <param name="windowSize">Number of recent scores to consider.</param>
    public AdaptiveDifficulty(
        double initialDifficulty = 0.5,
        double targetPerformance = 0.5,
        int windowSize = 10)
    {
        _difficulty = Math.Clamp(initialDifficulty, 0, 1);
        _targetPerformance = Math.Clamp(targetPerformance, 0.05, 0.95);
        _windowSize = Math.Max(1, windowSize);
        _recentScores = new Queue<double>();
    }

    /// <summary>
    /// Record a player performance score for a completed round/encounter.
    /// Score should be normalized to [0, 1] where 1 = player performed perfectly.
    /// </summary>
    /// <param name="performanceScore">Player performance [0, 1].</param>
    public void RecordPerformance(double performanceScore)
    {
        double clamped = Math.Clamp(performanceScore, 0, 1);
        _recentScores.Enqueue(clamped);

        while (_recentScores.Count > _windowSize)
            _recentScores.Dequeue();

        UpdateDifficulty();
    }

    /// <summary>
    /// Apply the current difficulty to an AI agent's parameters.
    /// Returns a dictionary of suggested parameter adjustments.
    /// </summary>
    /// <returns>
    /// Dictionary with keys:
    ///   "reactionDelay" — seconds of delay before AI reacts (lower difficulty = more delay)
    ///   "accuracyScale" — multiplier for AI aiming/control accuracy [0.3–1.0]
    ///   "aggressiveness" — how aggressively AI pursues (used for pursuit gain scaling)
    ///   "explorationNoise" — noise to add to AI actions (higher = more erratic/easier)
    /// </returns>
    public Dictionary<string, double> GetAIParameters()
    {
        double d = _difficulty;

        return new Dictionary<string, double>
        {
            // At d=0: 0.5s delay, at d=1: 0s delay
            ["reactionDelay"] = 0.5 * (1.0 - d),

            // At d=0: 30% accuracy, at d=1: 100% accuracy
            ["accuracyScale"] = 0.3 + 0.7 * d,

            // At d=0: 0.2 pursuit gain, at d=1: 1.0 pursuit gain
            ["aggressiveness"] = 0.2 + 0.8 * d,

            // At d=0: 0.5 noise, at d=1: 0.02 noise
            ["explorationNoise"] = 0.5 * (1.0 - d) + 0.02
        };
    }

    /// <summary>
    /// Scale a value by the current difficulty.
    /// At difficulty 0, returns minScale × value.
    /// At difficulty 1, returns value unchanged.
    /// </summary>
    public double Scale(double value, double minScale = 0.3)
    {
        double scale = minScale + (1.0 - minScale) * _difficulty;
        return value * scale;
    }

    /// <summary>
    /// Reset the difficulty tracker to initial state.
    /// </summary>
    public void Reset(double difficulty = 0.5)
    {
        _difficulty = Math.Clamp(difficulty, 0, 1);
        _recentScores.Clear();
    }

    private void UpdateDifficulty()
    {
        if (_recentScores.Count == 0) return;

        double avgPerf = _recentScores.Average();

        // If player is performing above target → increase difficulty
        // If player is performing below target → decrease difficulty
        double error = avgPerf - _targetPerformance;
        _difficulty += AdaptationRate * error;
        _difficulty = Math.Clamp(_difficulty, MinDifficulty, MaxDifficulty);
    }
}
