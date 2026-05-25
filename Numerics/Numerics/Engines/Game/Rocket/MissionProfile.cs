using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Defines a sequence of flight phases from launch to orbit.
/// Each phase has entry/exit conditions and guidance mode.
/// </summary>
public class MissionProfile
{
    private readonly List<MissionPhase> _phases = new List<MissionPhase>();
    private int _activePhaseIndex;

    /// <summary>All defined mission phases in order.</summary>
    public IReadOnlyList<MissionPhase> Phases => _phases;

    /// <summary>Index of the currently active phase.</summary>
    public int ActivePhaseIndex => _activePhaseIndex;

    /// <summary>Currently active phase, or null if mission is complete.</summary>
    public MissionPhase ActivePhase =>
        _activePhaseIndex < _phases.Count ? _phases[_activePhaseIndex] : null;

    /// <summary>Whether all phases have been completed.</summary>
    public bool IsComplete => _activePhaseIndex >= _phases.Count;

    /// <summary>Adds a phase to the mission profile.</summary>
    public void AddPhase(MissionPhase phase)
    {
        _phases.Add(phase);
    }

    /// <summary>
    /// Checks if the current phase's exit condition is met and advances to the next phase.
    /// </summary>
    /// <param name="time">Current simulation time (s).</param>
    /// <param name="altitude">Current altitude (m).</param>
    /// <param name="velocity">Current speed (m/s).</param>
    /// <returns>True if a phase transition occurred.</returns>
    public bool Update(double time, double altitude, double velocity)
    {
        if (IsComplete) return false;

        var phase = ActivePhase;
        if (phase.ExitCondition(time, altitude, velocity))
        {
            phase.IsCompleted = true;
            _activePhaseIndex++;
            if (!IsComplete)
                ActivePhase.StartTime = time;
            return true;
        }
        return false;
    }

    /// <summary>Resets the mission profile to the first phase.</summary>
    public void Reset()
    {
        _activePhaseIndex = 0;
        foreach (var phase in _phases)
        {
            phase.IsCompleted = false;
            phase.StartTime = 0;
        }
    }

    /// <summary>
    /// Creates a standard ascent-to-orbit mission profile.
    /// </summary>
    public static MissionProfile StandardAscentToOrbit(double targetAltitude = 400000)
    {
        var profile = new MissionProfile();

        profile.AddPhase(new MissionPhase("Vertical Ascent",
            MissionPhaseType.VerticalAscent,
            (t, alt, v) => alt > 1000)); // Clear tower

        profile.AddPhase(new MissionPhase("Gravity Turn",
            MissionPhaseType.GravityTurn,
            (t, alt, v) => alt > 80000)); // Above atmosphere

        profile.AddPhase(new MissionPhase("Coast to Apoapsis",
            MissionPhaseType.Coast,
            (t, alt, v) => alt > targetAltitude * 0.95));

        profile.AddPhase(new MissionPhase("Circularization",
            MissionPhaseType.Circularization,
            (t, alt, v) => v > Math.Sqrt(3.986004418e14 / (6378137 + targetAltitude)) * 0.99));

        profile.AddPhase(new MissionPhase("Orbit",
            MissionPhaseType.Orbit,
            (t, alt, v) => false)); // Never exits

        return profile;
    }
}

/// <summary>
/// A single phase of a mission profile.
/// </summary>
public class MissionPhase
{
    /// <summary>Human-readable name of the phase.</summary>
    public string Name { get; }

    /// <summary>Phase type enum for guidance selection.</summary>
    public MissionPhaseType Type { get; }

    /// <summary>Condition that triggers transition to the next phase.</summary>
    public Func<double, double, double, bool> ExitCondition { get; }

    /// <summary>Time at which this phase started.</summary>
    public double StartTime { get; set; }

    /// <summary>Whether this phase has been completed.</summary>
    public bool IsCompleted { get; set; }

    public MissionPhase(string name, MissionPhaseType type, Func<double, double, double, bool> exitCondition)
    {
        Name = name;
        Type = type;
        ExitCondition = exitCondition;
    }
}

/// <summary>
/// Types of mission phases for guidance mode selection.
/// </summary>
public enum MissionPhaseType
{
    VerticalAscent,
    GravityTurn,
    Coast,
    Circularization,
    Orbit
}
