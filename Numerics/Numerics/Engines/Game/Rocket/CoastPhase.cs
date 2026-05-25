using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Represents a ballistic coast phase between stage separation and next stage ignition.
/// During coast, no thrust is produced — only gravity and aerodynamic drag act on the vehicle.
/// 
/// Coast phases are common in multi-stage rockets to allow the vehicle to reach
/// optimal altitude/velocity for the next burn, or to allow propellant settling.
/// </summary>
public class CoastPhase
{
    /// <summary>Duration of the coast phase in seconds.</summary>
    public double Duration { get; }

    /// <summary>Time at which coast began (set when Start is called).</summary>
    public double StartTime { get; private set; } = -1;

    /// <summary>True once the coast phase has been started.</summary>
    public bool IsStarted => StartTime >= 0;

    /// <summary>True if the coast phase is currently active (started and not yet complete).</summary>
    public bool IsActive(double currentTime)
    {
        if (!IsStarted) return false;
        return (currentTime - StartTime) < Duration;
    }

    /// <summary>True if the coast phase has completed.</summary>
    public bool IsComplete(double currentTime)
    {
        if (!IsStarted) return false;
        return (currentTime - StartTime) >= Duration;
    }

    /// <summary>
    /// Fraction of coast completed (0–1).
    /// </summary>
    public double Progress(double currentTime)
    {
        if (!IsStarted) return 0;
        return Math.Clamp((currentTime - StartTime) / Duration, 0, 1);
    }

    /// <summary>
    /// Creates a coast phase with the given duration.
    /// </summary>
    /// <param name="duration">Coast duration in seconds.</param>
    public CoastPhase(double duration)
    {
        if (duration <= 0) throw new ArgumentOutOfRangeException(nameof(duration));
        Duration = duration;
    }

    /// <summary>
    /// Starts the coast phase at the given simulation time.
    /// </summary>
    /// <param name="time">Current simulation time in seconds.</param>
    public void Start(double time)
    {
        StartTime = time;
    }

    /// <summary>
    /// Resets the coast phase to initial state.
    /// </summary>
    public void Reset()
    {
        StartTime = -1;
    }
}
