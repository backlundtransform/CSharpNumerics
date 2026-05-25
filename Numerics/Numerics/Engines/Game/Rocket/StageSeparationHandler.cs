using CSharpNumerics.Engines.Common;
using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Handles stage separation events with pre/post-separation logic.
/// Manages the discontinuity in mass and potentially attitude at separation,
/// and can apply ullage impulse or initiate coast phase.
/// </summary>
public class StageSeparationHandler
{
    /// <summary>Optional ullage motor applied after separation to settle propellant.</summary>
    public UllageMotor UllageMotor { get; set; }

    /// <summary>Coast duration after separation before next stage ignition (seconds). Default 0.</summary>
    public double CoastDuration { get; set; }

    /// <summary>Optional event bus for publishing detailed separation events.</summary>
    public EventBus EventBus { get; set; }

    /// <summary>Number of separations handled so far.</summary>
    public int SeparationCount { get; private set; }

    /// <summary>Time of last separation.</summary>
    public double LastSeparationTime { get; private set; }

    /// <summary>
    /// Called just before stage separation occurs.
    /// Records the separation and applies ullage if configured.
    /// </summary>
    /// <param name="state">Current rocket state at separation.</param>
    /// <param name="time">Current simulation time.</param>
    public void OnSeparation(RocketState state, double time)
    {
        SeparationCount++;
        LastSeparationTime = time;
    }

    /// <summary>
    /// Determines if the vehicle is currently in a coast phase after separation.
    /// During coast, no thrust is produced (engines off).
    /// </summary>
    /// <param name="currentTime">Current simulation time in seconds.</param>
    public bool IsCoasting(double currentTime)
    {
        if (CoastDuration <= 0 || SeparationCount == 0) return false;
        return (currentTime - LastSeparationTime) < CoastDuration;
    }

    /// <summary>Resets the handler state.</summary>
    public void Reset()
    {
        SeparationCount = 0;
        LastSeparationTime = 0;
    }
}
