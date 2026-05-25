using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Monitors flight conditions and detects when staging events should occur.
/// Uses step-halving (bisection) to locate the exact moment a trigger condition
/// is crossed, providing sub-timestep precision for staging events.
/// </summary>
public class EventDetector
{
    /// <summary>Tolerance for bisection event location (seconds). Default 0.01s.</summary>
    public double Tolerance { get; set; } = 0.01;

    /// <summary>Maximum bisection iterations. Default 10 (~1/1024 of original step).</summary>
    public int MaxIterations { get; set; } = 10;

    /// <summary>
    /// Locates the precise time a condition transitions from false to true
    /// within a time interval [t0, t1] using bisection.
    /// </summary>
    /// <param name="condition">Function that evaluates the trigger condition at a given time fraction (0–1).</param>
    /// <param name="dt">Full time step duration in seconds.</param>
    /// <returns>Fraction of dt at which the event occurs (0–1), or -1 if no event detected.</returns>
    public double LocateEvent(Func<double, bool> condition, double dt)
    {
        // Check if event occurs at all during this step
        bool startState = condition(0);
        bool endState = condition(1);

        if (startState || !endState)
            return startState ? 0 : -1;

        // Bisect to find the transition point
        double lo = 0, hi = 1;
        for (int i = 0; i < MaxIterations; i++)
        {
            double mid = (lo + hi) * 0.5;
            if (condition(mid))
                hi = mid;
            else
                lo = mid;

            if ((hi - lo) * dt < Tolerance)
                break;
        }

        return (lo + hi) * 0.5;
    }

    /// <summary>
    /// Checks a fuel depletion condition with interpolated propellant mass.
    /// </summary>
    /// <param name="propellantMass">Current propellant mass (kg).</param>
    /// <param name="massFlowRate">Mass flow rate (kg/s).</param>
    /// <param name="dt">Time step (seconds).</param>
    /// <returns>Fraction of dt at which fuel depletes, or -1 if not during this step.</returns>
    public double LocateFuelDepletion(double propellantMass, double massFlowRate, double dt)
    {
        if (massFlowRate <= 0 || propellantMass <= 0) return -1;

        double burnTime = propellantMass / massFlowRate;
        if (burnTime >= dt) return -1;

        return burnTime / dt;
    }

    /// <summary>
    /// Checks an altitude threshold crossing.
    /// </summary>
    /// <param name="altitudeStart">Altitude at start of step (m).</param>
    /// <param name="altitudeEnd">Altitude at end of step (m).</param>
    /// <param name="threshold">Altitude threshold (m).</param>
    /// <param name="dt">Time step (seconds).</param>
    /// <returns>Fraction of dt at which threshold is crossed, or -1 if not crossed.</returns>
    public double LocateAltitudeCrossing(double altitudeStart, double altitudeEnd, double threshold, double dt)
    {
        if (altitudeStart >= threshold || altitudeEnd < threshold) return -1;

        // Linear interpolation for crossing point
        double fraction = (threshold - altitudeStart) / (altitudeEnd - altitudeStart);
        return Math.Clamp(fraction, 0, 1);
    }
}
