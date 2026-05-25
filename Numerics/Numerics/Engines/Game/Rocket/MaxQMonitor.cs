using CSharpNumerics.Engines.Common;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;
using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Monitors dynamic pressure during ascent and detects Max-Q.
/// Publishes a <see cref="MaxQEvent"/> via the EventBus when peak dynamic pressure is reached.
/// Provides structural load factor (q / q_limit) for throttle decisions.
/// </summary>
public class MaxQMonitor
{
    private double _previousQ;
    private bool _maxQReached;

    /// <summary>Peak dynamic pressure recorded so far (Pa).</summary>
    public double PeakQ { get; private set; }

    /// <summary>Altitude at which Max-Q occurred (m).</summary>
    public double MaxQAltitude { get; private set; }

    /// <summary>Time at which Max-Q occurred (s).</summary>
    public double MaxQTime { get; private set; }

    /// <summary>Current dynamic pressure (Pa).</summary>
    public double CurrentQ { get; private set; }

    /// <summary>True once Max-Q has been detected.</summary>
    public bool MaxQReached => _maxQReached;

    /// <summary>Structural load limit in Pascals. Default 35 kPa (typical for orbital rockets).</summary>
    public double QLimit { get; set; } = 35000;

    /// <summary>
    /// Current structural load factor: q / q_limit. Values above 1.0 exceed structural limits.
    /// </summary>
    public double LoadFactor => QLimit > 0 ? CurrentQ / QLimit : 0;

    /// <summary>Optional event bus for publishing Max-Q event.</summary>
    public EventBus EventBus { get; set; }

    /// <summary>
    /// Updates the monitor with current flight conditions.
    /// Call once per simulation step.
    /// </summary>
    /// <param name="altitude">Current altitude in metres.</param>
    /// <param name="velocity">Current airspeed in m/s.</param>
    /// <param name="time">Current simulation time in seconds.</param>
    public void Update(double altitude, double velocity, double time)
    {
        CurrentQ = DynamicPressure.FromAltitudeAndVelocity(altitude, velocity);

        if (CurrentQ > PeakQ)
        {
            PeakQ = CurrentQ;
            MaxQAltitude = altitude;
            MaxQTime = time;
        }

        if (!_maxQReached && _previousQ > 0 && CurrentQ < _previousQ)
        {
            _maxQReached = true;
            EventBus?.Publish(new MaxQEvent(time, altitude, velocity, PeakQ));
        }

        _previousQ = CurrentQ;
    }

    /// <summary>
    /// Resets the monitor to initial state.
    /// </summary>
    public void Reset()
    {
        _previousQ = 0;
        _maxQReached = false;
        PeakQ = 0;
        MaxQAltitude = 0;
        MaxQTime = 0;
        CurrentQ = 0;
    }
}

/// <summary>
/// Event published when Max-Q (peak dynamic pressure) is reached.
/// </summary>
public class MaxQEvent
{
    /// <summary>Time of Max-Q in seconds.</summary>
    public double Time { get; }

    /// <summary>Altitude at Max-Q in metres.</summary>
    public double Altitude { get; }

    /// <summary>Velocity at Max-Q in m/s.</summary>
    public double Velocity { get; }

    /// <summary>Peak dynamic pressure in Pascals.</summary>
    public double PeakQ { get; }

    public MaxQEvent(double time, double altitude, double velocity, double peakQ)
    {
        Time = time;
        Altitude = altitude;
        Velocity = velocity;
        PeakQ = peakQ;
    }
}
