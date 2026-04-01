using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Runtime configuration for a wildfire spread simulation.
/// </summary>
public class WildfireParameters
{
    /// <summary>
    /// Ignition points as (ix, iy) cell indices.
    /// </summary>
    public IReadOnlyList<(int ix, int iy)> IgnitionPoints { get; }

    /// <summary>Mid-flame wind speed in m/s.</summary>
    public double MidflameWindSpeed { get; }

    /// <summary>
    /// Horizontal wind direction vector (only x, y components used).
    /// The vector points in the direction the wind is blowing towards.
    /// </summary>
    public Vector WindDirection { get; }

    /// <summary>
    /// How long (seconds) a cell stays in <c>Burning</c> state before
    /// transitioning to <c>Burned</c>. Default 600 s (10 min).
    /// </summary>
    public double BurnDuration { get; }

    /// <summary>Enable spot-fire ember transport (future, default false).</summary>
    public bool SpotFireEnabled { get; }

    /// <summary>
    /// Creates wildfire simulation parameters.
    /// </summary>
    public WildfireParameters(
        IReadOnlyList<(int ix, int iy)> ignitionPoints,
        double midflameWindSpeed = 0,
        Vector windDirection = default,
        double burnDuration = 600,
        bool spotFireEnabled = false)
    {
        if (ignitionPoints == null || ignitionPoints.Count == 0)
            throw new ArgumentException("At least one ignition point is required.", nameof(ignitionPoints));

        IgnitionPoints = ignitionPoints;
        MidflameWindSpeed = midflameWindSpeed;
        WindDirection = windDirection;
        BurnDuration = burnDuration;
        SpotFireEnabled = spotFireEnabled;
    }
}
