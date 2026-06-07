using System;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Defines stochastic parameter variation ranges for Monte Carlo wildfire simulations.
/// Each nullable property pair describes the range from which the MC runner will sample.
/// When null, the corresponding parameter is held at its deterministic baseline value.
/// </summary>
public class WildfireVariation
{
    public double? WindSpeedMin { get; private set; }
    public double? WindSpeedMax { get; private set; }
    public double? WindDirectionJitterDeg { get; private set; }
    public double? MoistureMin { get; private set; }
    public double? MoistureMax { get; private set; }
    public double? IgnitionOffsetRadius { get; private set; }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent setters
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Uniform random wind speed range (m/s).</summary>
    public WildfireVariation WindSpeed(double min, double max)
    {
        if (min > max) throw new ArgumentException("min must be <= max.");
        WindSpeedMin = min;
        WindSpeedMax = max;
        return this;
    }

    /// <summary>Gaussian jitter on wind direction (std dev in degrees).</summary>
    public WildfireVariation WindDirectionJitter(double stdDevDegrees)
    {
        WindDirectionJitterDeg = stdDevDegrees;
        return this;
    }

    /// <summary>Uniform random dead fuel moisture content range (fraction).</summary>
    public WildfireVariation Moisture(double min, double max)
    {
        if (min > max) throw new ArgumentException("min must be <= max.");
        MoistureMin = min;
        MoistureMax = max;
        return this;
    }

    /// <summary>
    /// Random ignition point offset radius in cell units.
    /// Ignition point is shifted randomly within this radius.
    /// </summary>
    public WildfireVariation IgnitionOffset(double radiusCells)
    {
        IgnitionOffsetRadius = radiusCells;
        return this;
    }

    /// <summary>True when at least one variation range is configured.</summary>
    public bool HasVariation =>
        WindSpeedMin.HasValue ||
        (WindDirectionJitterDeg.HasValue && WindDirectionJitterDeg.Value > 0) ||
        MoistureMin.HasValue ||
        IgnitionOffsetRadius.HasValue;
}
