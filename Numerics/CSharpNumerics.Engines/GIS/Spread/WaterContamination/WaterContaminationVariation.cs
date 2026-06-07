using System;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Defines stochastic parameter variation ranges for Monte Carlo
/// water contamination simulations.
/// Each nullable property pair describes the range from which the MC runner will sample.
/// When null, the corresponding parameter is held at its deterministic baseline value.
/// </summary>
public class WaterContaminationVariation
{
    /// <summary>Minimum base discharge (m³/s) for sampling. Null = use baseline.</summary>
    public double? DischargeMin { get; private set; }
    /// <summary>Maximum base discharge (m³/s) for sampling.</summary>
    public double? DischargeMax { get; private set; }
    /// <summary>Minimum source concentration (mg/L) for sampling. Null = use baseline.</summary>
    public double? SourceConcentrationMin { get; private set; }
    /// <summary>Maximum source concentration (mg/L) for sampling.</summary>
    public double? SourceConcentrationMax { get; private set; }
    /// <summary>Minimum Manning's n roughness coefficient for sampling. Null = use baseline.</summary>
    public double? ManningNMin { get; private set; }
    /// <summary>Maximum Manning's n roughness coefficient for sampling.</summary>
    public double? ManningNMax { get; private set; }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent setters
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Uniform random baseline discharge range (m³/s).</summary>
    public WaterContaminationVariation Discharge(double min, double max)
    {
        if (min > max) throw new ArgumentException("min must be <= max.");
        DischargeMin = min;
        DischargeMax = max;
        return this;
    }

    /// <summary>Uniform random source concentration multiplier range (mg/L).</summary>
    public WaterContaminationVariation SourceConcentration(double min, double max)
    {
        if (min > max) throw new ArgumentException("min must be <= max.");
        SourceConcentrationMin = min;
        SourceConcentrationMax = max;
        return this;
    }

    /// <summary>Uniform random Manning's n range.</summary>
    public WaterContaminationVariation ManningN(double min, double max)
    {
        if (min > max) throw new ArgumentException("min must be <= max.");
        ManningNMin = min;
        ManningNMax = max;
        return this;
    }

    /// <summary>True when at least one variation range is configured.</summary>
    public bool HasVariation =>
        DischargeMin.HasValue ||
        SourceConcentrationMin.HasValue ||
        ManningNMin.HasValue;
}
