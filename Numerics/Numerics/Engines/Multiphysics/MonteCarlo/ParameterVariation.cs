using CSharpNumerics.Statistics.Random;
using System;

namespace CSharpNumerics.Engines.Multiphysics.MonteCarlo;

/// <summary>
/// Defines ranges for stochastic variation of simulation parameters.
/// Used by <see cref="MultiphysicsMonteCarloModel"/> to sample from
/// uniform or Gaussian distributions per trial.
/// <para>
/// Usage:
/// <code>
/// var variation = new ParameterVariation()
///     .ThermalConductivity(40, 60)
///     .BoundaryTemperature(90, 110)
///     .Load(800, 1200);
/// </code>
/// </para>
/// </summary>
public class ParameterVariation
{
    // ── Material property ranges ─────────────────────────────────

    /// <summary>Thermal conductivity k range [W/(m·K)].</summary>
    public double? ThermalConductivityMin { get; private set; }
    public double? ThermalConductivityMax { get; private set; }

    /// <summary>Young's modulus E range [Pa].</summary>
    public double? YoungsModulusMin { get; private set; }
    public double? YoungsModulusMax { get; private set; }

    /// <summary>Dynamic viscosity μ range [Pa·s].</summary>
    public double? ViscosityMin { get; private set; }
    public double? ViscosityMax { get; private set; }

    /// <summary>Density ρ range [kg/m³].</summary>
    public double? DensityMin { get; private set; }
    public double? DensityMax { get; private set; }

    // ── Boundary condition ranges ────────────────────────────────

    /// <summary>Top/primary BC value range.</summary>
    public double? BoundaryValueMin { get; private set; }
    public double? BoundaryValueMax { get; private set; }

    /// <summary>Pressure gradient dP/dx range [Pa/m].</summary>
    public double? PressureGradientMin { get; private set; }
    public double? PressureGradientMax { get; private set; }

    // ── Load ranges ──────────────────────────────────────────────

    /// <summary>Point or distributed load range [N or N/m].</summary>
    public double? LoadMin { get; private set; }
    public double? LoadMax { get; private set; }

    /// <summary>Source intensity range [W/m³ or C/m³].</summary>
    public double? SourceIntensityMin { get; private set; }
    public double? SourceIntensityMax { get; private set; }

    /// <summary>Whether any variation is configured.</summary>
    public bool HasVariation =>
        ThermalConductivityMin.HasValue || YoungsModulusMin.HasValue ||
        ViscosityMin.HasValue || DensityMin.HasValue ||
        BoundaryValueMin.HasValue || PressureGradientMin.HasValue ||
        LoadMin.HasValue || SourceIntensityMin.HasValue;

    // ═══════════════════════════════════════════════════════════════
    //  Fluent setters
    // ═══════════════════════════════════════════════════════════════

    public ParameterVariation ThermalConductivity(double min, double max)
    {
        ThermalConductivityMin = min; ThermalConductivityMax = max;
        return this;
    }

    public ParameterVariation YoungsModulus(double min, double max)
    {
        YoungsModulusMin = min; YoungsModulusMax = max;
        return this;
    }

    public ParameterVariation Viscosity(double min, double max)
    {
        ViscosityMin = min; ViscosityMax = max;
        return this;
    }

    public ParameterVariation Density(double min, double max)
    {
        DensityMin = min; DensityMax = max;
        return this;
    }

    public ParameterVariation BoundaryTemperature(double min, double max)
    {
        BoundaryValueMin = min; BoundaryValueMax = max;
        return this;
    }

    public ParameterVariation PressureGradient(double min, double max)
    {
        PressureGradientMin = min; PressureGradientMax = max;
        return this;
    }

    public ParameterVariation Load(double min, double max)
    {
        LoadMin = min; LoadMax = max;
        return this;
    }

    public ParameterVariation SourceIntensity(double min, double max)
    {
        SourceIntensityMin = min; SourceIntensityMax = max;
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Sampling
    // ═══════════════════════════════════════════════════════════════

    internal double SampleOrDefault(RandomGenerator rng, double? min, double? max, double defaultValue)
    {
        if (min.HasValue && max.HasValue)
            return rng.NextUniform(min.Value, max.Value);
        return defaultValue;
    }
}
