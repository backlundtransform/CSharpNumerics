using CSharpNumerics.Physics.Materials.Water;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Runtime configuration for a water contamination spread simulation.
/// </summary>
public class WaterContaminationParameters
{
    /// <summary>
    /// Point-source injection sites. Each tuple specifies: cell (ix, iy),
    /// injection concentration (mg/L), and duration of injection (seconds).
    /// After the duration expires the source stops injecting.
    /// </summary>
    public IReadOnlyList<(int ix, int iy, double concentrationMgL, double durationSeconds)> Sources { get; }

    /// <summary>
    /// Contaminant descriptor (decay half-life, partition coefficient, toxicity thresholds).
    /// </summary>
    public AquaticContaminant Contaminant { get; }

    /// <summary>
    /// Reference river discharge Q (m³/s). Used for tributary mixing mass balance.
    /// </summary>
    public double BaseDischargeM3s { get; }

    /// <summary>
    /// Bed sediment porosity n (fraction, 0–1). Used for retardation factor.
    /// Default 0.4.
    /// </summary>
    public double BedPorosity { get; }

    /// <summary>
    /// Bulk sediment density ρs (kg/m³). Used for retardation factor.
    /// Default 1600.
    /// </summary>
    public double BedBulkDensity { get; }

    /// <summary>
    /// Creates water contamination simulation parameters.
    /// </summary>
    public WaterContaminationParameters(
        IReadOnlyList<(int ix, int iy, double concentrationMgL, double durationSeconds)> sources,
        AquaticContaminant contaminant,
        double baseDischargeM3s = 10.0,
        double bedPorosity = 0.4,
        double bedBulkDensity = 1600.0)
    {
        if (sources == null || sources.Count == 0)
            throw new ArgumentException("At least one contamination source is required.", nameof(sources));

        Sources = sources;
        Contaminant = contaminant;
        BaseDischargeM3s = baseDischargeM3s;
        BedPorosity = bedPorosity;
        BedBulkDensity = bedBulkDensity;
    }
}
