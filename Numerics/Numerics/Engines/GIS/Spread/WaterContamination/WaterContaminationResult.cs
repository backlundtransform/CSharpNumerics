using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Result of a single deterministic water contamination simulation.
/// Contains the time series of <see cref="SpreadSnapshot"/> frames and
/// convenience accessors for contamination-specific metrics.
/// </summary>
public class WaterContaminationResult
{
    /// <summary>Ordered list of per-time-step snapshots.</summary>
    public IReadOnlyList<SpreadSnapshot> Snapshots { get; }

    /// <summary>The spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>
    /// Maximum concentration (mg/L) observed across all cells and time steps.
    /// </summary>
    public double MaxConcentration
    {
        get
        {
            double max = 0;
            foreach (var snap in Snapshots)
            {
                double m = snap.MaxConcentration;
                if (m > max) max = m;
            }
            return max;
        }
    }

    /// <summary>
    /// Time (seconds) at which the maximum concentration first appears at
    /// the downstream-most contaminated cell.
    /// </summary>
    public double PeakArrivalTimeSeconds
    {
        get
        {
            double peak = 0;
            double arrivalTime = 0;
            foreach (var snap in Snapshots)
            {
                double m = snap.MaxConcentration;
                if (m > peak)
                {
                    peak = m;
                    arrivalTime = snap.Time;
                }
            }
            return arrivalTime;
        }
    }

    /// <summary>
    /// Total length (km) of river reach affected above the toxicity threshold
    /// at the final time step.
    /// </summary>
    public double TotalAffectedReachKm
    {
        get
        {
            if (Snapshots.Count == 0) return 0;
            return Snapshots[Snapshots.Count - 1].AffectedReachLengthKm(_threshold);
        }
    }

    private readonly double _threshold;

    /// <summary>
    /// Creates a water contamination result from a list of snapshots.
    /// </summary>
    /// <param name="snapshots">Per-time-step snapshots.</param>
    /// <param name="grid">The spatial grid.</param>
    /// <param name="toxicityThreshold">Concentration threshold for affected-reach metrics.</param>
    public WaterContaminationResult(
        IReadOnlyList<SpreadSnapshot> snapshots,
        GeoGrid grid,
        double toxicityThreshold = 0.01)
    {
        Snapshots = snapshots ?? throw new ArgumentNullException(nameof(snapshots));
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        _threshold = toxicityThreshold;
    }

    /// <summary>
    /// Total time (seconds) that any cell exceeds the given threshold.
    /// Summed across all cells and time steps, weighted by step duration.
    /// </summary>
    public double ExceedanceDurationSeconds(double threshold)
    {
        if (Snapshots.Count < 2) return 0;
        double dt = Snapshots.Count >= 2
            ? Snapshots[1].Time - Snapshots[0].Time
            : 0;

        double total = 0;
        foreach (var snap in Snapshots)
        {
            var conc = snap.Snapshot.GetLayer("concentration");
            for (int i = 0; i < conc.Length; i++)
            {
                if (conc[i] > threshold)
                {
                    total += dt;
                    break; // count per time step, not per cell
                }
            }
        }
        return total;
    }

    /// <summary>
    /// Generates a contamination extent polygon at the given time-step index
    /// using the concentration layer.
    /// </summary>
    public ExposurePolygon GenerateContaminationExtent(int timeIndex, double threshold)
    {
        if (timeIndex < 0 || timeIndex >= Snapshots.Count)
            throw new ArgumentOutOfRangeException(nameof(timeIndex));

        var snap = Snapshots[timeIndex];
        var gridSnapshots = new List<GridSnapshot> { snap.Snapshot };
        return ExposurePolygonGenerator.PeakExposure(gridSnapshots, threshold, "concentration");
    }
}
