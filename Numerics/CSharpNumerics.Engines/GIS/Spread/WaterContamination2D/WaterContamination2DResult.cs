using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination2D;

/// <summary>
/// Result of a single deterministic 2D water contamination simulation.
/// Contains the time series of <see cref="SpreadSnapshot"/> frames and
/// convenience accessors for contamination-specific metrics.
/// </summary>
public class WaterContamination2DResult
{
    /// <summary>Ordered list of per-time-step snapshots.</summary>
    public IReadOnlyList<SpreadSnapshot> Snapshots { get; }

    /// <summary>The spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>Peak concentration (mg/L) across all cells and time steps.</summary>
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

    /// <summary>Time (seconds) at which maximum concentration first occurs.</summary>
    public double PeakArrivalTimeSeconds
    {
        get
        {
            double peak = 0, arrival = 0;
            foreach (var snap in Snapshots)
            {
                double m = snap.MaxConcentration;
                if (m > peak) { peak = m; arrival = snap.Time; }
            }
            return arrival;
        }
    }

    /// <summary>Number of cells contaminated above the threshold at the final time step.</summary>
    public int AffectedCellCount
    {
        get
        {
            if (Snapshots.Count == 0) return 0;
            return Snapshots[Snapshots.Count - 1].ContaminatedCellCount(_threshold);
        }
    }

    /// <summary>Area (m²) of cells contaminated above threshold at the final time step.</summary>
    public double AffectedAreaM2 => AffectedCellCount * Grid.Step * Grid.Step;

    /// <summary>True if a CFL stability violation was detected.</summary>
    public bool CflViolationDetected { get; }

    private readonly double _threshold;

    /// <summary>Creates a 2D water contamination result.</summary>
    public WaterContamination2DResult(
        IReadOnlyList<SpreadSnapshot> snapshots,
        GeoGrid grid,
        double toxicityThreshold = 0.01,
        bool cflViolationDetected = false)
    {
        Snapshots = snapshots ?? throw new ArgumentNullException(nameof(snapshots));
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        _threshold = toxicityThreshold;
        CflViolationDetected = cflViolationDetected;
    }

    /// <summary>
    /// Total time (seconds) that any cell exceeds the given threshold,
    /// summed across all cells and time steps.
    /// </summary>
    public double ExceedanceDurationSeconds(double threshold)
    {
        if (Snapshots.Count < 2) return 0;
        double dt = Snapshots[1].Time - Snapshots[0].Time;
        double total = 0;
        foreach (var snap in Snapshots)
            total += snap.ContaminatedCellCount(threshold) * dt;
        return total;
    }

    /// <summary>
    /// Generates an exposure polygon for the contaminated area at a given time step.
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
