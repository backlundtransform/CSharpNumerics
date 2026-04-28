using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Result of a single deterministic wildfire simulation.
/// Contains the time series of <see cref="SpreadSnapshot"/> frames and
/// convenience accessors for fire-specific metrics.
/// </summary>
public class WildfireScenarioResult
{
    /// <summary>Ordered list of per-time-step snapshots.</summary>
    public IReadOnlyList<SpreadSnapshot> Snapshots { get; }

    /// <summary>The spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>
    /// Total burned area (hectares) at the final time step.
    /// Includes cells in both Burning and Burned state.
    /// </summary>
    public double FinalBurnedArea => Snapshots.Count > 0
        ? Snapshots[Snapshots.Count - 1].BurnedAreaHectares
        : 0;

    /// <summary>
    /// Maximum flame length (metres) observed across all cells and time steps.
    /// </summary>
    public double MaxFlameLength
    {
        get
        {
            double max = 0;
            foreach (var snap in Snapshots)
            {
                var fl = snap.Snapshot.GetLayer("flameLength");
                for (int i = 0; i < fl.Length; i++)
                    if (fl[i] > max) max = fl[i];
            }
            return max;
        }
    }

    /// <summary>
    /// Creates a wildfire scenario result from a list of snapshots.
    /// </summary>
    public WildfireScenarioResult(IReadOnlyList<SpreadSnapshot> snapshots, GeoGrid grid)
    {
        Snapshots = snapshots ?? throw new ArgumentNullException(nameof(snapshots));
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
    }

    /// <summary>
    /// Generates a fire perimeter polygon at the given time-step index.
    /// Uses marching squares on the burnState layer (threshold = 0.5 to capture
    /// Burning + Burned cells). Firebreak cells are excluded so that the
    /// perimeter does not extend over non-burnable areas such as water.
    /// </summary>
    public ExposurePolygon GenerateFirePerimeter(int timeIndex)
    {
        if (timeIndex < 0 || timeIndex >= Snapshots.Count)
            throw new ArgumentOutOfRangeException(nameof(timeIndex));

        var snap = Snapshots[timeIndex];
        var gridSnapshots = new List<GridSnapshot> { snap.Snapshot };
        var exclude = new HashSet<double> { (double)CellBurnState.Firebreak };
        return ExposurePolygonGenerator.PeakExposure(gridSnapshots, 0.5, "burnState", exclude);
    }

    /// <summary>
    /// Generates fire perimeter polygons for all time steps.
    /// </summary>
    public IReadOnlyList<ExposurePolygon> FirePerimeters
    {
        get
        {
            var perimeters = new List<ExposurePolygon>();
            for (int i = 0; i < Snapshots.Count; i++)
                perimeters.Add(GenerateFirePerimeter(i));
            return perimeters.AsReadOnly();
        }
    }
}
