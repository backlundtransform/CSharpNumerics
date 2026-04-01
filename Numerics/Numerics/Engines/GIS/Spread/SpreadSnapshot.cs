using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread;

/// <summary>
/// A single time-step snapshot of a spread simulation. Wraps a <see cref="GridSnapshot"/>
/// and provides fire-specific named layers and convenience queries.
/// <para>
/// Standard layers: <c>"burnState"</c>, <c>"flameLength"</c>,
/// <c>"rateOfSpread"</c>, <c>"burnTime"</c>.
/// </para>
/// </summary>
public class SpreadSnapshot
{
    /// <summary>The underlying grid snapshot (primary values = burn state as double).</summary>
    public GridSnapshot Snapshot { get; }

    /// <summary>Simulation time in seconds.</summary>
    public double Time => Snapshot.Time;

    /// <summary>Zero-based time-step index.</summary>
    public int TimeIndex => Snapshot.TimeIndex;

    /// <summary>The grid this snapshot is defined over.</summary>
    public GeoGrid Grid => Snapshot.Grid;

    /// <summary>
    /// Creates a spread snapshot from pre-built arrays.
    /// All arrays must have length = grid.Nx × grid.Ny (2-D, Nz=1).
    /// </summary>
    public SpreadSnapshot(
        GeoGrid grid,
        double time,
        int timeIndex,
        double[] burnState,
        double[] flameLength,
        double[] rateOfSpread,
        double[] burnTime)
    {
        Snapshot = new GridSnapshot(grid, burnState, time, timeIndex);
        Snapshot.SetLayer("burnState", burnState);
        Snapshot.SetLayer("flameLength", flameLength);
        Snapshot.SetLayer("rateOfSpread", rateOfSpread);
        Snapshot.SetLayer("burnTime", burnTime);
    }

    /// <summary>Number of cells currently in <see cref="CellBurnState.Burning"/> state.</summary>
    public int BurningCellCount
    {
        get
        {
            var bs = Snapshot.GetLayer("burnState");
            int count = 0;
            for (int i = 0; i < bs.Length; i++)
                if ((int)bs[i] == (int)CellBurnState.Burning) count++;
            return count;
        }
    }

    /// <summary>Number of cells in <see cref="CellBurnState.Burned"/> state.</summary>
    public int BurnedCellCount
    {
        get
        {
            var bs = Snapshot.GetLayer("burnState");
            int count = 0;
            for (int i = 0; i < bs.Length; i++)
                if ((int)bs[i] == (int)CellBurnState.Burned) count++;
            return count;
        }
    }

    /// <summary>
    /// Total area (hectares) of cells that are Burning or Burned.
    /// </summary>
    public double BurnedAreaHectares
    {
        get
        {
            var bs = Snapshot.GetLayer("burnState");
            int affectedCount = 0;
            for (int i = 0; i < bs.Length; i++)
            {
                int state = (int)bs[i];
                if (state == (int)CellBurnState.Burning || state == (int)CellBurnState.Burned)
                    affectedCount++;
            }
            double cellAreaM2 = Grid.Step * Grid.Step;
            return affectedCount * cellAreaM2 / 10000.0;
        }
    }

    /// <summary>Gets the burn state for cell (ix, iy).</summary>
    public CellBurnState GetBurnState(int ix, int iy)
    {
        var bs = Snapshot.GetLayer("burnState");
        return (CellBurnState)(int)bs[iy * Grid.Nx + ix];
    }
}
