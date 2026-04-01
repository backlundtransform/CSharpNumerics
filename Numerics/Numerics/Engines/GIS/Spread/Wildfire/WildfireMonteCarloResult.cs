using CSharpNumerics.Engines.GIS.Grid;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Result of a Monte Carlo wildfire ensemble.
/// Contains per-cell burn probability and per-iteration burned area statistics.
/// </summary>
public class WildfireMonteCarloResult
{
    /// <summary>Number of Monte Carlo iterations.</summary>
    public int Iterations { get; }

    /// <summary>The spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>
    /// Per-cell burn probability (fraction 0–1).
    /// Length = Nx × Ny. Value = fraction of iterations where the cell burned.
    /// </summary>
    public double[] BurnProbability { get; }

    /// <summary>Burned area (hectares) per iteration.</summary>
    public double[] BurnedAreas { get; }

    /// <summary>Mean burned area across all iterations (hectares).</summary>
    public double MeanBurnedArea
    {
        get
        {
            double sum = 0;
            for (int i = 0; i < BurnedAreas.Length; i++) sum += BurnedAreas[i];
            return BurnedAreas.Length > 0 ? sum / BurnedAreas.Length : 0;
        }
    }

    /// <summary>Maximum burned area across all iterations (hectares).</summary>
    public double MaxBurnedArea
    {
        get
        {
            double max = 0;
            for (int i = 0; i < BurnedAreas.Length; i++)
                if (BurnedAreas[i] > max) max = BurnedAreas[i];
            return max;
        }
    }

    /// <summary>
    /// Per-iteration snapshot lists for optional downstream analysis (clustering).
    /// </summary>
    public IReadOnlyList<IReadOnlyList<SpreadSnapshot>> AllSnapshots { get; }

    /// <summary>
    /// Creates a Monte Carlo wildfire result.
    /// </summary>
    public WildfireMonteCarloResult(
        int iterations,
        GeoGrid grid,
        double[] burnProbability,
        double[] burnedAreas,
        IReadOnlyList<IReadOnlyList<SpreadSnapshot>> allSnapshots)
    {
        Iterations = iterations;
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        BurnProbability = burnProbability ?? throw new ArgumentNullException(nameof(burnProbability));
        BurnedAreas = burnedAreas ?? throw new ArgumentNullException(nameof(burnedAreas));
        AllSnapshots = allSnapshots;
    }
}
