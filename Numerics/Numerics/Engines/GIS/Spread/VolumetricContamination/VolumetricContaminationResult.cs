using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;

/// <summary>
/// Result of a single deterministic 3D volumetric contamination simulation.
/// Provides time series of snapshots and depth-aware analysis accessors.
/// </summary>
public class VolumetricContaminationResult
{
    /// <summary>Ordered list of per-time-step snapshots.</summary>
    public IReadOnlyList<SpreadSnapshot> Snapshots { get; }

    /// <summary>The spatial grid (Nx × Ny × Nz).</summary>
    public GeoGrid Grid { get; }

    /// <summary>True if a CFL stability violation was detected.</summary>
    public bool CflViolationDetected { get; }

    private readonly double _threshold;

    /// <summary>Creates a 3D volumetric contamination result.</summary>
    public VolumetricContaminationResult(
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

    /// <summary>Number of cells contaminated above threshold at the final time step.</summary>
    public int AffectedCellCount
    {
        get
        {
            if (Snapshots.Count == 0) return 0;
            return Snapshots[Snapshots.Count - 1].ContaminatedCellCount(_threshold);
        }
    }

    /// <summary>
    /// Per-cell peak concentration over the entire time series.
    /// Flat array of length Nx × Ny × Nz in grid index order.
    /// </summary>
    public double[] PeakConcentration3D()
    {
        int cellCount = Grid.CellCount;
        var peak = new double[cellCount];
        foreach (var snap in Snapshots)
        {
            var conc = snap.Snapshot.GetLayer("concentration");
            for (int i = 0; i < cellCount; i++)
                if (conc[i] > peak[i]) peak[i] = conc[i];
        }
        return peak;
    }

    /// <summary>
    /// Concentration-vs-depth profile at horizontal position (ix, iy)
    /// for a given time step.
    /// </summary>
    public DepthProfile GetDepthProfile(int ix, int iy, int timeIndex = -1)
    {
        if (ix < 0 || ix >= Grid.Nx) throw new ArgumentOutOfRangeException(nameof(ix));
        if (iy < 0 || iy >= Grid.Ny) throw new ArgumentOutOfRangeException(nameof(iy));

        int ti = timeIndex < 0 ? Snapshots.Count - 1 : timeIndex;
        if (ti < 0 || ti >= Snapshots.Count)
            throw new ArgumentOutOfRangeException(nameof(timeIndex));

        int nz = Grid.Nz;
        var depths = new double[nz];
        var concs = new double[nz];
        var conc = Snapshots[ti].Snapshot.GetLayer("concentration");

        for (int iz = 0; iz < nz; iz++)
        {
            depths[iz] = Grid.ZMin + iz * Grid.Step;
            concs[iz] = conc[Grid.Index(ix, iy, iz)];
        }

        return new DepthProfile(depths, concs);
    }

    /// <summary>
    /// Surface (iz = Nz − 1) concentration slice at a given time step.
    /// Returns flat array of length Nx × Ny.
    /// </summary>
    public double[] SurfaceSlice(int timeIndex)
    {
        return HorizontalSlice(timeIndex, Grid.Nz - 1);
    }

    /// <summary>
    /// Bottom (iz = 0) concentration slice at a given time step.
    /// Returns flat array of length Nx × Ny.
    /// </summary>
    public double[] BottomSlice(int timeIndex)
    {
        return HorizontalSlice(timeIndex, 0);
    }

    /// <summary>
    /// Horizontal concentration slice at a given depth layer and time step.
    /// Returns flat array of length Nx × Ny in row-major order (ix fastest).
    /// </summary>
    public double[] HorizontalSlice(int timeIndex, int iz)
    {
        if (timeIndex < 0 || timeIndex >= Snapshots.Count)
            throw new ArgumentOutOfRangeException(nameof(timeIndex));
        if (iz < 0 || iz >= Grid.Nz)
            throw new ArgumentOutOfRangeException(nameof(iz));

        int nx = Grid.Nx;
        int ny = Grid.Ny;
        var slice = new double[nx * ny];
        var conc = Snapshots[timeIndex].Snapshot.GetLayer("concentration");

        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
                slice[iy * nx + ix] = conc[Grid.Index(ix, iy, iz)];

        return slice;
    }

    /// <summary>
    /// Deepest z-coordinate (metres) where peak concentration exceeds the
    /// given threshold at any time step. Returns <see cref="Grid.ZMin"/>
    /// if no cell exceeds the threshold.
    /// </summary>
    public double MaxAffectedDepth(double threshold)
    {
        var peak = PeakConcentration3D();
        int nx = Grid.Nx;
        int ny = Grid.Ny;
        int nz = Grid.Nz;
        int deepestIz = -1;

        for (int iz = 0; iz < nz; iz++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = Grid.Index(ix, iy, iz);
                    if (peak[idx] > threshold && iz > deepestIz)
                        deepestIz = iz;
                }

        if (deepestIz < 0) return Grid.ZMin;
        return Grid.ZMin + deepestIz * Grid.Step;
    }

    /// <summary>
    /// Generates a contamination extent polygon for a horizontal slice
    /// at the surface layer at a given time step.
    /// </summary>
    public ExposurePolygon GenerateContaminationExtent(int timeIndex, double threshold)
    {
        if (timeIndex < 0 || timeIndex >= Snapshots.Count)
            throw new ArgumentOutOfRangeException(nameof(timeIndex));

        // Build a 2D GridSnapshot from the surface slice
        int nx = Grid.Nx;
        int ny = Grid.Ny;
        int surfaceIz = Grid.Nz - 1;
        var surfaceConc = new double[nx * ny];
        var conc = Snapshots[timeIndex].Snapshot.GetLayer("concentration");

        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
                surfaceConc[iy * nx + ix] = conc[Grid.Index(ix, iy, surfaceIz)];

        // Create a 2D grid for the polygon generator
        var grid2d = new GeoGrid(Grid.XMin, Grid.XMax, Grid.YMin, Grid.YMax, 0, 0, Grid.Step);
        var gs = new GridSnapshot(grid2d, surfaceConc, Snapshots[timeIndex].Time, timeIndex);
        gs.SetLayer("concentration", surfaceConc);

        var gridSnapshots = new List<GridSnapshot> { gs };
        return ExposurePolygonGenerator.PeakExposure(gridSnapshots, threshold, "concentration");
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
}
