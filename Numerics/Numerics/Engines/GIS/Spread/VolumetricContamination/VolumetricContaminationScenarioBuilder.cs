using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Terrain;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;

/// <summary>
/// Fluent builder for 3D volumetric water contamination spread scenarios.
/// <para>
/// Usage:
/// <code>
/// RiskScenario
///     .ForVolumetricContamination()
///     .WithSource(5, 5, 3, 100, double.MaxValue)
///     .WithDiffusivity(1.0, 0.01)
///     .WithCurrent(0.1, 0, 0)
///     .OverGrid(grid)     // Must have Nz &gt; 1
///     .OverTime(0, 86400, 60)
///     .RunSingle();
/// </code>
/// </para>
/// </summary>
public class VolumetricContaminationScenarioBuilder
{
    private readonly List<(int ix, int iy, int iz, double concentrationMgL, double durationSeconds)> _sources = new();
    private double _horizontalDiffusivity = 1.0;
    private double _verticalDiffusivity = 0.01;
    private Func<int, int, int, (double vx, double vy, double vz)> _velocityField;
    private double _decayRate;
    private bool[] _landMask;
    private double _toxicityThreshold = 0.01;
    private GeoGrid _grid;
    private TerrainGrid _terrain;
    private TimeFrame _timeFrame;

    internal VolumetricContaminationScenarioBuilder() { }

    /// <summary>Adds a contamination point-source at (ix, iy, iz).</summary>
    public VolumetricContaminationScenarioBuilder WithSource(
        int ix, int iy, int iz, double concentrationMgL, double durationSeconds)
    {
        _sources.Add((ix, iy, iz, concentrationMgL, durationSeconds));
        return this;
    }

    /// <summary>
    /// Sets the diffusivity coefficients (m²/s).
    /// Horizontal applies to x and y; vertical applies to z.
    /// </summary>
    public VolumetricContaminationScenarioBuilder WithDiffusivity(
        double horizontal, double vertical = 0.01)
    {
        _horizontalDiffusivity = horizontal;
        _verticalDiffusivity = vertical;
        return this;
    }

    /// <summary>
    /// Sets a prescribed 3D current velocity field. The function receives grid
    /// indices (ix, iy, iz) and should return velocity (vx, vy, vz) in m/s.
    /// </summary>
    public VolumetricContaminationScenarioBuilder WithCurrentField(
        Func<int, int, int, (double vx, double vy, double vz)> field)
    {
        _velocityField = field;
        return this;
    }

    /// <summary>Sets a uniform current velocity (m/s) across the entire domain.</summary>
    public VolumetricContaminationScenarioBuilder WithCurrent(double vx, double vy, double vz)
    {
        _velocityField = (_, __, ___) => (vx, vy, vz);
        return this;
    }

    /// <summary>Sets the first-order decay rate (s⁻¹). Zero = conservative tracer.</summary>
    public VolumetricContaminationScenarioBuilder WithDecay(double rateSec)
    {
        _decayRate = rateSec;
        return this;
    }

    /// <summary>Directly sets a land/seabed mask array.</summary>
    public VolumetricContaminationScenarioBuilder WithLandMask(bool[] landMask)
    {
        _landMask = landMask;
        return this;
    }

    /// <summary>
    /// Generates a land mask from terrain: cells where elevation &gt; seaLevel are land.
    /// The mask is computed per (ix, iy) column and applied to all depths.
    /// </summary>
    public VolumetricContaminationScenarioBuilder WithLandMaskFromTerrain(
        TerrainGrid terrain, double seaLevel)
    {
        if (_grid == null)
            throw new InvalidOperationException("Call OverGrid() before WithLandMaskFromTerrain().");

        int nx = _grid.Nx;
        int ny = _grid.Ny;
        int nz = _grid.Nz;
        int nxy = nx * ny;
        _landMask = new bool[_grid.CellCount];

        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                double elev = terrain[ix, iy];
                if (elev > seaLevel)
                {
                    for (int iz = 0; iz < nz; iz++)
                        _landMask[iz * nxy + iy * nx + ix] = true;
                }
            }

        return this;
    }

    /// <summary>Sets the toxicity threshold (mg/L) for result metrics.</summary>
    public VolumetricContaminationScenarioBuilder WithToxicityThreshold(double threshold)
    {
        _toxicityThreshold = threshold;
        return this;
    }

    /// <summary>Sets the terrain (optional).</summary>
    public VolumetricContaminationScenarioBuilder WithTerrain(TerrainGrid terrain)
    {
        _terrain = terrain;
        return this;
    }

    /// <summary>Sets the spatial grid. Must have Nz &gt; 1 for volumetric simulation.</summary>
    public VolumetricContaminationScenarioBuilder OverGrid(GeoGrid grid)
    {
        _grid = grid ?? throw new ArgumentNullException(nameof(grid));
        return this;
    }

    /// <summary>Sets the simulation time frame.</summary>
    public VolumetricContaminationScenarioBuilder OverTime(
        double startSeconds, double endSeconds, double stepSeconds)
    {
        _timeFrame = new TimeFrame(startSeconds, endSeconds, stepSeconds);
        return this;
    }

    /// <summary>Runs a single deterministic 3D volumetric contamination simulation.</summary>
    public VolumetricContaminationResult RunSingle()
    {
        Validate();

        var parameters = new VolumetricContaminationParameters(
            _sources.AsReadOnly(),
            _horizontalDiffusivity,
            _verticalDiffusivity,
            _velocityField,
            _decayRate,
            _landMask,
            _toxicityThreshold);

        var simulator = new VolumetricContaminationSimulator(parameters);
        var terrain = _terrain ?? TerrainGrid.FromFunction(_grid, (x, y) => 0);
        var snapshots = simulator.Run(_grid, terrain, new FuelMap(_grid), _timeFrame);

        return new VolumetricContaminationResult(
            snapshots, _grid, _toxicityThreshold,
            simulator.CflViolationDetected);
    }

    private void Validate()
    {
        if (_grid == null)
            throw new InvalidOperationException("Grid not set. Call OverGrid().");
        if (_timeFrame == null)
            throw new InvalidOperationException("Time frame not set. Call OverTime().");
        if (_sources.Count == 0)
            throw new InvalidOperationException("No sources defined. Call WithSource().");
    }
}
