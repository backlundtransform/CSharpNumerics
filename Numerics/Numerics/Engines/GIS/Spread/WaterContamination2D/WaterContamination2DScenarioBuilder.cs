using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Materials.Water;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination2D;

/// <summary>
/// Fluent builder for 2D water contamination spread scenarios on open water bodies.
/// <para>
/// Usage:
/// <code>
/// RiskScenario
///     .ForWaterContamination2D()
///     .WithSource(5, 5, 100, double.MaxValue)
///     .WithContaminant(ContaminantLibrary.Get("Benzene"))
///     .WithDiffusion(1.0, 1.0)
///     .WithCurrentField((ix, iy) => (0.1, 0))
///     .OverGrid(grid)
///     .OverTime(0, 86400, 60)
///     .RunSingle();
/// </code>
/// </para>
/// </summary>
public class WaterContamination2DScenarioBuilder
{
    private readonly List<(int ix, int iy, double concentrationMgL, double durationSeconds)> _sources = new();
    private AquaticContaminant _contaminant;
    private double _diffusionX = 1.0;
    private double _diffusionY = -1;
    private Func<int, int, (double vx, double vy)> _velocityField;
    private bool[] _landMask;
    private GeoGrid _grid;
    private TerrainGrid _terrain;
    private TimeFrame _timeFrame;

    internal WaterContamination2DScenarioBuilder() { }

    /// <summary>Adds a contamination point-source.</summary>
    public WaterContamination2DScenarioBuilder WithSource(
        int ix, int iy, double concentrationMgL, double durationSeconds)
    {
        _sources.Add((ix, iy, concentrationMgL, durationSeconds));
        return this;
    }

    /// <summary>Sets the contaminant descriptor.</summary>
    public WaterContamination2DScenarioBuilder WithContaminant(AquaticContaminant contaminant)
    {
        _contaminant = contaminant;
        return this;
    }

    /// <summary>
    /// Sets the diffusion coefficients (m²/s). If only Dx is given, Dy = Dx (isotropic).
    /// </summary>
    public WaterContamination2DScenarioBuilder WithDiffusion(double dx, double dy = -1)
    {
        _diffusionX = dx;
        _diffusionY = dy;
        return this;
    }

    /// <summary>
    /// Sets a prescribed 2D current velocity field. The function receives grid indices (ix, iy)
    /// and should return velocity (vx, vy) in m/s.
    /// </summary>
    public WaterContamination2DScenarioBuilder WithCurrentField(
        Func<int, int, (double vx, double vy)> field)
    {
        _velocityField = field;
        return this;
    }

    /// <summary>
    /// Sets a uniform current velocity (m/s) across the entire domain.
    /// </summary>
    public WaterContamination2DScenarioBuilder WithCurrent(double vx, double vy)
    {
        _velocityField = (_, __) => (vx, vy);
        return this;
    }

    /// <summary>
    /// Sets a land mask. Cells where landMask[iy * Nx + ix] = true are treated
    /// as impermeable barriers.
    /// </summary>
    public WaterContamination2DScenarioBuilder WithLandMask(bool[] landMask)
    {
        _landMask = landMask;
        return this;
    }

    /// <summary>Sets the terrain (optional, unused by 2D solver but kept for API consistency).</summary>
    public WaterContamination2DScenarioBuilder WithTerrain(TerrainGrid terrain)
    {
        _terrain = terrain;
        return this;
    }

    /// <summary>Sets the spatial grid.</summary>
    public WaterContamination2DScenarioBuilder OverGrid(GeoGrid grid)
    {
        _grid = grid ?? throw new ArgumentNullException(nameof(grid));
        return this;
    }

    /// <summary>Sets the simulation time frame.</summary>
    public WaterContamination2DScenarioBuilder OverTime(
        double startSeconds, double endSeconds, double stepSeconds)
    {
        _timeFrame = new TimeFrame(startSeconds, endSeconds, stepSeconds);
        return this;
    }

    /// <summary>Runs a single deterministic 2D water contamination simulation.</summary>
    public WaterContamination2DResult RunSingle()
    {
        Validate();

        var parameters = new WaterContamination2DParameters(
            _sources.AsReadOnly(), _contaminant,
            _diffusionX, _diffusionY,
            _velocityField, _landMask);

        var simulator = new WaterContamination2DSimulator(parameters);
        var terrain = _terrain ?? TerrainGrid.FromFunction(_grid, (x, y) => 0);
        var snapshots = simulator.Run(_grid, terrain, new FuelMap(_grid), _timeFrame);

        return new WaterContamination2DResult(
            snapshots, _grid,
            _contaminant.ToxicityThresholdMgL,
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
        if (_contaminant.Name == null)
            throw new InvalidOperationException("Contaminant not set. Call WithContaminant().");
    }
}
