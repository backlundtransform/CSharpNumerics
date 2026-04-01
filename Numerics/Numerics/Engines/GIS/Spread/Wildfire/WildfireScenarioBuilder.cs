using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using CSharpNumerics.Statistics.Random;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Fluent builder for wildfire spread scenarios. Follows the same pattern
/// as <c>RiskScenarioBuilder</c> (plume). Supports single deterministic runs
/// and Monte Carlo ensembles.
/// <para>
/// Usage:
/// <code>
/// RiskScenario
///     .ForWildfire()
///     .WithTerrain(terrain)
///     .WithFuel(fuelMap)
///     .WithIgnition(cx, cy)
///     .WithWind(5.0, new Vector(1, 0, 0))
///     .WithMoisture(0.08)
///     .OverGrid(grid)
///     .OverTime(0, 7200, 60)
///     .RunSingle();
/// </code>
/// </para>
/// </summary>
public class WildfireScenarioBuilder
{
    // ── Terrain & fuel ───────────────────────────────────────────
    private TerrainGrid _terrain;
    private FuelMap _fuelMap;

    // ── Ignition ─────────────────────────────────────────────────
    private readonly List<(int ix, int iy)> _ignitionPoints = new List<(int, int)>();

    // ── Weather ──────────────────────────────────────────────────
    private double _windSpeed;
    private Vector _windDirection;
    private double _burnDuration = 600;

    // ── Variation ────────────────────────────────────────────────
    private WildfireVariation _variation;

    // ── Grid & time ──────────────────────────────────────────────
    private GeoGrid _grid;
    private TimeFrame _timeFrame;

    internal WildfireScenarioBuilder() { }

    // ═══════════════════════════════════════════════════════════════
    //  Configuration steps (return self)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Sets the terrain elevation grid.</summary>
    public WildfireScenarioBuilder WithTerrain(TerrainGrid terrain)
    {
        _terrain = terrain ?? throw new ArgumentNullException(nameof(terrain));
        return this;
    }

    /// <summary>Sets the fuel map (per-cell fuel model + moisture).</summary>
    public WildfireScenarioBuilder WithFuel(FuelMap fuelMap)
    {
        _fuelMap = fuelMap ?? throw new ArgumentNullException(nameof(fuelMap));
        return this;
    }

    /// <summary>Adds an ignition point at cell (ix, iy).</summary>
    public WildfireScenarioBuilder WithIgnition(int ix, int iy)
    {
        _ignitionPoints.Add((ix, iy));
        return this;
    }

    /// <summary>Sets mid-flame wind speed and direction.</summary>
    public WildfireScenarioBuilder WithWind(double speed, Vector direction)
    {
        _windSpeed = speed;
        _windDirection = direction;
        return this;
    }

    /// <summary>
    /// Sets uniform dead fuel moisture content (fraction) across the entire fuel map.
    /// </summary>
    public WildfireScenarioBuilder WithMoisture(double moisture)
    {
        if (_fuelMap == null) throw new InvalidOperationException("Call WithFuel() before WithMoisture().");
        _fuelMap.SetUniformMoisture(moisture);
        return this;
    }

    /// <summary>Sets how long (seconds) a cell stays in Burning state. Default 600.</summary>
    public WildfireScenarioBuilder WithBurnDuration(double seconds)
    {
        _burnDuration = seconds;
        return this;
    }

    /// <summary>Configures stochastic variation for Monte Carlo runs.</summary>
    public WildfireScenarioBuilder WithVariation(Action<WildfireVariation> configure)
    {
        _variation = new WildfireVariation();
        configure(_variation);
        return this;
    }

    /// <summary>Sets the spatial grid.</summary>
    public WildfireScenarioBuilder OverGrid(GeoGrid grid)
    {
        _grid = grid ?? throw new ArgumentNullException(nameof(grid));
        return this;
    }

    /// <summary>Sets the simulation time frame.</summary>
    public WildfireScenarioBuilder OverTime(double startSeconds, double endSeconds, double stepSeconds)
    {
        _timeFrame = new TimeFrame(startSeconds, endSeconds, stepSeconds);
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Terminal: Run single deterministic scenario
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs a single deterministic wildfire simulation.
    /// </summary>
    public WildfireScenarioResult RunSingle()
    {
        Validate();

        var parameters = new WildfireParameters(
            _ignitionPoints.AsReadOnly(),
            _windSpeed,
            _windDirection,
            _burnDuration);

        var simulator = new WildfireSimulator(parameters);
        var snapshots = simulator.Run(_grid, _terrain, _fuelMap, _timeFrame);

        return new WildfireScenarioResult(snapshots, _grid);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Terminal: Run Monte Carlo ensemble
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs a Monte Carlo ensemble of wildfire simulations with stochastic
    /// parameter variation. Returns per-cell burn probability and area statistics.
    /// </summary>
    public WildfireMonteCarloResult RunMonteCarlo(int iterations, int? seed = null)
    {
        Validate();
        if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

        var rng = seed.HasValue ? new RandomGenerator(seed.Value) : new RandomGenerator();
        int cellCount = _grid.Nx * _grid.Ny;
        var burnCount = new int[cellCount]; // per-cell count of iterations where cell burned
        var burnedAreas = new double[iterations];
        var allSnapshots = new List<IReadOnlyList<SpreadSnapshot>>(iterations);

        for (int i = 0; i < iterations; i++)
        {
            // Sample parameters
            double windSpeed = _windSpeed;
            Vector windDir = _windDirection;
            double moisture = -1; // sentinel: no override
            var ignitionPoints = new List<(int ix, int iy)>(_ignitionPoints);

            if (_variation != null && _variation.HasVariation)
            {
                if (_variation.WindSpeedMin.HasValue)
                    windSpeed = rng.NextUniform(_variation.WindSpeedMin.Value, _variation.WindSpeedMax.Value);

                if (_variation.WindDirectionJitterDeg.HasValue && _variation.WindDirectionJitterDeg.Value > 0)
                {
                    double jitterRad = rng.NextGaussian(0, _variation.WindDirectionJitterDeg.Value) * Math.PI / 180.0;
                    double cos = Math.Cos(jitterRad);
                    double sin = Math.Sin(jitterRad);
                    windDir = new Vector(
                        windDir.x * cos - windDir.y * sin,
                        windDir.x * sin + windDir.y * cos,
                        0);
                }

                if (_variation.MoistureMin.HasValue)
                    moisture = rng.NextUniform(_variation.MoistureMin.Value, _variation.MoistureMax.Value);

                if (_variation.IgnitionOffsetRadius.HasValue && _variation.IgnitionOffsetRadius.Value > 0)
                {
                    double radius = _variation.IgnitionOffsetRadius.Value;
                    for (int p = 0; p < ignitionPoints.Count; p++)
                    {
                        var (ix, iy) = ignitionPoints[p];
                        int dx = (int)Math.Round(rng.NextUniform(-radius, radius));
                        int dy = (int)Math.Round(rng.NextUniform(-radius, radius));
                        int newIx = Math.Max(0, Math.Min(_grid.Nx - 1, ix + dx));
                        int newIy = Math.Max(0, Math.Min(_grid.Ny - 1, iy + dy));
                        ignitionPoints[p] = (newIx, newIy);
                    }
                }
            }

            // Apply moisture override for this iteration
            FuelMap iterFuelMap = _fuelMap;
            if (moisture >= 0)
            {
                // Create a shallow copy of moisture by overriding uniform
                // This is acceptable because we reset it each iteration
                _fuelMap.SetUniformMoisture(moisture);
            }

            var parameters = new WildfireParameters(
                ignitionPoints.AsReadOnly(),
                windSpeed,
                windDir,
                _burnDuration);

            var simulator = new WildfireSimulator(parameters);
            var snapshots = simulator.Run(_grid, _terrain, _fuelMap, _timeFrame);
            allSnapshots.Add(snapshots);

            // Accumulate burn statistics from last snapshot
            if (snapshots.Count > 0)
            {
                var last = snapshots[snapshots.Count - 1];
                burnedAreas[i] = last.BurnedAreaHectares;

                var bs = last.Snapshot.GetLayer("burnState");
                for (int c = 0; c < cellCount; c++)
                {
                    int state = (int)bs[c];
                    if (state == (int)CellBurnState.Burning || state == (int)CellBurnState.Burned)
                        burnCount[c]++;
                }
            }
        }

        // Compute burn probability
        var burnProb = new double[cellCount];
        for (int c = 0; c < cellCount; c++)
            burnProb[c] = (double)burnCount[c] / iterations;

        return new WildfireMonteCarloResult(iterations, _grid, burnProb, burnedAreas, allSnapshots);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Validation
    // ═══════════════════════════════════════════════════════════════

    private void Validate()
    {
        if (_grid == null) throw new InvalidOperationException("Grid not set. Call OverGrid().");
        if (_timeFrame == null) throw new InvalidOperationException("Time frame not set. Call OverTime().");
        if (_terrain == null) throw new InvalidOperationException("Terrain not set. Call WithTerrain().");
        if (_fuelMap == null) throw new InvalidOperationException("Fuel map not set. Call WithFuel().");
        if (_ignitionPoints.Count == 0) throw new InvalidOperationException("No ignition points. Call WithIgnition().");
    }
}
