using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Materials.Water;
using CSharpNumerics.Statistics.Random;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Fluent builder for water contamination spread scenarios.
/// Follows the same pattern as <see cref="Wildfire.WildfireScenarioBuilder"/>.
/// Supports single deterministic runs and Monte Carlo ensembles.
/// <para>
/// Usage:
/// <code>
/// RiskScenario
///     .ForWaterContamination()
///     .WithRiverNetwork(riverNetwork)
///     .WithChannels(channelMap)
///     .WithTerrain(terrainGrid)
///     .WithSource(5, 0, 100, double.MaxValue)
///     .WithContaminant(ContaminantLibrary.Get("Benzene"))
///     .WithDischarge(10)
///     .OverGrid(grid)
///     .OverTime(0, 86400, 300)
///     .RunSingle();
/// </code>
/// </para>
/// </summary>
public class WaterContaminationScenarioBuilder
{
    // ── River topology ───────────────────────────────────────────
    private RiverNetwork _riverNetwork;
    private ChannelMap _channelMap;
    private TerrainGrid _terrain;

    // ── Contamination ────────────────────────────────────────────
    private readonly List<(int ix, int iy, double concentrationMgL, double durationSeconds)> _sources =
        new List<(int, int, double, double)>();
    private AquaticContaminant _contaminant;
    private double _baseDischarge = 10.0;
    private double _bedPorosity = 0.4;
    private double _bedBulkDensity = 1600.0;

    // ── Variation ────────────────────────────────────────────────
    private WaterContaminationVariation _variation;

    // ── Grid & time ──────────────────────────────────────────────
    private GeoGrid _grid;
    private TimeFrame _timeFrame;

    internal WaterContaminationScenarioBuilder() { }

    // ═══════════════════════════════════════════════════════════════
    //  Configuration steps (return self)
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Sets the river network flow graph.</summary>
    public WaterContaminationScenarioBuilder WithRiverNetwork(RiverNetwork network)
    {
        _riverNetwork = network ?? throw new ArgumentNullException(nameof(network));
        return this;
    }

    /// <summary>Sets the channel hydraulic properties.</summary>
    public WaterContaminationScenarioBuilder WithChannels(ChannelMap channels)
    {
        _channelMap = channels ?? throw new ArgumentNullException(nameof(channels));
        return this;
    }

    /// <summary>Sets the terrain elevation grid.</summary>
    public WaterContaminationScenarioBuilder WithTerrain(TerrainGrid terrain)
    {
        _terrain = terrain ?? throw new ArgumentNullException(nameof(terrain));
        return this;
    }

    /// <summary>Adds a contamination point-source.</summary>
    public WaterContaminationScenarioBuilder WithSource(
        int ix, int iy, double concentrationMgL, double durationSeconds)
    {
        _sources.Add((ix, iy, concentrationMgL, durationSeconds));
        return this;
    }

    /// <summary>Sets the contaminant descriptor.</summary>
    public WaterContaminationScenarioBuilder WithContaminant(AquaticContaminant contaminant)
    {
        _contaminant = contaminant;
        return this;
    }

    /// <summary>Sets the baseline river discharge (m³/s).</summary>
    public WaterContaminationScenarioBuilder WithDischarge(double dischargeM3s)
    {
        _baseDischarge = dischargeM3s;
        return this;
    }

    /// <summary>Sets bed sediment properties for retardation calculation.</summary>
    public WaterContaminationScenarioBuilder WithBedProperties(double porosity, double bulkDensityKgM3)
    {
        _bedPorosity = porosity;
        _bedBulkDensity = bulkDensityKgM3;
        return this;
    }

    /// <summary>Configures stochastic variation for Monte Carlo runs.</summary>
    public WaterContaminationScenarioBuilder WithVariation(Action<WaterContaminationVariation> configure)
    {
        _variation = new WaterContaminationVariation();
        configure(_variation);
        return this;
    }

    /// <summary>Sets the spatial grid.</summary>
    public WaterContaminationScenarioBuilder OverGrid(GeoGrid grid)
    {
        _grid = grid ?? throw new ArgumentNullException(nameof(grid));
        return this;
    }

    /// <summary>Sets the simulation time frame.</summary>
    public WaterContaminationScenarioBuilder OverTime(double startSeconds, double endSeconds, double stepSeconds)
    {
        _timeFrame = new TimeFrame(startSeconds, endSeconds, stepSeconds);
        return this;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Terminal: Run single deterministic scenario
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs a single deterministic water contamination simulation.
    /// </summary>
    public WaterContaminationResult RunSingle()
    {
        Validate();

        var parameters = new WaterContaminationParameters(
            _sources.AsReadOnly(), _contaminant,
            _baseDischarge, _bedPorosity, _bedBulkDensity);

        var simulator = new WaterContaminationSimulator(parameters, _riverNetwork, _channelMap);
        var snapshots = simulator.Run(_grid, _terrain, new FuelMap(_grid), _timeFrame);

        return new WaterContaminationResult(snapshots, _grid, _contaminant.ToxicityThresholdMgL);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Terminal: Run Monte Carlo ensemble
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs a Monte Carlo ensemble of water contamination simulations with
    /// stochastic parameter variation. Returns per-cell exceedance probability
    /// and peak concentration statistics.
    /// </summary>
    public WaterContaminationMonteCarloResult RunMonteCarlo(int iterations, int? seed = null)
    {
        Validate();
        if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

        var rng = seed.HasValue ? new RandomGenerator(seed.Value) : new RandomGenerator();
        int cellCount = _grid.Nx * _grid.Ny;
        double toxThreshold = _contaminant.ToxicityThresholdMgL;

        var exceedCount = new int[cellCount];
        var peakConcentrations = new double[iterations];
        double earliestArrival = double.MaxValue;
        var allSnapshots = new List<IReadOnlyList<SpreadSnapshot>>(iterations);

        for (int i = 0; i < iterations; i++)
        {
            // Sample parameters
            double discharge = _baseDischarge;
            double concMultiplier = 1.0;
            double manningN = -1; // sentinel: no override

            if (_variation != null && _variation.HasVariation)
            {
                if (_variation.DischargeMin.HasValue)
                    discharge = rng.NextUniform(_variation.DischargeMin.Value, _variation.DischargeMax.Value);

                if (_variation.SourceConcentrationMin.HasValue)
                    concMultiplier = rng.NextUniform(
                        _variation.SourceConcentrationMin.Value,
                        _variation.SourceConcentrationMax.Value) / _sources[0].concentrationMgL;

                if (_variation.ManningNMin.HasValue)
                    manningN = rng.NextUniform(_variation.ManningNMin.Value, _variation.ManningNMax.Value);
            }

            // Build iteration-specific sources with concentration multiplier
            var iterSources = new List<(int, int, double, double)>(_sources.Count);
            foreach (var (ix, iy, conc, dur) in _sources)
                iterSources.Add((ix, iy, conc * concMultiplier, dur));

            var parameters = new WaterContaminationParameters(
                iterSources.AsReadOnly(), _contaminant,
                discharge, _bedPorosity, _bedBulkDensity);

            // Apply Manning's n override if varied
            ChannelMap iterChannels = _channelMap;
            if (manningN >= 0)
            {
                iterChannels = new ChannelMap(_grid);
                iterChannels.SetUniformChannel(
                    _channelMap.GetWidth(0, 0),
                    _channelMap.GetDepth(0, 0),
                    manningN);
            }

            var simulator = new WaterContaminationSimulator(parameters, _riverNetwork, iterChannels);
            var snapshots = simulator.Run(_grid, _terrain, new FuelMap(_grid), _timeFrame);
            allSnapshots.Add(snapshots);

            // Accumulate statistics
            double iterPeak = 0;
            bool exceeded = false;
            foreach (var snap in snapshots)
            {
                double m = snap.MaxConcentration;
                if (m > iterPeak) iterPeak = m;

                if (!exceeded && m > toxThreshold)
                {
                    exceeded = true;
                    if (snap.Time < earliestArrival)
                        earliestArrival = snap.Time;
                }

                var conc = snap.Snapshot.GetLayer("concentration");
                for (int c = 0; c < cellCount; c++)
                {
                    if (conc[c] > toxThreshold)
                        exceedCount[c]++;
                }
            }

            peakConcentrations[i] = iterPeak;
        }

        // Convert exceedance counts to probability (any time step in any iteration)
        // Re-compute as binary per-iteration exceedance
        var binaryExceed = new int[cellCount];
        for (int i = 0; i < iterations; i++)
        {
            var peaked = new bool[cellCount];
            foreach (var snap in allSnapshots[i])
            {
                var conc = snap.Snapshot.GetLayer("concentration");
                for (int c = 0; c < cellCount; c++)
                    if (conc[c] > toxThreshold) peaked[c] = true;
            }
            for (int c = 0; c < cellCount; c++)
                if (peaked[c]) binaryExceed[c]++;
        }

        var exceedProb = new double[cellCount];
        for (int c = 0; c < cellCount; c++)
            exceedProb[c] = (double)binaryExceed[c] / iterations;

        if (earliestArrival == double.MaxValue)
            earliestArrival = 0;

        return new WaterContaminationMonteCarloResult(
            iterations, _grid, exceedProb, peakConcentrations,
            earliestArrival, allSnapshots, _timeFrame);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Validation
    // ═══════════════════════════════════════════════════════════════

    private void Validate()
    {
        if (_grid == null) throw new InvalidOperationException("Grid not set. Call OverGrid().");
        if (_timeFrame == null) throw new InvalidOperationException("Time frame not set. Call OverTime().");
        if (_terrain == null) throw new InvalidOperationException("Terrain not set. Call WithTerrain().");
        if (_riverNetwork == null) throw new InvalidOperationException("River network not set. Call WithRiverNetwork().");
        if (_channelMap == null) throw new InvalidOperationException("Channel map not set. Call WithChannels().");
        if (_sources.Count == 0) throw new InvalidOperationException("No sources. Call WithSource().");
        if (_contaminant.Name == null) throw new InvalidOperationException("Contaminant not set. Call WithContaminant().");
    }
}
