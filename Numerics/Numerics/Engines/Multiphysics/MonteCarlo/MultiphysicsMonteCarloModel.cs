using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Engineering;
using CSharpNumerics.Statistics.MonteCarlo;
using CSharpNumerics.Statistics.Random;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.MonteCarlo;

/// <summary>
/// Monte Carlo adapter for multiphysics simulations. Samples parameters
/// from a <see cref="ParameterVariation"/> and runs the simulation per trial.
/// <para>
/// Provides two modes:
/// <list type="bullet">
///   <item><see cref="RunBatch"/> — runs N scenarios and returns the full
///   scenario matrix (rows = iterations, columns = output cells).</item>
///   <item><see cref="IMonteCarloModel.Evaluate"/> — returns a scalar summary
///   (peak output value) compatible with <see cref="Statistics.MonteCarlo.MonteCarloSimulator"/>.</item>
/// </list>
/// </para>
/// </summary>
public class MultiphysicsMonteCarloModel : IMonteCarloModel
{
    private readonly SimulationBuilder _baseline;
    private readonly ParameterVariation _variation;

    /// <summary>
    /// Creates a Monte Carlo model from a fully configured simulation builder
    /// and a set of parameter variation ranges.
    /// </summary>
    /// <param name="baseline">Baseline simulation config (used as template for each trial).</param>
    /// <param name="variation">Ranges for stochastic variation of parameters.</param>
    public MultiphysicsMonteCarloModel(SimulationBuilder baseline, ParameterVariation variation)
    {
        _baseline = baseline ?? throw new ArgumentNullException(nameof(baseline));
        _variation = variation ?? throw new ArgumentNullException(nameof(variation));
    }

    // ═══════════════════════════════════════════════════════════════
    //  IMonteCarloModel — scalar summary
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs one stochastic trial and returns the peak output value as a scalar.
    /// </summary>
    public double Evaluate(RandomGenerator rng)
    {
        var result = RunOneScenario(rng);
        return Math.Max(Math.Abs(result.MaxValue), Math.Abs(result.MinValue));
    }

    // ═══════════════════════════════════════════════════════════════
    //  Batch runner — full scenario matrix
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Runs N stochastic scenarios and returns a <see cref="MultiphysicsScenarioResult"/>
    /// containing the scenario matrix plus individual results.
    /// </summary>
    /// <param name="iterations">Number of Monte Carlo iterations.</param>
    /// <param name="seed">Optional seed for reproducibility.</param>
    public MultiphysicsScenarioResult RunBatch(int iterations, int? seed = null)
    {
        if (iterations <= 0) throw new ArgumentException("iterations must be positive.");

        var rng = seed.HasValue ? new RandomGenerator(seed.Value) : new RandomGenerator();

        // Run first to determine column count
        var firstResult = RunOneScenario(rng);
        int cols = GetFeatureCount(firstResult);

        var scenarioMatrix = new double[iterations, cols];
        var allResults = new List<SimulationResult>(iterations);

        // Fill first row
        FlattenResult(firstResult, scenarioMatrix, 0, cols);
        allResults.Add(firstResult);

        // Run remaining
        for (int i = 1; i < iterations; i++)
        {
            var result = RunOneScenario(rng);
            FlattenResult(result, scenarioMatrix, i, cols);
            allResults.Add(result);
        }

        return new MultiphysicsScenarioResult(
            new Matrix(scenarioMatrix),
            allResults,
            _baseline.Type);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Internal — single scenario with sampled parameters
    // ═══════════════════════════════════════════════════════════════

    private SimulationResult RunOneScenario(RandomGenerator rng)
    {
        var mat = _baseline.Material.Value;

        // Sample material properties
        double k = _variation.SampleOrDefault(rng,
            _variation.ThermalConductivityMin, _variation.ThermalConductivityMax,
            mat.ThermalConductivity);
        double E = _variation.SampleOrDefault(rng,
            _variation.YoungsModulusMin, _variation.YoungsModulusMax,
            mat.YoungsModulus);
        double mu = _variation.SampleOrDefault(rng,
            _variation.ViscosityMin, _variation.ViscosityMax,
            mat.DynamicViscosity);
        double rho = _variation.SampleOrDefault(rng,
            _variation.DensityMin, _variation.DensityMax,
            mat.Density);

        var sampledMaterial = new EngineeringMaterial(
            mat.Name, k, mat.SpecificHeat, rho,
            mu, mat.ElectricPermittivity, E, mat.PoissonsRatio);

        // Build a new simulation with sampled parameters
        var builder = SimulationType.Create(_baseline.Type)
            .WithMaterial(sampledMaterial);

        // Copy geometry
        if (_baseline.Nx > 0 && _baseline.Ny > 0)
            builder.WithGeometry(_baseline.GeomWidth, _baseline.GeomHeight, _baseline.Nx, _baseline.Ny);
        else if (_baseline.GeomRadius > 0)
            builder.WithGeometry(_baseline.GeomLength, _baseline.GeomRadius, _baseline.Nodes);
        else if (_baseline.Nodes > 0)
            builder.WithGeometry(_baseline.GeomLength, _baseline.Nodes);

        // Copy cross-section (beam)
        if (_baseline.SecondMomentOverride > 0)
            builder.WithSecondMoment(_baseline.SecondMomentOverride);
        else if (_baseline.SectionRadius > 0)
            builder.WithCrossSection(_baseline.SectionRadius);
        else if (_baseline.SectionWidth > 0 && _baseline.SectionHeight > 0)
            builder.WithCrossSection(_baseline.SectionWidth, _baseline.SectionHeight);

        // Copy/sample boundary conditions
        if (_baseline.Support.HasValue)
        {
            builder.WithBoundary(_baseline.Support.Value);
        }
        else if (_baseline.PressureGradient != 0)
        {
            double pg = _variation.SampleOrDefault(rng,
                _variation.PressureGradientMin, _variation.PressureGradientMax,
                _baseline.PressureGradient);
            builder.WithBoundary(pg);
        }
        else
        {
            double topBC = _variation.BoundaryValueMin.HasValue
                ? _variation.SampleOrDefault(rng, _variation.BoundaryValueMin, _variation.BoundaryValueMax, _baseline.TopBC)
                : _baseline.TopBC;
            builder.WithBoundary(topBC, _baseline.BottomBC, _baseline.LeftBC, _baseline.RightBC);
        }

        // Copy initial condition
        if (_baseline.HasInitialFunc)
            builder.WithInitialCondition(_baseline.InitialFunc);
        else
            builder.WithInitialCondition(_baseline.InitialValue);

        // Copy/sample sources
        foreach (var (ix, iy, value) in _baseline.Sources2D)
        {
            double v = _variation.SourceIntensityMin.HasValue
                ? _variation.SampleOrDefault(rng, _variation.SourceIntensityMin, _variation.SourceIntensityMax, value)
                : value;
            builder.AddSource(ix, iy, v);
        }

        foreach (var (pos, value) in _baseline.Sources1D)
        {
            double v = _variation.LoadMin.HasValue
                ? _variation.SampleOrDefault(rng, _variation.LoadMin, _variation.LoadMax, value)
                : value;
            builder.AddSource(pos, v);
        }

        if (_baseline.UniformLoad != 0)
        {
            double load = _variation.LoadMin.HasValue
                ? _variation.SampleOrDefault(rng, _variation.LoadMin, _variation.LoadMax, _baseline.UniformLoad)
                : _baseline.UniformLoad;
            builder.WithSource(load);
        }

        return builder.Solve(
            dt: _baseline.Dt,
            steps: _baseline.Steps,
            maxIterations: _baseline.MaxIterations,
            tolerance: _baseline.Tolerance);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers — flatten result into scenario matrix row
    // ═══════════════════════════════════════════════════════════════

    private static int GetFeatureCount(SimulationResult result)
    {
        if (result.Field != null)
            return result.Field.GetLength(0) * result.Field.GetLength(1);
        if (result.Values != null)
            return result.Values.Length;
        return 0;
    }

    private static void FlattenResult(SimulationResult result, double[,] matrix, int row, int cols)
    {
        if (result.Field != null)
        {
            int nx = result.Field.GetLength(0);
            int ny = result.Field.GetLength(1);
            int idx = 0;
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                    matrix[row, idx++] = result.Field[ix, iy];
        }
        else if (result.Values != null)
        {
            for (int i = 0; i < result.Values.Length && i < cols; i++)
                matrix[row, i] = result.Values[i];
        }
    }
}
