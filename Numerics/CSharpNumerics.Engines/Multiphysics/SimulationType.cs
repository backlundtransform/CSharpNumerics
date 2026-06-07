using CSharpNumerics.Engines.Multiphysics.Enums;

namespace CSharpNumerics.Engines.Multiphysics;

/// <summary>
/// Static entry point for the multiphysics fluent simulation pipeline.
/// <para>
/// Usage:
/// <code>
/// var result = SimulationType.Create(MultiphysicsType.HeatPlate)
///     .WithMaterial(EngineeringLibrary.Steel)
///     .WithGeometry(width: 0.1, height: 0.1, nx: 50, ny: 50)
///     .WithBoundary(top: 100, bottom: 0, left: 0, right: 0)
///     .WithInitialCondition(20.0)
///     .AddSource(25, 25, 1000.0)
///     .Solve(dt: 0.01, steps: 100);
///
/// var beam = SimulationType.Create(MultiphysicsType.BeamStress)
///     .WithMaterial(EngineeringLibrary.Steel)
///     .WithGeometry(length: 2.0, nodes: 100)
///     .WithCrossSection(width: 0.05, height: 0.1)
///     .WithBoundary(BeamSupport.Cantilever)
///     .AddSource(load: 1000.0, position: 2.0)
///     .Solve();
/// </code>
/// </para>
/// </summary>
public static class SimulationType
{
    /// <summary>
    /// Create a new multiphysics simulation builder for the given type.
    /// </summary>
    public static SimulationBuilder Create(MultiphysicsType type)
        => new SimulationBuilder(type);
}
