namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Internal contract for a multiphysics solver.
/// Each simulation type has its own implementation.
/// </summary>
internal interface IMultiphysicsSolver
{
    SimulationResult Solve(SimulationBuilder config);
}
