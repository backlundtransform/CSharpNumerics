namespace CSharpNumerics.Physics.Electromagnetism.Interfaces;

/// <summary>
/// Abstraction for electrostatic physics used by electric field solvers.
/// Provides vacuum permittivity and charge-density RHS computation.
/// </summary>
public interface IElectrostaticModel
{
    /// <summary>
    /// Vacuum permittivity ε₀ in F/m.
    /// </summary>
    double VacuumPermittivity { get; }

    /// <summary>
    /// Computes effective permittivity ε = ε₀ · ε_r.
    /// </summary>
    double EffectivePermittivity(double relativePermittivity);

    /// <summary>
    /// Computes the Poisson RHS term −ρ / ε for a given charge density.
    /// </summary>
    double PoissonRhs(double chargeDensity, double permittivity);
}
