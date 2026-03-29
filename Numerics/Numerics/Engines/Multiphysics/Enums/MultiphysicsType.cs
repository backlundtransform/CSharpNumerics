namespace CSharpNumerics.Engines.Multiphysics.Enums;

/// <summary>
/// The four simulation types supported by the multiphysics engine.
/// </summary>
public enum MultiphysicsType
{
    /// <summary>2D heat conduction: ∂T/∂t = α∇²T + source.</summary>
    HeatPlate,

    /// <summary>1D transient pipe flow (Hagen–Poiseuille): ∂v/∂t = ν∇²v − (1/ρ)dP/dx.</summary>
    PipeFlow,

    /// <summary>2D electrostatics: ∇²φ = −ρ/ε (Poisson/Laplace).</summary>
    ElectricField,

    /// <summary>1D static beam analysis (Euler–Bernoulli): EIu⁗ = q.</summary>
    BeamStress
}
