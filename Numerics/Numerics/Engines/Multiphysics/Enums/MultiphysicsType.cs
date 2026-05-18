namespace CSharpNumerics.Engines.Multiphysics.Enums;

/// <summary>
/// The simulation types supported by the multiphysics engine.
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
    BeamStress,

    /// <summary>2D incompressible Navier–Stokes around a cylinder (Chorin projection method).</summary>
    CylinderFlow,

    /// <summary>2D transient Navier–Stokes on a rectangular domain (Chorin projection method).</summary>
    FluidFlow2D,

    /// <summary>2D magnetostatics: ∇²A = −μJ (vector potential formulation).</summary>
    MagneticField,

    /// <summary>2D plane-stress elasticity: displacement and stress tensor fields.</summary>
    PlaneStress,

    /// <summary>3D heat conduction: ∂T/∂t = α∇²T + source (volumetric).</summary>
    HeatBlock3D,

    /// <summary>3D advection-diffusion: ∂c/∂t = D∇²c − v·∇c + source.</summary>
    FluidDiffusion3D,

    /// <summary>3D incompressible Navier–Stokes around a cylinder (Chorin projection, periodic z).</summary>
    CylinderFlow3D,

    /// <summary>2D potential flow past an airfoil (Hess-Smith panel method).</summary>
    AirfoilFlow
}
