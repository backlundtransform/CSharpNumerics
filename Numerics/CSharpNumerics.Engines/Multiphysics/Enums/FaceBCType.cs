namespace CSharpNumerics.Engines.Multiphysics.Enums;

/// <summary>
/// Type of thermal boundary condition applied to a face.
/// </summary>
public enum FaceBCType
{
    /// <summary>Fixed temperature: T = T_wall.</summary>
    Dirichlet,

    /// <summary>Convective (Robin): −k ∂T/∂n = h(T − T∞).</summary>
    Convection
}
