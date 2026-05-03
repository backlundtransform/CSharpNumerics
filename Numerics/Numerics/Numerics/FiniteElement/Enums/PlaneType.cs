namespace CSharpNumerics.Numerics.FiniteElement.Enums;

/// <summary>
/// Constitutive assumption for 2D elasticity problems.
/// </summary>
public enum PlaneType
{
    /// <summary>Thin bodies where out-of-plane stress is neglected.</summary>
    PlaneStress,

    /// <summary>Thick bodies where out-of-plane strain is constrained.</summary>
    PlaneStrain
}