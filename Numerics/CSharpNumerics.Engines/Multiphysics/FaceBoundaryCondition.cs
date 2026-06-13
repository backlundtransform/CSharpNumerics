using CSharpNumerics.Engines.Multiphysics.Enums;

namespace CSharpNumerics.Engines.Multiphysics;

/// <summary>
/// Describes the thermal boundary condition on a single face.
/// Either Dirichlet (fixed temperature) or Robin/convection (−k ∂T/∂n = h(T − T∞)).
/// </summary>
public sealed class FaceBoundaryCondition
{
    /// <summary>Type of boundary condition.</summary>
    public FaceBCType Type { get; }

    /// <summary>Fixed temperature for Dirichlet, or ambient temperature T∞ for convection (K).</summary>
    public double Temperature { get; }

    /// <summary>Convective heat transfer coefficient h in W/(m²·K). Only used for <see cref="FaceBCType.Convection"/>.</summary>
    public double HeatTransferCoefficient { get; }

    private FaceBoundaryCondition(FaceBCType type, double temperature, double h)
    {
        Type = type;
        Temperature = temperature;
        HeatTransferCoefficient = h;
    }

    /// <summary>Create a Dirichlet (fixed temperature) boundary condition.</summary>
    public static FaceBoundaryCondition Dirichlet(double temperature)
        => new FaceBoundaryCondition(FaceBCType.Dirichlet, temperature, 0);

    /// <summary>
    /// Create a Robin/convection boundary condition: −k ∂T/∂n = h(T − T∞).
    /// </summary>
    /// <param name="heatTransferCoefficient">Convective coefficient h in W/(m²·K).</param>
    /// <param name="ambientTemperature">Far-field fluid temperature T∞ in K.</param>
    public static FaceBoundaryCondition Convection(double heatTransferCoefficient, double ambientTemperature)
        => new FaceBoundaryCondition(FaceBCType.Convection, ambientTemperature, heatTransferCoefficient);
}
