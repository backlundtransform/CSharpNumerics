namespace CSharpNumerics.Physics.Thermodynamics.Interfaces;

/// <summary>
/// Abstraction for heat transfer physics used by thermal solvers.
/// Provides thermal diffusivity and volumetric heat source conversion.
/// </summary>
public interface IHeatTransferModel
{
    /// <summary>
    /// Computes thermal diffusivity α = k / (ρ · c_p).
    /// </summary>
    double ThermalDiffusivity(double conductivity, double density, double specificHeat);

    /// <summary>
    /// Converts a volumetric heat source power Q (W/m³) to the
    /// temperature rate contribution Q / (ρ · c_p), in K/s.
    /// </summary>
    double HeatSourceRate(double power, double density, double specificHeat);
}
