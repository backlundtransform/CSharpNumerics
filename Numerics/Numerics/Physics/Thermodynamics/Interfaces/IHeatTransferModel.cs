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

    /// <summary>
    /// Computes the convective (Robin) boundary correction to the temperature rate
    /// at a boundary cell: −h / (ρ · c_p · dx) · (T_cell − T_ambient).
    /// This term is added to α·∇²T at faces with convective cooling.
    /// </summary>
    /// <param name="h">Convective heat transfer coefficient in W/(m²·K).</param>
    /// <param name="density">Material density ρ in kg/m³.</param>
    /// <param name="specificHeat">Specific heat capacity c_p in J/(kg·K).</param>
    /// <param name="dx">Grid spacing normal to the boundary face in metres.</param>
    /// <param name="cellTemperature">Temperature at the boundary cell in K.</param>
    /// <param name="ambientTemperature">Far-field fluid temperature T∞ in K.</param>
    double ConvectiveBoundaryRate(double h, double density, double specificHeat,
        double dx, double cellTemperature, double ambientTemperature);
}
