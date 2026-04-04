using CSharpNumerics.Physics.Thermodynamics.Interfaces;

namespace CSharpNumerics.Physics.Thermodynamics;

/// <summary>
/// Default heat transfer model backed by <see cref="HeatExtensions"/>.
/// </summary>
public class HeatTransferModel : IHeatTransferModel
{
    public double ThermalDiffusivity(double conductivity, double density, double specificHeat)
    {
        return (density > 0 && specificHeat > 0)
            ? conductivity / (density * specificHeat)
            : 0.0;
    }

    public double HeatSourceRate(double power, double density, double specificHeat)
    {
        double rhoCp = density * specificHeat;
        return rhoCp > 0 ? power / rhoCp : 0;
    }
}
