using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Electromagnetism.Interfaces;

namespace CSharpNumerics.Physics.Electromagnetism;

/// <summary>
/// Default electrostatic model using constants from <see cref="PhysicsConstants"/>.
/// </summary>
public class ElectrostaticModel : IElectrostaticModel
{
    public double VacuumPermittivity => PhysicsConstants.VacuumPermittivity;

    public double EffectivePermittivity(double relativePermittivity)
    {
        double epsR = relativePermittivity > 0 ? relativePermittivity : 1.0;
        return PhysicsConstants.VacuumPermittivity * epsR;
    }

    public double PoissonRhs(double chargeDensity, double permittivity)
    {
        return -chargeDensity / permittivity;
    }
}
