namespace CSharpNumerics.Physics.Materials.Engineering;

/// <summary>
/// Immutable material descriptor for multiphysics simulations.
/// Carries thermal, fluid, electrical, and structural properties.
/// </summary>
public readonly struct EngineeringMaterial
{
    /// <summary>Material name (e.g. "Steel").</summary>
    public string Name { get; }

    /// <summary>Thermal conductivity k in W/(m·K).</summary>
    public double ThermalConductivity { get; }

    /// <summary>Specific heat capacity c_p in J/(kg·K).</summary>
    public double SpecificHeat { get; }

    /// <summary>Density ρ in kg/m³.</summary>
    public double Density { get; }

    /// <summary>Dynamic viscosity μ in Pa·s.</summary>
    public double DynamicViscosity { get; }

    /// <summary>Relative electric permittivity ε_r (dimensionless).</summary>
    public double ElectricPermittivity { get; }

    /// <summary>Young's modulus E in Pa.</summary>
    public double YoungsModulus { get; }

    /// <summary>Poisson's ratio ν (dimensionless, typically 0–0.5).</summary>
    public double PoissonsRatio { get; }

    /// <summary>Thermal diffusivity α = k / (ρ · c_p) in m²/s.</summary>
    public double ThermalDiffusivity =>
        (Density > 0 && SpecificHeat > 0)
            ? ThermalConductivity / (Density * SpecificHeat)
            : 0.0;

    /// <summary>Kinematic viscosity ν = μ / ρ in m²/s.</summary>
    public double KinematicViscosity =>
        Density > 0 ? DynamicViscosity / Density : 0.0;

    public EngineeringMaterial(
        string name,
        double thermalConductivity,
        double specificHeat,
        double density,
        double dynamicViscosity,
        double electricPermittivity,
        double youngsModulus,
        double poissonsRatio)
    {
        Name = name;
        ThermalConductivity = thermalConductivity;
        SpecificHeat = specificHeat;
        Density = density;
        DynamicViscosity = dynamicViscosity;
        ElectricPermittivity = electricPermittivity;
        YoungsModulus = youngsModulus;
        PoissonsRatio = poissonsRatio;
    }
}
