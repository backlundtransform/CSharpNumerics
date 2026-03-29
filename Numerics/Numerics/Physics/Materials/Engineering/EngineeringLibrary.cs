namespace CSharpNumerics.Physics.Materials.Engineering;

/// <summary>
/// Pre-defined engineering materials with typical room-temperature properties.
/// </summary>
public static class EngineeringLibrary
{
    // Metals
    public static EngineeringMaterial Steel => new(
        "Steel", thermalConductivity: 50.0, specificHeat: 500.0, density: 7850.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 200e9, poissonsRatio: 0.30);

    public static EngineeringMaterial Aluminum => new(
        "Aluminum", thermalConductivity: 237.0, specificHeat: 897.0, density: 2700.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 69e9, poissonsRatio: 0.33);

    public static EngineeringMaterial Copper => new(
        "Copper", thermalConductivity: 401.0, specificHeat: 385.0, density: 8960.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 120e9, poissonsRatio: 0.34);

    // Fluids
    public static EngineeringMaterial Water => new(
        "Water", thermalConductivity: 0.606, specificHeat: 4186.0, density: 998.0,
        dynamicViscosity: 1.002e-3, electricPermittivity: 80.1, youngsModulus: 0, poissonsRatio: 0.5);

    public static EngineeringMaterial Air => new(
        "Air", thermalConductivity: 0.026, specificHeat: 1005.0, density: 1.225,
        dynamicViscosity: 1.81e-5, electricPermittivity: 1.0006, youngsModulus: 0, poissonsRatio: 0);

    // Construction
    public static EngineeringMaterial Concrete => new(
        "Concrete", thermalConductivity: 1.7, specificHeat: 880.0, density: 2400.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 30e9, poissonsRatio: 0.20);

    public static EngineeringMaterial Glass => new(
        "Glass", thermalConductivity: 1.05, specificHeat: 840.0, density: 2500.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 70e9, poissonsRatio: 0.22);
}
