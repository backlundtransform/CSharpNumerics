namespace CSharpNumerics.Physics.Materials.Engineering;

/// <summary>
/// Pre-defined engineering materials with typical room-temperature properties.
/// </summary>
public static class EngineeringLibrary
{
    // ── Metals ───────────────────────────────────────────────────

    public static EngineeringMaterial Steel => new(
        "Steel", thermalConductivity: 50.0, specificHeat: 500.0, density: 7850.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 200e9, poissonsRatio: 0.30,
        magneticPermeability: 100.0);

    public static EngineeringMaterial Aluminum => new(
        "Aluminum", thermalConductivity: 237.0, specificHeat: 897.0, density: 2700.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 69e9, poissonsRatio: 0.33,
        magneticPermeability: 1.000022);

    public static EngineeringMaterial Copper => new(
        "Copper", thermalConductivity: 401.0, specificHeat: 385.0, density: 8960.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 120e9, poissonsRatio: 0.34,
        magneticPermeability: 0.999994);

    public static EngineeringMaterial Titanium => new(
        "Titanium", thermalConductivity: 21.9, specificHeat: 523.0, density: 4507.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 116e9, poissonsRatio: 0.34,
        magneticPermeability: 1.00018);

    public static EngineeringMaterial Brass => new(
        "Brass", thermalConductivity: 109.0, specificHeat: 380.0, density: 8500.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 100e9, poissonsRatio: 0.34,
        magneticPermeability: 1.0);

    public static EngineeringMaterial StainlessSteel => new(
        "Stainless Steel", thermalConductivity: 16.0, specificHeat: 500.0, density: 8000.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 193e9, poissonsRatio: 0.29,
        magneticPermeability: 1.02);

    // ── Fluids ───────────────────────────────────────────────────

    public static EngineeringMaterial Water => new(
        "Water", thermalConductivity: 0.606, specificHeat: 4186.0, density: 998.0,
        dynamicViscosity: 1.002e-3, electricPermittivity: 80.1, youngsModulus: 0, poissonsRatio: 0.5,
        magneticPermeability: 0.999992);

    public static EngineeringMaterial Air => new(
        "Air", thermalConductivity: 0.026, specificHeat: 1005.0, density: 1.225,
        dynamicViscosity: 1.81e-5, electricPermittivity: 1.0006, youngsModulus: 0, poissonsRatio: 0,
        magneticPermeability: 1.00000037);

    public static EngineeringMaterial Oil => new(
        "Oil", thermalConductivity: 0.145, specificHeat: 2000.0, density: 870.0,
        dynamicViscosity: 0.03, electricPermittivity: 2.2, youngsModulus: 0, poissonsRatio: 0.5,
        magneticPermeability: 1.0);

    public static EngineeringMaterial Glycerin => new(
        "Glycerin", thermalConductivity: 0.285, specificHeat: 2430.0, density: 1261.0,
        dynamicViscosity: 1.412, electricPermittivity: 42.5, youngsModulus: 0, poissonsRatio: 0.5,
        magneticPermeability: 1.0);

    // ── Construction / Structural ────────────────────────────────

    public static EngineeringMaterial Concrete => new(
        "Concrete", thermalConductivity: 1.7, specificHeat: 880.0, density: 2400.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 30e9, poissonsRatio: 0.20,
        magneticPermeability: 1.0);

    public static EngineeringMaterial Glass => new(
        "Glass", thermalConductivity: 1.05, specificHeat: 840.0, density: 2500.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 70e9, poissonsRatio: 0.22,
        magneticPermeability: 1.0);

    public static EngineeringMaterial Wood => new(
        "Wood", thermalConductivity: 0.15, specificHeat: 1700.0, density: 600.0,
        dynamicViscosity: 0, electricPermittivity: 0, youngsModulus: 12e9, poissonsRatio: 0.35,
        magneticPermeability: 1.0);

    public static EngineeringMaterial Rubber => new(
        "Rubber", thermalConductivity: 0.16, specificHeat: 2010.0, density: 1100.0,
        dynamicViscosity: 0, electricPermittivity: 7.0, youngsModulus: 0.01e9, poissonsRatio: 0.49,
        magneticPermeability: 1.0);

    public static EngineeringMaterial Plastic => new(
        "Plastic (HDPE)", thermalConductivity: 0.50, specificHeat: 1900.0, density: 960.0,
        dynamicViscosity: 0, electricPermittivity: 2.3, youngsModulus: 1.1e9, poissonsRatio: 0.42,
        magneticPermeability: 1.0);
}
