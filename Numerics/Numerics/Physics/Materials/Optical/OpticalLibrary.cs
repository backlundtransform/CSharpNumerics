using CSharpNumerics.Physics.Optics;

namespace CSharpNumerics.Physics.Materials.Optical;

/// <summary>
/// Pre-defined optical materials for use in optics simulations.
/// Each material carries a refractive index, absorption coefficient and Abbe number.
/// </summary>
public static class OpticalLibrary
{
    // Gases
    public static OpticalMedium Vacuum => new("Vacuum", 1.0);
    public static OpticalMedium Air => new("Air", 1.000293);

    // Liquids
    public static OpticalMedium Water => new("Water", 1.333, abbeNumber: 55.8);
    public static OpticalMedium Ethanol => new("Ethanol", 1.361);

    // Common glasses
    public static OpticalMedium CrownGlass => new("Crown Glass (BK7)", 1.5168, abbeNumber: 64.17);
    public static OpticalMedium FlintGlass => new("Flint Glass (SF11)", 1.7847, abbeNumber: 25.76);
    public static OpticalMedium FusedSilica => new("Fused Silica", 1.4585, abbeNumber: 67.82);
    public static OpticalMedium BorosilicateGlass => new("Borosilicate Glass", 1.515, abbeNumber: 65.0);

    // Crystals
    public static OpticalMedium Diamond => new("Diamond", 2.417, abbeNumber: 55.3);
    public static OpticalMedium Sapphire => new("Sapphire", 1.770, abbeNumber: 72.2);
    public static OpticalMedium Quartz => new("Quartz", 1.544, abbeNumber: 69.8);
    public static OpticalMedium CalciumFluoride => new("Calcium Fluoride", 1.434, abbeNumber: 95.1);

    // Polymers
    public static OpticalMedium Acrylic => new("Acrylic (PMMA)", 1.491, abbeNumber: 57.4);
    public static OpticalMedium Polycarbonate => new("Polycarbonate", 1.585, abbeNumber: 30.0);

    // Specialty
    public static OpticalMedium Ice => new("Ice", 1.31);
    public static OpticalMedium SodiumChloride => new("Rock Salt (NaCl)", 1.544, abbeNumber: 42.9);

    /// <summary>
    /// Looks up an optical medium by name (case-insensitive).
    /// Supports common names: "vacuum", "air", "water", "crownglass", "diamond", etc.
    /// </summary>
    public static OpticalMedium Get(string name)
    {
        switch (name?.Trim().ToLowerInvariant().Replace(" ", ""))
        {
            case "vacuum": return Vacuum;
            case "air": return Air;
            case "water": return Water;
            case "ethanol": return Ethanol;
            case "crownglass":
            case "bk7": return CrownGlass;
            case "flintglass":
            case "sf11": return FlintGlass;
            case "fusedsilica": return FusedSilica;
            case "borosilicateglass": return BorosilicateGlass;
            case "diamond": return Diamond;
            case "sapphire": return Sapphire;
            case "quartz": return Quartz;
            case "calciumfluoride":
            case "caf2": return CalciumFluoride;
            case "acrylic":
            case "pmma": return Acrylic;
            case "polycarbonate": return Polycarbonate;
            case "ice": return Ice;
            case "sodiumchloride":
            case "nacl":
            case "rocksalt": return SodiumChloride;
            default:
                throw new System.ArgumentException($"Unknown optical material: '{name}'.", nameof(name));
        }
    }
}
