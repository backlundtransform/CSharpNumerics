namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Pre-defined optical media with standard refractive indices at 589 nm (sodium D-line).
/// </summary>
public static class OpticalMaterialLibrary
{
    // Gases
    public static OpticalMedium Vacuum => new("Vacuum", 1.0);
    public static OpticalMedium Air => new("Air", 1.000293);

    // Liquids
    public static OpticalMedium Water => new("Water", 1.333, abbeNumber: 55.8);
    public static OpticalMedium Ethanol => new("Ethanol", 1.361);
    public static OpticalMedium GlycerineOil => new("Glycerine", 1.473);

    // Common glasses
    public static OpticalMedium CrownGlass => new("Crown Glass (BK7)", 1.5168, abbeNumber: 64.17);
    public static OpticalMedium FlintGlass => new("Flint Glass (SF11)", 1.7847, abbeNumber: 25.76);
    public static OpticalMedium FusedSilica => new("Fused Silica", 1.4585, abbeNumber: 67.82);
    public static OpticalMedium BorosilicateGlass => new("Borosilicate Glass", 1.515, abbeNumber: 65.0);

    // Crystals
    public static OpticalMedium Diamond => new("Diamond", 2.417, abbeNumber: 55.3);
    public static OpticalMedium Sapphire => new("Sapphire", 1.770, abbeNumber: 72.2);
    public static OpticalMedium Quartz => new("Quartz", 1.544, abbeNumber: 69.8);
    public static OpticalMedium CalciteFluorite => new("Calcium Fluoride", 1.434, abbeNumber: 95.1);

    // Polymers
    public static OpticalMedium Acrylic => new("Acrylic (PMMA)", 1.491, abbeNumber: 57.4);
    public static OpticalMedium Polycarbonate => new("Polycarbonate", 1.585, abbeNumber: 30.0);

    // Specialty
    public static OpticalMedium Ice => new("Ice", 1.31);
    public static OpticalMedium SodiumChloride => new("Rock Salt (NaCl)", 1.544, abbeNumber: 42.9);
}
