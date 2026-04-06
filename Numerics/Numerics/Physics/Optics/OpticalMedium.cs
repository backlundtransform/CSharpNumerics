namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Describes the optical properties of a medium through which light travels.
/// </summary>
public readonly struct OpticalMedium
{
    /// <summary>Material name.</summary>
    public string Name { get; }

    /// <summary>
    /// Index of refraction n (dimensionless). Vacuum = 1.0.
    /// May be a function of wavelength; this value is the reference at 589 nm (sodium D-line).
    /// </summary>
    public double RefractiveIndex { get; }

    /// <summary>
    /// Absorption coefficient α in 1/m (Beer–Lambert law).
    /// Zero for perfectly transparent media.
    /// </summary>
    public double AbsorptionCoefficient { get; }

    /// <summary>
    /// Abbe number V_d characterising chromatic dispersion.
    /// Higher values mean lower dispersion. Zero if not specified.
    /// </summary>
    public double AbbeNumber { get; }

    public OpticalMedium(string name, double refractiveIndex,
        double absorptionCoefficient = 0.0, double abbeNumber = 0.0)
    {
        Name = name;
        RefractiveIndex = refractiveIndex;
        AbsorptionCoefficient = absorptionCoefficient;
        AbbeNumber = abbeNumber;
    }

    /// <summary>
    /// Returns an approximate refractive index for a given wavelength using the Cauchy
    /// dispersion formula: n(λ) ≈ n_d + (n_d − 1) / V_d · (C₁/λ² − C₂) where
    /// C₁ and C₂ are calibrated to the sodium D-line.
    /// Falls back to <see cref="RefractiveIndex"/> when Abbe number is zero.
    /// </summary>
    public double RefractiveIndexAt(double wavelengthNm)
    {
        if (AbbeNumber <= 0) return RefractiveIndex;

        // Cauchy-like approximation using Abbe number
        const double lambdaD = 589.3; // nm (sodium D-line)
        const double lambdaF = 486.1; // nm (hydrogen F-line)
        const double lambdaC = 656.3; // nm (hydrogen C-line)

        double nF_nC = (RefractiveIndex - 1.0) / AbbeNumber;
        double invLF2 = 1.0 / (lambdaF * lambdaF);
        double invLC2 = 1.0 / (lambdaC * lambdaC);
        double B = nF_nC / (invLF2 - invLC2);
        double A = RefractiveIndex - B / (lambdaD * lambdaD);

        return A + B / (wavelengthNm * wavelengthNm);
    }

    public override string ToString() =>
        $"{Name} (n={RefractiveIndex:F4}, α={AbsorptionCoefficient:G3}, V={AbbeNumber:F1})";
}
