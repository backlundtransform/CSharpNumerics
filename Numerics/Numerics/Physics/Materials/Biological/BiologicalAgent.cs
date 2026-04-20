using System;
using System.Globalization;

namespace CSharpNumerics.Physics.Materials.Biological;

/// <summary>
/// Broad biological aerosol classes supported by the dispersion pipeline.
/// </summary>
public enum BiologicalAgentClass
{
    Virus,
    Bacteria,
    Spore
}

/// <summary>
/// Immutable descriptor of a biological material that can be attached to a
/// dispersion scenario. The descriptor stores enough information to derive
/// screening-level biological layers from a transported mass concentration.
/// </summary>
public readonly struct BiologicalAgent : IEquatable<BiologicalAgent>
{
    /// <summary>Short identifier used for lookup, for example "Virus".</summary>
    public string Code { get; }

    /// <summary>Human-readable name.</summary>
    public string Name { get; }

    /// <summary>High-level aerosol class.</summary>
    public BiologicalAgentClass Classification { get; }

    /// <summary>Representative aerodynamic diameter in microns.</summary>
    public double TypicalDiameterMicrons { get; }

    /// <summary>
    /// Representative mass of one biological unit in kilograms.
    /// Used to convert kg/m^3 into units/m^3.
    /// </summary>
    public double UnitMassKg { get; }

    /// <summary>
    /// Screening half-life for viability in air, in seconds.
    /// </summary>
    public double ViabilityHalfLifeSeconds { get; }

    /// <summary>True when viability is treated as non-decaying.</summary>
    public bool IsPersistent => double.IsPositiveInfinity(ViabilityHalfLifeSeconds);

    /// <summary>
    /// Creates a biological agent descriptor.
    /// </summary>
    public BiologicalAgent(
        string code,
        string name,
        BiologicalAgentClass classification,
        double typicalDiameterMicrons,
        double unitMassKg,
        double viabilityHalfLifeSeconds)
    {
        if (code == null) throw new ArgumentNullException(nameof(code));
        if (name == null) throw new ArgumentNullException(nameof(name));
        if (typicalDiameterMicrons <= 0) throw new ArgumentOutOfRangeException(nameof(typicalDiameterMicrons));
        if (unitMassKg <= 0) throw new ArgumentOutOfRangeException(nameof(unitMassKg));
        if (viabilityHalfLifeSeconds <= 0 && !double.IsPositiveInfinity(viabilityHalfLifeSeconds))
            throw new ArgumentOutOfRangeException(nameof(viabilityHalfLifeSeconds));

        Code = code;
        Name = name;
        Classification = classification;
        TypicalDiameterMicrons = typicalDiameterMicrons;
        UnitMassKg = unitMassKg;
        ViabilityHalfLifeSeconds = viabilityHalfLifeSeconds;
    }

    /// <summary>
    /// Converts a mass concentration in kg/m^3 into biological units per m^3.
    /// </summary>
    public double KgM3ToUnitsPerM3(double kgPerM3)
    {
        return kgPerM3 / UnitMassKg;
    }

    /// <summary>
    /// Converts biological units per m^3 back to kg/m^3.
    /// </summary>
    public double UnitsPerM3ToKgM3(double unitsPerM3)
    {
        return unitsPerM3 * UnitMassKg;
    }

    /// <summary>
    /// Returns the surviving viable fraction at the given elapsed time.
    /// </summary>
    public double ViabilityFractionAt(double timeSeconds)
    {
        if (timeSeconds <= 0 || IsPersistent)
            return 1.0;

        return Math.Exp(-Math.Log(2.0) * timeSeconds / ViabilityHalfLifeSeconds);
    }

    /// <summary>Generic viral aerosol for screening scenarios.</summary>
    public static readonly BiologicalAgent GenericVirus = new BiologicalAgent(
        code: "Virus",
        name: "Generic Viral Aerosol",
        classification: BiologicalAgentClass.Virus,
        typicalDiameterMicrons: 0.12,
        unitMassKg: 1e-18,
        viabilityHalfLifeSeconds: 6 * 3600);

    /// <summary>Generic bacterial aerosol for screening scenarios.</summary>
    public static readonly BiologicalAgent GenericBacteria = new BiologicalAgent(
        code: "Bacteria",
        name: "Generic Bacterial Aerosol",
        classification: BiologicalAgentClass.Bacteria,
        typicalDiameterMicrons: 1.5,
        unitMassKg: 1e-15,
        viabilityHalfLifeSeconds: 12 * 3600);

    /// <summary>Generic biological spore aerosol for screening scenarios.</summary>
    public static readonly BiologicalAgent GenericSpore = new BiologicalAgent(
        code: "Spore",
        name: "Generic Biological Spore",
        classification: BiologicalAgentClass.Spore,
        typicalDiameterMicrons: 2.5,
        unitMassKg: 3e-15,
        viabilityHalfLifeSeconds: 7 * 24 * 3600);

    public bool Equals(BiologicalAgent other) =>
        string.Equals(Code, other.Code, StringComparison.OrdinalIgnoreCase);

    public override bool Equals(object obj) =>
        obj is BiologicalAgent other && Equals(other);

    public override int GetHashCode() =>
        StringComparer.OrdinalIgnoreCase.GetHashCode(Code ?? string.Empty);

    public static bool operator ==(BiologicalAgent left, BiologicalAgent right) => left.Equals(right);
    public static bool operator !=(BiologicalAgent left, BiologicalAgent right) => !left.Equals(right);

    public override string ToString()
    {
        string halfLife = IsPersistent
            ? "persistent"
            : string.Format(CultureInfo.InvariantCulture, "t1/2={0:F1} h", ViabilityHalfLifeSeconds / 3600.0);

        return string.Format(
            CultureInfo.InvariantCulture,
            "{0} ({1}), class={2}, {3}",
            Name,
            Code,
            Classification,
            halfLife);
    }
}
