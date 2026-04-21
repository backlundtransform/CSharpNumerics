using System;
using System.Globalization;
using CSharpNumerics.Physics.Materials.Water.Enums;

namespace CSharpNumerics.Physics.Materials.Water;

/// <summary>
/// Immutable descriptor of an aquatic contaminant carrying the physical,
/// chemical, and toxicological properties needed for waterway transport modelling.
/// <para>
/// Common contaminants are available as static fields:
/// <c>AquaticContaminant.Benzene</c>, <c>AquaticContaminant.Cs137</c>, etc.
/// </para>
/// </summary>
public readonly struct AquaticContaminant : IEquatable<AquaticContaminant>
{
    /// <summary>Short identifier (e.g. "Cs137", "Benzene", "EColi").</summary>
    public string Name { get; }

    /// <summary>Primary hazard classification.</summary>
    public ContaminantType Type { get; }

    /// <summary>
    /// Half-life in seconds. Zero or <see cref="double.PositiveInfinity"/>
    /// indicates a conservative (non-decaying) tracer.
    /// </summary>
    public double HalfLifeSeconds { get; }

    /// <summary>
    /// Solid–water partition coefficient Kd (L/kg).
    /// Controls sorption to bed sediment. Zero means no adsorption.
    /// </summary>
    public double PartitionCoefficient { get; }

    /// <summary>
    /// Toxicity threshold concentration (mg/L) — e.g. drinking water limit.
    /// Cells above this are considered contaminated.
    /// </summary>
    public double ToxicityThresholdMgL { get; }

    /// <summary>
    /// Lethal threshold concentration (mg/L) — acute lethality level.
    /// </summary>
    public double LethalThresholdMgL { get; }

    /// <summary>
    /// Creates an aquatic contaminant descriptor.
    /// </summary>
    public AquaticContaminant(
        string name,
        ContaminantType type,
        double halfLifeSeconds,
        double partitionCoefficient,
        double toxicityThresholdMgL,
        double lethalThresholdMgL)
    {
        Name = name ?? throw new ArgumentNullException(nameof(name));
        Type = type;
        HalfLifeSeconds = halfLifeSeconds;
        PartitionCoefficient = partitionCoefficient;
        ToxicityThresholdMgL = toxicityThresholdMgL;
        LethalThresholdMgL = lethalThresholdMgL;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Derived helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// First-order decay constant λ = ln(2) / t½ (1/s).
    /// Returns 0 for conservative tracers (half-life ≤ 0 or infinite).
    /// </summary>
    public double DecayConstant =>
        HalfLifeSeconds > 0 && !double.IsPositiveInfinity(HalfLifeSeconds)
            ? Math.Log(2) / HalfLifeSeconds
            : 0.0;

    /// <summary>
    /// Returns true if this contaminant has no decay (conservative tracer).
    /// </summary>
    public bool IsConservative =>
        HalfLifeSeconds <= 0 || double.IsPositiveInfinity(HalfLifeSeconds);

    // ═══════════════════════════════════════════════════════════════
    //  Built-in contaminants
    // ═══════════════════════════════════════════════════════════════

    // ── Radioactive ──────────────────────────────────────────────

    /// <summary>Caesium-137. t½ ≈ 30.17 yr. Strong sediment sorption (Kd ≈ 1000 L/kg).</summary>
    public static readonly AquaticContaminant Cs137 = new AquaticContaminant(
        "Cs137", ContaminantType.Radioactive,
        halfLifeSeconds: 30.17 * 365.25 * 86400,   // ≈ 9.52 × 10⁸ s
        partitionCoefficient: 1000,                  // L/kg — clay/silt affinity
        toxicityThresholdMgL: 0.0001,               // very low (radiation hazard)
        lethalThresholdMgL: 0.01);

    /// <summary>Strontium-90. t½ ≈ 28.8 yr. Moderate sediment sorption.</summary>
    public static readonly AquaticContaminant Sr90 = new AquaticContaminant(
        "Sr90", ContaminantType.Radioactive,
        halfLifeSeconds: 28.8 * 365.25 * 86400,
        partitionCoefficient: 100,
        toxicityThresholdMgL: 0.0001,
        lethalThresholdMgL: 0.01);

    /// <summary>Iodine-131. t½ ≈ 8.02 days. Low sorption, highly mobile.</summary>
    public static readonly AquaticContaminant I131 = new AquaticContaminant(
        "I131", ContaminantType.Radioactive,
        halfLifeSeconds: 8.02 * 86400,
        partitionCoefficient: 5,
        toxicityThresholdMgL: 0.0001,
        lethalThresholdMgL: 0.001);

    // ── Chemical ─────────────────────────────────────────────────

    /// <summary>Benzene. WHO drinking water guideline: 0.01 mg/L.</summary>
    public static readonly AquaticContaminant Benzene = new AquaticContaminant(
        "Benzene", ContaminantType.Chemical,
        halfLifeSeconds: 10 * 86400,                // ≈ 10 days biodegradation in water
        partitionCoefficient: 0.6,                   // low sorption
        toxicityThresholdMgL: 0.01,                  // WHO guideline
        lethalThresholdMgL: 500);                    // acute oral LD50-equivalent

    /// <summary>Toluene. WHO drinking water guideline: 0.7 mg/L.</summary>
    public static readonly AquaticContaminant Toluene = new AquaticContaminant(
        "Toluene", ContaminantType.Chemical,
        halfLifeSeconds: 7 * 86400,
        partitionCoefficient: 1.0,
        toxicityThresholdMgL: 0.7,
        lethalThresholdMgL: 600);

    /// <summary>Cyanide (free CN⁻). WHO guideline: 0.07 mg/L.</summary>
    public static readonly AquaticContaminant Cyanide = new AquaticContaminant(
        "Cyanide", ContaminantType.Chemical,
        halfLifeSeconds: 0,                          // conservative — photolysis/volatilization neglected
        partitionCoefficient: 0,
        toxicityThresholdMgL: 0.07,
        lethalThresholdMgL: 3.0);

    /// <summary>Mercury (inorganic Hg). WHO guideline: 0.006 mg/L. High Kd.</summary>
    public static readonly AquaticContaminant Mercury = new AquaticContaminant(
        "Mercury", ContaminantType.Chemical,
        halfLifeSeconds: 0,                          // conservative (metal, no decay)
        partitionCoefficient: 500,                   // strong sediment binding
        toxicityThresholdMgL: 0.006,
        lethalThresholdMgL: 10);

    /// <summary>Lead (Pb²⁺). WHO guideline: 0.01 mg/L.</summary>
    public static readonly AquaticContaminant Lead = new AquaticContaminant(
        "Lead", ContaminantType.Chemical,
        halfLifeSeconds: 0,
        partitionCoefficient: 900,
        toxicityThresholdMgL: 0.01,
        lethalThresholdMgL: 50);

    /// <summary>Arsenic (As). WHO guideline: 0.01 mg/L.</summary>
    public static readonly AquaticContaminant Arsenic = new AquaticContaminant(
        "Arsenic", ContaminantType.Chemical,
        halfLifeSeconds: 0,
        partitionCoefficient: 200,
        toxicityThresholdMgL: 0.01,
        lethalThresholdMgL: 20);

    // ── Biological ───────────────────────────────────────────────

    /// <summary>Escherichia coli. Die-off t½ ≈ 2 days in surface water.</summary>
    public static readonly AquaticContaminant EColi = new AquaticContaminant(
        "EColi", ContaminantType.Biological,
        halfLifeSeconds: 2 * 86400,
        partitionCoefficient: 50,                    // moderate attachment to particles
        toxicityThresholdMgL: 0.001,                 // ~1 CFU/mL proxy
        lethalThresholdMgL: 1.0);

    /// <summary>Enterococcus. Die-off t½ ≈ 3 days in surface water.</summary>
    public static readonly AquaticContaminant Enterococcus = new AquaticContaminant(
        "Enterococcus", ContaminantType.Biological,
        halfLifeSeconds: 3 * 86400,
        partitionCoefficient: 50,
        toxicityThresholdMgL: 0.001,
        lethalThresholdMgL: 1.0);

    // ── Thermal ──────────────────────────────────────────────────

    /// <summary>
    /// Generic thermal pollution (heated effluent). No decay, no adsorption.
    /// Thresholds in °C delta above ambient.
    /// </summary>
    public static readonly AquaticContaminant GenericHeat = new AquaticContaminant(
        "GenericHeat", ContaminantType.Thermal,
        halfLifeSeconds: 0,
        partitionCoefficient: 0,
        toxicityThresholdMgL: 3.0,                   // 3 °C above ambient (fish stress)
        lethalThresholdMgL: 10.0);                   // 10 °C above ambient

    // ═══════════════════════════════════════════════════════════════
    //  Equality & display
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Determines equality by case-insensitive name comparison.</summary>
    public bool Equals(AquaticContaminant other) =>
        string.Equals(Name, other.Name, StringComparison.OrdinalIgnoreCase);

    public override bool Equals(object obj) =>
        obj is AquaticContaminant ac && Equals(ac);

    public override int GetHashCode() =>
        StringComparer.OrdinalIgnoreCase.GetHashCode(Name ?? string.Empty);

    public static bool operator ==(AquaticContaminant a, AquaticContaminant b) => a.Equals(b);
    public static bool operator !=(AquaticContaminant a, AquaticContaminant b) => !a.Equals(b);

    public override string ToString() =>
        string.Format(CultureInfo.InvariantCulture,
            "AquaticContaminant {0} ({1}, t½={2:G3} s, Kd={3:G3} L/kg, tox={4:G3} mg/L)",
            Name, Type, HalfLifeSeconds, PartitionCoefficient, ToxicityThresholdMgL);
}
