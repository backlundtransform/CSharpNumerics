using System;
using System.Globalization;

namespace CSharpNumerics.Physics.Materials.Nuclear.Isotopes
{
    /// <summary>
    /// Specifies the primary decay mode of a radioactive isotope.
    /// </summary>
    public enum DecayMode
    {
        /// <summary>Alpha decay (emission of a helium-4 nucleus).</summary>
        Alpha,

        /// <summary>Beta-minus decay (neutron → proton + electron + antineutrino).</summary>
        BetaMinus,

        /// <summary>Beta-plus / electron capture.</summary>
        BetaPlus,

        /// <summary>Gamma emission (isomeric transition).</summary>
        Gamma,

        /// <summary>Stable isotope (no decay).</summary>
        Stable
    }

    /// <summary>
    /// An immutable descriptor of a radioactive (or stable) isotope.
    /// Carries nuclear properties needed for decay, activity, and dose calculations.
    /// <para>
    /// Common isotopes are available as static properties:
    /// <c>Isotope.Cs137</c>, <c>Isotope.I131</c>, etc.
    /// </para>
    /// </summary>
    public readonly struct Isotope : IEquatable<Isotope>
    {
        /// <summary>Isotope symbol including mass number (e.g. "Cs137").</summary>
        public string Name { get; }

        /// <summary>Atomic number Z (number of protons).</summary>
        public int AtomicNumber { get; }

        /// <summary>Mass number A (protons + neutrons).</summary>
        public int MassNumber { get; }

        /// <summary>Physical half-life in seconds.</summary>
        public double HalfLife { get; }

        /// <summary>Specific activity in Bq/kg.</summary>
        public double SpecificActivity { get; }

        /// <summary>Primary decay mode.</summary>
        public DecayMode DecayMode { get; }

        /// <summary>Principal gamma energy in MeV (0 if no gamma).</summary>
        public double GammaEnergy { get; }

        /// <summary>
        /// ICRP effective dose coefficient for inhalation (Sv/Bq).
        /// Public access, default 0.
        /// </summary>
        public double InhalationDoseCoeff { get; }

        /// <summary>
        /// External dose-rate coefficient for ground-shine (Sv·m²/Bq·s).
        /// </summary>
        public double GroundShineDoseCoeff { get; }

        /// <summary>Whether this isotope is stable (infinite half-life).</summary>
        public bool IsStable => DecayMode == DecayMode.Stable || double.IsPositiveInfinity(HalfLife);

        /// <summary>Decay constant λ = ln(2) / t½  (s⁻¹). Zero for stable isotopes.</summary>
        public double Lambda => IsStable ? 0.0 : Math.Log(2) / HalfLife;

        /// <summary>
        /// Creates an isotope descriptor.
        /// </summary>
        public Isotope(
            string name,
            int atomicNumber,
            int massNumber,
            double halfLife,
            double specificActivity,
            DecayMode decayMode,
            double gammaEnergy = 0,
            double inhalationDoseCoeff = 0,
            double groundShineDoseCoeff = 0)
        {
            Name = name ?? throw new ArgumentNullException(nameof(name));
            AtomicNumber = atomicNumber;
            MassNumber = massNumber;
            HalfLife = halfLife;
            SpecificActivity = specificActivity;
            DecayMode = decayMode;
            GammaEnergy = gammaEnergy;
            InhalationDoseCoeff = inhalationDoseCoeff;
            GroundShineDoseCoeff = groundShineDoseCoeff;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Common isotopes (ICRP-72 / UNSCEAR data)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Caesium-137. t½ = 30.17 yr, β⁻, 0.662 MeV γ.</summary>
        public static readonly Isotope Cs137 = new Isotope(
            "Cs137", 55, 137,
            halfLife: 30.17 * 365.25 * 86400,       // ~9.51e8 s
            specificActivity: 3.215e12,              // Bq/kg
            decayMode: DecayMode.BetaMinus,
            gammaEnergy: 0.662,
            inhalationDoseCoeff: 3.9e-8,             // Sv/Bq (ICRP-72, adult)
            groundShineDoseCoeff: 2.1e-16);          // Sv·m²/(Bq·s)

        /// <summary>Iodine-131. t½ = 8.02 d, β⁻, 0.364 MeV γ.</summary>
        public static readonly Isotope I131 = new Isotope(
            "I131", 53, 131,
            halfLife: 8.02 * 86400,                  // ~6.93e5 s
            specificActivity: 4.6e15,
            decayMode: DecayMode.BetaMinus,
            gammaEnergy: 0.364,
            inhalationDoseCoeff: 7.4e-9,
            groundShineDoseCoeff: 1.5e-16);

        /// <summary>Strontium-90. t½ = 28.8 yr, β⁻, no γ.</summary>
        public static readonly Isotope Sr90 = new Isotope(
            "Sr90", 38, 90,
            halfLife: 28.8 * 365.25 * 86400,
            specificActivity: 5.11e12,
            decayMode: DecayMode.BetaMinus,
            gammaEnergy: 0,
            inhalationDoseCoeff: 1.6e-7,
            groundShineDoseCoeff: 0);

        /// <summary>Cobalt-60. t½ = 5.27 yr, β⁻, 1.173 + 1.333 MeV γ.</summary>
        public static readonly Isotope Co60 = new Isotope(
            "Co60", 27, 60,
            halfLife: 5.2714 * 365.25 * 86400,
            specificActivity: 4.19e13,
            decayMode: DecayMode.BetaMinus,
            gammaEnergy: 1.253,                      // average of two gammas
            inhalationDoseCoeff: 1.7e-8,
            groundShineDoseCoeff: 9.2e-16);

        /// <summary>Americium-241. t½ = 432.2 yr, α, 0.060 MeV γ.</summary>
        public static readonly Isotope Am241 = new Isotope(
            "Am241", 95, 241,
            halfLife: 432.2 * 365.25 * 86400,
            specificActivity: 1.27e11,
            decayMode: DecayMode.Alpha,
            gammaEnergy: 0.060,
            inhalationDoseCoeff: 9.6e-5,
            groundShineDoseCoeff: 1.2e-17);

        /// <summary>Barium-137m (metastable). t½ = 153 s, γ (isomeric transition), 0.662 MeV.</summary>
        public static readonly Isotope Ba137m = new Isotope(
            "Ba137m", 56, 137,
            halfLife: 153.12,
            specificActivity: 3.13e18,
            decayMode: DecayMode.Gamma,
            gammaEnergy: 0.662);

        /// <summary>Xenon-131 (stable).</summary>
        public static readonly Isotope Xe131 = new Isotope(
            "Xe131", 54, 131,
            halfLife: double.PositiveInfinity,
            specificActivity: 0,
            decayMode: DecayMode.Stable);

        /// <summary>Barium-137 (stable).</summary>
        public static readonly Isotope Ba137 = new Isotope(
            "Ba137", 56, 137,
            halfLife: double.PositiveInfinity,
            specificActivity: 0,
            decayMode: DecayMode.Stable);

        // ═══════════════════════════════════════════════════════════════
        //  Equality
        // ═══════════════════════════════════════════════════════════════

        public bool Equals(Isotope other) =>
            Name == other.Name && AtomicNumber == other.AtomicNumber && MassNumber == other.MassNumber;

        public override bool Equals(object obj) => obj is Isotope iso && Equals(iso);

        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 31 + (Name?.GetHashCode() ?? 0);
                hash = hash * 31 + AtomicNumber;
                hash = hash * 31 + MassNumber;
                return hash;
            }
        }

        public static bool operator ==(Isotope a, Isotope b) => a.Equals(b);
        public static bool operator !=(Isotope a, Isotope b) => !a.Equals(b);

        public override string ToString() =>
            string.Format(CultureInfo.InvariantCulture,
                "{0} (Z={1}, A={2}, t½={3})",
                Name, AtomicNumber, MassNumber,
                IsStable ? "stable" : FormatHalfLife(HalfLife));

        private static string FormatHalfLife(double seconds)
        {
            if (seconds < 60) return $"{seconds:F1}s";
            if (seconds < 3600) return $"{seconds / 60:F1}min";
            if (seconds < 86400) return $"{seconds / 3600:F1}h";
            if (seconds < 365.25 * 86400) return $"{seconds / 86400:F1}d";
            return string.Format(CultureInfo.InvariantCulture, "{0:F2}yr", seconds / (365.25 * 86400));
        }
    }
}
