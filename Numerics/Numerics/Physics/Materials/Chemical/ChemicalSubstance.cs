using System;
using System.Globalization;

namespace CSharpNumerics.Physics.Materials.Chemical
{
    /// <summary>
    /// Physical state of a chemical substance at standard conditions (20 °C, 1 atm).
    /// </summary>
    public enum PhaseAtSTP
    {
        /// <summary>Gas at standard temperature and pressure.</summary>
        Gas,

        /// <summary>Liquid at standard temperature and pressure (may form gas cloud on release).</summary>
        Liquid,

        /// <summary>Liquefied gas stored under pressure.</summary>
        LiquefiedGas
    }

    /// <summary>
    /// An immutable descriptor of a hazardous chemical substance carrying
    /// the physical and toxicological properties needed for atmospheric
    /// dispersion and exposure modelling.
    /// <para>
    /// Common substances are available as static properties:
    /// <c>ChemicalSubstance.Chlorine</c>, <c>ChemicalSubstance.Ammonia</c>, etc.
    /// </para>
    /// </summary>
    public readonly struct ChemicalSubstance : IEquatable<ChemicalSubstance>
    {
        /// <summary>Short identifier (e.g. "Cl2", "NH3").</summary>
        public string Formula { get; }

        /// <summary>Common name (e.g. "Chlorine", "Ammonia").</summary>
        public string Name { get; }

        /// <summary>CAS registry number (e.g. "7782-50-5").</summary>
        public string CasNumber { get; }

        /// <summary>Molar mass in g/mol.</summary>
        public double MolarMass { get; }

        /// <summary>Vapour density relative to air (air = 1.0).</summary>
        public double VapourDensity { get; }

        /// <summary>Boiling point in °C at 1 atm.</summary>
        public double BoilingPointC { get; }

        /// <summary>Phase at standard temperature and pressure.</summary>
        public PhaseAtSTP Phase { get; }

        // ── Toxicological thresholds (ppm) ───────────────────────────

        /// <summary>
        /// IDLH — Immediately Dangerous to Life or Health (ppm).
        /// NIOSH value. 30-minute exposure limit.
        /// </summary>
        public double IDLH { get; }

        /// <summary>
        /// ERPG-2 — Emergency Response Planning Guideline Level 2 (ppm).
        /// Maximum airborne concentration below which nearly all persons
        /// could be exposed for up to 1 hour without experiencing
        /// irreversible or serious health effects.
        /// </summary>
        public double ERPG2 { get; }

        /// <summary>
        /// ERPG-3 — Emergency Response Planning Guideline Level 3 (ppm).
        /// Maximum airborne concentration below which nearly all persons
        /// could be exposed for up to 1 hour without experiencing
        /// life-threatening health effects.
        /// </summary>
        public double ERPG3 { get; }

        /// <summary>
        /// LC50 — Lethal Concentration 50% (ppm) for 1-hour inhalation in rats.
        /// </summary>
        public double LC50 { get; }

        /// <summary>
        /// TLV-TWA — Threshold Limit Value, Time-Weighted Average (ppm).
        /// 8-hour occupational limit.
        /// </summary>
        public double TLV_TWA { get; }

        /// <summary>
        /// TLV-STEL — Threshold Limit Value, Short-Term Exposure Limit (ppm).
        /// 15-minute ceiling.
        /// </summary>
        public double TLV_STEL { get; }

        // ── Conversion helpers ───────────────────────────────────────

        /// <summary>
        /// Converts a mass concentration (kg/m³) to ppm (parts per million by volume)
        /// at 20 °C, 1 atm using the molar volume of an ideal gas (24.04 L/mol at 20 °C).
        /// </summary>
        public double KgM3ToPpm(double kgPerM3)
        {
            // ppm = (C_kg/m³ × 10⁶ × 24.04e-3 m³/mol) / (MolarMass × 10⁻³ kg/mol)
            // Simplified: ppm = C × 24.04e6 / MolarMass
            return kgPerM3 * 24.04e6 / MolarMass;
        }

        /// <summary>
        /// Converts ppm (parts per million by volume) to kg/m³
        /// at 20 °C, 1 atm.
        /// </summary>
        public double PpmToKgM3(double ppm)
        {
            return ppm * MolarMass / 24.04e6;
        }

        /// <summary>
        /// Creates a chemical substance descriptor.
        /// </summary>
        public ChemicalSubstance(
            string formula,
            string name,
            string casNumber,
            double molarMass,
            double vapourDensity,
            double boilingPointC,
            PhaseAtSTP phase,
            double idlh,
            double erpg2,
            double erpg3,
            double lc50,
            double tlvTwa,
            double tlvStel)
        {
            Formula = formula ?? throw new ArgumentNullException(nameof(formula));
            Name = name ?? throw new ArgumentNullException(nameof(name));
            CasNumber = casNumber;
            MolarMass = molarMass;
            VapourDensity = vapourDensity;
            BoilingPointC = boilingPointC;
            Phase = phase;
            IDLH = idlh;
            ERPG2 = erpg2;
            ERPG3 = erpg3;
            LC50 = lc50;
            TLV_TWA = tlvTwa;
            TLV_STEL = tlvStel;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Built-in substances
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Chlorine (Cl₂) — toxic yellow-green gas, heavier than air.</summary>
        public static readonly ChemicalSubstance Chlorine = new ChemicalSubstance(
            formula: "Cl2",
            name: "Chlorine",
            casNumber: "7782-50-5",
            molarMass: 70.906,
            vapourDensity: 2.49,
            boilingPointC: -34.04,
            phase: PhaseAtSTP.Gas,
            idlh: 10,
            erpg2: 3,
            erpg3: 20,
            lc50: 293,
            tlvTwa: 0.5,
            tlvStel: 1.0);

        /// <summary>Ammonia (NH₃) — toxic colourless gas, lighter than air.</summary>
        public static readonly ChemicalSubstance Ammonia = new ChemicalSubstance(
            formula: "NH3",
            name: "Ammonia",
            casNumber: "7664-41-7",
            molarMass: 17.031,
            vapourDensity: 0.59,
            boilingPointC: -33.34,
            phase: PhaseAtSTP.Gas,
            idlh: 300,
            erpg2: 200,
            erpg3: 1000,
            lc50: 2000,
            tlvTwa: 25,
            tlvStel: 35);

        /// <summary>Hydrogen sulfide (H₂S) — highly toxic, heavier than air, rotten-egg odour.</summary>
        public static readonly ChemicalSubstance HydrogenSulfide = new ChemicalSubstance(
            formula: "H2S",
            name: "Hydrogen Sulfide",
            casNumber: "7783-06-4",
            molarMass: 34.081,
            vapourDensity: 1.19,
            boilingPointC: -60.28,
            phase: PhaseAtSTP.Gas,
            idlh: 50,
            erpg2: 30,
            erpg3: 100,
            lc50: 444,
            tlvTwa: 1,
            tlvStel: 5);

        /// <summary>Methane (CH₄) — flammable, lighter than air, primary component of natural gas.</summary>
        public static readonly ChemicalSubstance Methane = new ChemicalSubstance(
            formula: "CH4",
            name: "Methane",
            casNumber: "74-82-8",
            molarMass: 16.043,
            vapourDensity: 0.55,
            boilingPointC: -161.5,
            phase: PhaseAtSTP.Gas,
            idlh: double.PositiveInfinity,   // simple asphyxiant — no IDLH for toxicity
            erpg2: 5000,                     // lower explosive limit concern
            erpg3: 50000,
            lc50: double.PositiveInfinity,   // not acutely toxic
            tlvTwa: 1000,
            tlvStel: double.PositiveInfinity);

        /// <summary>Propane (C₃H₈) — flammable liquefied gas, heavier than air (LPG).</summary>
        public static readonly ChemicalSubstance Propane = new ChemicalSubstance(
            formula: "C3H8",
            name: "Propane",
            casNumber: "74-98-6",
            molarMass: 44.096,
            vapourDensity: 1.52,
            boilingPointC: -42.1,
            phase: PhaseAtSTP.LiquefiedGas,
            idlh: 2100,
            erpg2: 5000,
            erpg3: 17000,
            lc50: double.PositiveInfinity,   // simple asphyxiant
            tlvTwa: 1000,
            tlvStel: double.PositiveInfinity);

        // ═══════════════════════════════════════════════════════════════
        //  Equality & display
        // ═══════════════════════════════════════════════════════════════

        public bool Equals(ChemicalSubstance other) =>
            string.Equals(Formula, other.Formula, StringComparison.OrdinalIgnoreCase);

        public override bool Equals(object obj) =>
            obj is ChemicalSubstance cs && Equals(cs);

        public override int GetHashCode() =>
            StringComparer.OrdinalIgnoreCase.GetHashCode(Formula ?? "");

        public static bool operator ==(ChemicalSubstance a, ChemicalSubstance b) => a.Equals(b);
        public static bool operator !=(ChemicalSubstance a, ChemicalSubstance b) => !a.Equals(b);

        public override string ToString() =>
            string.Format(CultureInfo.InvariantCulture,
                "{0} ({1}), MW={2:F1} g/mol, IDLH={3} ppm",
                Name, Formula, MolarMass,
                double.IsPositiveInfinity(IDLH) ? "N/A" : IDLH.ToString(CultureInfo.InvariantCulture));
    }
}
