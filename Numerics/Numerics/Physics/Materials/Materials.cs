using CSharpNumerics.Physics.Materials.Biological;
using CSharpNumerics.Physics.Materials.Chemical;
using CSharpNumerics.Physics.Materials.Nuclear.DecayChains;
using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using CSharpNumerics.Physics.Materials.Optical;
using CSharpNumerics.Physics.Optics;
using System;

namespace CSharpNumerics.Physics.Materials
{
    /// <summary>
    /// Describes a material attached to a dispersion scenario.
    /// Can represent either a radioactive isotope (with optional decay chain)
    /// or a hazardous chemical substance — or both.
    /// <para>
    /// Created via the static factory: <c>Materials.Radioisotope("Cs137")</c>
    /// or <c>Materials.Chemical("Cl2")</c>.
    /// </para>
    /// </summary>
    public class MaterialDescriptor
    {
        /// <summary>The primary isotope being dispersed. Null for chemical-only materials.</summary>
        public Isotope Isotope { get; }

        /// <summary>
        /// Optional decay chain. When set, activity is computed for all daughters as well.
        /// </summary>
        public DecayChain Chain { get; }

        /// <summary>
        /// The chemical substance being dispersed. Null for nuclear-only materials.
        /// When set, the simulator computes ppm concentrations and toxic-dose layers.
        /// </summary>
        public ChemicalSubstance? Substance { get; }

        /// <summary>The biological agent being dispersed. Null for non-biological materials.
        /// When set, the simulator computes bioUnits, viableBioUnits, and infectiousDose layers.
        /// </summary>
        public BiologicalAgent? BiologicalAgent { get; }

        /// <summary>Whether this descriptor contains a radioactive isotope.</summary>
        public bool IsNuclear => Isotope.Name != null;

        /// <summary>Whether this descriptor contains a chemical substance.</summary>
        public bool IsChemical => Substance.HasValue;

        /// <summary>Whether this descriptor contains a biological agent.</summary>
        public bool IsBiological => BiologicalAgent.HasValue;

        internal MaterialDescriptor(Isotope isotope, DecayChain chain = null)
        {
            Isotope = isotope;
            Chain = chain;
        }

        internal MaterialDescriptor(ChemicalSubstance substance)
        {
            Substance = substance;
        }

        internal MaterialDescriptor(Isotope isotope, DecayChain chain, ChemicalSubstance substance)
        {
            Isotope = isotope;
            Chain = chain;
            Substance = substance;
        }

        internal MaterialDescriptor(BiologicalAgent agent)
        {
            BiologicalAgent = agent;
        }
    }

    /// <summary>
    /// Static factory for creating <see cref="MaterialDescriptor"/> instances to
    /// attach to dispersion scenarios.
    /// <para>
    /// Usage: <c>Materials.Radioisotope("Cs137")</c>
    /// </para>
    /// </summary>
    public static class Materials
    {
        /// <summary>
        /// Creates a material descriptor for the named radioisotope.
        /// Looks up the isotope in <see cref="IsotopeLibrary"/> and
        /// attaches a well-known decay chain if available.
        /// </summary>
        /// <param name="name">Isotope name (e.g. "Cs137", "I131"). Case-insensitive.</param>
        public static MaterialDescriptor Radioisotope(string name)
        {
            if (name == null) throw new ArgumentNullException(nameof(name));

            var isotope = IsotopeLibrary.Get(name);
            var chain = TryGetKnownChain(isotope);
            return new MaterialDescriptor(isotope, chain);
        }

        /// <summary>
        /// Creates a material descriptor from an explicit isotope instance.
        /// </summary>
        public static MaterialDescriptor Radioisotope(Isotope isotope)
        {
            var chain = TryGetKnownChain(isotope);
            return new MaterialDescriptor(isotope, chain);
        }

        /// <summary>
        /// Creates a material descriptor with a custom decay chain.
        /// </summary>
        public static MaterialDescriptor Radioisotope(Isotope isotope, DecayChain chain)
        {
            return new MaterialDescriptor(isotope, chain);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Chemical substances
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Creates a material descriptor for the named chemical substance.
        /// Looks up the substance in <see cref="ChemicalLibrary"/>.
        /// </summary>
        /// <param name="formula">Chemical formula (e.g. "Cl2", "NH3"). Case-insensitive.</param>
        public static MaterialDescriptor Chemical(string formula)
        {
            if (formula == null) throw new ArgumentNullException(nameof(formula));
            var substance = ChemicalLibrary.Get(formula);
            return new MaterialDescriptor(substance);
        }

        /// <summary>
        /// Creates a material descriptor from an explicit chemical substance instance.
        /// </summary>
        public static MaterialDescriptor Chemical(ChemicalSubstance substance)
        {
            return new MaterialDescriptor(substance);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Biological agents
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Creates a material descriptor for the named biological agent.
        /// Looks up the agent in <see cref="BiologicalLibrary"/>.
        /// </summary>
        /// <param name="name">Agent code or alias (e.g. "Virus", "Bacteria"). Case-insensitive.</param>
        public static MaterialDescriptor Biological(string name)
        {
            if (name == null) throw new ArgumentNullException(nameof(name));
            var agent = BiologicalLibrary.Get(name);
            return new MaterialDescriptor(agent);
        }

        /// <summary>
        /// Creates a material descriptor from an explicit biological agent instance.
        /// </summary>
        public static MaterialDescriptor Biological(BiologicalAgent agent)
        {
            return new MaterialDescriptor(agent);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Optical materials
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Returns an <see cref="OpticalMedium"/> by name from the optical library.
        /// </summary>
        /// <param name="name">Material name (e.g. "CrownGlass", "Diamond"). Case-insensitive.</param>
        public static OpticalMedium Optical(string name)
        {
            if (name == null) throw new ArgumentNullException(nameof(name));
            return OpticalLibrary.Get(name);
        }

        private static DecayChain TryGetKnownChain(Isotope isotope)
        {
            if (isotope.Name == "Cs137") return DecayChain.Cs137Chain();
            if (isotope.Name == "I131") return DecayChain.I131Chain();
            return null;
        }
    }
}
