using CSharpNumerics.Physics.Materials.Nuclear.DecayChains;
using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using System;

namespace CSharpNumerics.Physics.Materials
{
    /// <summary>
    /// Descriptor that couples an <see cref="Isotope"/> (and optional <see cref="DecayChain"/>)
    /// to a dispersion scenario so the simulator can compute activity and dose alongside concentration.
    /// <para>
    /// Created via the static factory: <c>Materials.Radioisotope("Cs137")</c>.
    /// </para>
    /// </summary>
    public class MaterialDescriptor
    {
        /// <summary>The primary isotope being dispersed.</summary>
        public Isotope Isotope { get; }

        /// <summary>
        /// Optional decay chain. When set, activity is computed for all daughters as well.
        /// </summary>
        public DecayChain Chain { get; }

        internal MaterialDescriptor(Isotope isotope, DecayChain chain = null)
        {
            Isotope = isotope;
            Chain = chain;
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

        private static DecayChain TryGetKnownChain(Isotope isotope)
        {
            if (isotope.Name == "Cs137") return DecayChain.Cs137Chain();
            if (isotope.Name == "I131") return DecayChain.I131Chain();
            return null;
        }
    }
}
