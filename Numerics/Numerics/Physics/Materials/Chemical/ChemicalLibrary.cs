using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Materials.Chemical
{
    /// <summary>
    /// Static registry of known <see cref="ChemicalSubstance"/> instances.
    /// Pre-registers Cl₂, NH₃, H₂S, CH₄, C₃H₈. Custom substances can
    /// be added at runtime via <see cref="Register"/>.
    /// </summary>
    public static class ChemicalLibrary
    {
        private static readonly Dictionary<string, ChemicalSubstance> _registry =
            new Dictionary<string, ChemicalSubstance>(StringComparer.OrdinalIgnoreCase);

        static ChemicalLibrary()
        {
            Register(ChemicalSubstance.Chlorine);
            Register(ChemicalSubstance.Ammonia);
            Register(ChemicalSubstance.HydrogenSulfide);
            Register(ChemicalSubstance.Methane);
            Register(ChemicalSubstance.Propane);
        }

        /// <summary>
        /// Retrieves a chemical substance by formula (e.g. "Cl2", "NH3").
        /// Case-insensitive. Throws <see cref="KeyNotFoundException"/> if not found.
        /// </summary>
        public static ChemicalSubstance Get(string formula)
        {
            if (formula == null) throw new ArgumentNullException(nameof(formula));
            if (_registry.TryGetValue(formula, out var substance))
                return substance;
            throw new KeyNotFoundException($"Chemical substance '{formula}' not found in library.");
        }

        /// <summary>
        /// Tries to retrieve a chemical substance by formula.
        /// Returns true if found.
        /// </summary>
        public static bool TryGet(string formula, out ChemicalSubstance substance)
        {
            if (formula != null && _registry.TryGetValue(formula, out substance))
                return true;
            substance = default;
            return false;
        }

        /// <summary>Returns all registered chemical substances.</summary>
        public static IReadOnlyCollection<ChemicalSubstance> All() =>
            _registry.Values.ToList().AsReadOnly();

        /// <summary>
        /// Registers a custom chemical substance. Overwrites any existing
        /// entry with the same formula.
        /// </summary>
        public static void Register(ChemicalSubstance substance)
        {
            _registry[substance.Formula] = substance;
        }
    }
}
