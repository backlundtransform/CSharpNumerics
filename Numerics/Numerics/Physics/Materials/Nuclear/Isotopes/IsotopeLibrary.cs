using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Materials.Nuclear.Isotopes
{
    /// <summary>
    /// Static registry of known isotopes. Provides lookup by name and
    /// supports registering custom isotopes at runtime.
    /// </summary>
    public static class IsotopeLibrary
    {
        private static readonly Dictionary<string, Isotope> _registry =
            new Dictionary<string, Isotope>(StringComparer.OrdinalIgnoreCase);

        static IsotopeLibrary()
        {
            Register(Isotope.Cs137);
            Register(Isotope.I131);
            Register(Isotope.Sr90);
            Register(Isotope.Co60);
            Register(Isotope.Am241);
            Register(Isotope.Ba137m);
            Register(Isotope.Xe131);
            Register(Isotope.Ba137);
        }

        /// <summary>
        /// Looks up an isotope by name (case-insensitive).
        /// Throws <see cref="KeyNotFoundException"/> if not found.
        /// </summary>
        public static Isotope Get(string name)
        {
            if (name == null) throw new ArgumentNullException(nameof(name));
            if (_registry.TryGetValue(name, out var iso))
                return iso;
            throw new KeyNotFoundException($"Isotope '{name}' not found in the library.");
        }

        /// <summary>
        /// Tries to look up an isotope by name.
        /// Returns true if found.
        /// </summary>
        public static bool TryGet(string name, out Isotope isotope)
        {
            if (name == null) { isotope = default; return false; }
            return _registry.TryGetValue(name, out isotope);
        }

        /// <summary>
        /// Returns all registered isotopes.
        /// </summary>
        public static IEnumerable<Isotope> All() => _registry.Values;

        /// <summary>
        /// Returns all isotopes with the given atomic number Z.
        /// </summary>
        public static IEnumerable<Isotope> ByElement(int atomicNumber)
        {
            foreach (var iso in _registry.Values)
                if (iso.AtomicNumber == atomicNumber)
                    yield return iso;
        }

        /// <summary>
        /// Registers a custom isotope. Overwrites any existing entry with the same name.
        /// </summary>
        public static void Register(Isotope isotope)
        {
            _registry[isotope.Name] = isotope;
        }
    }
}
