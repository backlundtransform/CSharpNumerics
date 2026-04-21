using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Materials.Water;

/// <summary>
/// Static registry of known <see cref="AquaticContaminant"/> instances.
/// Pre-registers common radioactive, chemical, biological, and thermal contaminants.
/// Custom contaminants can be added at runtime via <see cref="Register"/>.
/// </summary>
public static class ContaminantLibrary
{
    private static readonly Dictionary<string, AquaticContaminant> _registry =
        new Dictionary<string, AquaticContaminant>(StringComparer.OrdinalIgnoreCase);

    static ContaminantLibrary()
    {
        // Radioactive
        Register(AquaticContaminant.Cs137);
        Register(AquaticContaminant.Sr90);
        Register(AquaticContaminant.I131);

        // Chemical
        Register(AquaticContaminant.Benzene);
        Register(AquaticContaminant.Toluene);
        Register(AquaticContaminant.Cyanide);
        Register(AquaticContaminant.Mercury);
        Register(AquaticContaminant.Lead);
        Register(AquaticContaminant.Arsenic);

        // Biological
        Register(AquaticContaminant.EColi);
        Register(AquaticContaminant.Enterococcus);

        // Thermal
        Register(AquaticContaminant.GenericHeat);
    }

    /// <summary>
    /// Retrieves an aquatic contaminant by name (e.g. "Cs137", "Benzene").
    /// Case-insensitive. Throws <see cref="KeyNotFoundException"/> if not found.
    /// </summary>
    public static AquaticContaminant Get(string name)
    {
        if (name == null) throw new ArgumentNullException(nameof(name));
        if (_registry.TryGetValue(name, out var contaminant))
            return contaminant;
        throw new KeyNotFoundException($"Aquatic contaminant '{name}' not found in library.");
    }

    /// <summary>
    /// Tries to retrieve an aquatic contaminant by name.
    /// Returns true if found.
    /// </summary>
    public static bool TryGet(string name, out AquaticContaminant contaminant)
    {
        if (name != null && _registry.TryGetValue(name, out contaminant))
            return true;
        contaminant = default;
        return false;
    }

    /// <summary>Returns all registered aquatic contaminants.</summary>
    public static IReadOnlyCollection<AquaticContaminant> All() =>
        _registry.Values.ToList().AsReadOnly();

    /// <summary>
    /// Registers a custom aquatic contaminant. Overwrites any existing
    /// entry with the same name.
    /// </summary>
    public static void Register(AquaticContaminant contaminant)
    {
        _registry[contaminant.Name] = contaminant;
    }
}
