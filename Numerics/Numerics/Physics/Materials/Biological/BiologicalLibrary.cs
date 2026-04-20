using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Materials.Biological;

/// <summary>
/// Static registry of known <see cref="BiologicalAgent"/> descriptors.
/// Provides generic virus, bacteria, and spore aerosols out of the box.
/// </summary>
public static class BiologicalLibrary
{
    private static readonly Dictionary<string, BiologicalAgent> _registry =
        new Dictionary<string, BiologicalAgent>(StringComparer.OrdinalIgnoreCase);

    static BiologicalLibrary()
    {
        Register(BiologicalAgent.GenericVirus, "virus", "viruses");
        Register(BiologicalAgent.GenericBacteria, "bacteria", "bacterium");
        Register(BiologicalAgent.GenericSpore, "spore", "spores");
    }

    /// <summary>
    /// Retrieves a biological agent by code or alias.
    /// </summary>
    public static BiologicalAgent Get(string name)
    {
        if (name == null) throw new ArgumentNullException(nameof(name));
        if (_registry.TryGetValue(name, out var agent))
            return agent;

        throw new KeyNotFoundException($"Biological agent '{name}' not found in library.");
    }

    /// <summary>
    /// Tries to retrieve a biological agent by code or alias.
    /// </summary>
    public static bool TryGet(string name, out BiologicalAgent agent)
    {
        if (name != null && _registry.TryGetValue(name, out agent))
            return true;

        agent = default;
        return false;
    }

    /// <summary>Returns all registered biological agents without alias duplicates.</summary>
    public static IReadOnlyCollection<BiologicalAgent> All()
    {
        return _registry.Values
            .GroupBy(agent => agent.Code, StringComparer.OrdinalIgnoreCase)
            .Select(group => group.First())
            .ToList()
            .AsReadOnly();
    }

    /// <summary>
    /// Registers a custom biological agent and optional aliases.
    /// Existing entries are overwritten.
    /// </summary>
    public static void Register(BiologicalAgent agent, params string[] aliases)
    {
        _registry[agent.Code] = agent;
        _registry[agent.Name] = agent;

        if (aliases == null)
            return;

        for (int i = 0; i < aliases.Length; i++)
        {
            if (!string.IsNullOrWhiteSpace(aliases[i]))
                _registry[aliases[i]] = agent;
        }
    }
}
