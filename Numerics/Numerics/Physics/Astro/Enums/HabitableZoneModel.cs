namespace CSharpNumerics.Physics.Astro.Enums;

/// <summary>
/// Pair of <see cref="HabitableZoneLimit"/> boundaries that define a habitable zone.
/// </summary>
public enum HabitableZoneModel
{
    /// <summary>Runaway Greenhouse → Maximum Greenhouse.</summary>
    Conservative,

    /// <summary>Recent Venus → Early Mars (used by the PHL Habitable Worlds Catalog).</summary>
    Optimistic
}
