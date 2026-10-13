namespace CSharpNumerics.Physics.Astro.Enums;

/// <summary>
/// Habitable zone boundaries from Kopparapu et al. (2014), ordered from the
/// innermost (highest stellar flux) to the outermost (lowest stellar flux).
/// </summary>
public enum HabitableZoneLimit
{
    /// <summary>Optimistic inner edge: Venus had no surface water ~1 Gyr ago.</summary>
    RecentVenus,

    /// <summary>Conservative inner edge: onset of a runaway greenhouse (1 Earth mass).</summary>
    RunawayGreenhouse,

    /// <summary>Conservative outer edge: CO₂ greenhouse warming reaches its maximum.</summary>
    MaximumGreenhouse,

    /// <summary>Optimistic outer edge: Mars may have had surface water ~3.8 Gyr ago.</summary>
    EarlyMars
}
