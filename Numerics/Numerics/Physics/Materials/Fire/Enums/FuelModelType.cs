namespace CSharpNumerics.Physics.Materials.Fire.Enums;

/// <summary>
/// Anderson 13 standard fire behaviour fuel models.
/// Values match the original NFFL model numbers (1–13).
/// </summary>
public enum FuelModelType
{
    /// <summary>No burnable fuel.</summary>
    NoFuel = 0,

    /// <summary>Model 1 — Short grass (1 ft).</summary>
    ShortGrass = 1,

    /// <summary>Model 2 — Timber (grass and understory).</summary>
    TimberGrassUnderstory = 2,

    /// <summary>Model 3 — Tall grass (2.5 ft).</summary>
    TallGrass = 3,

    /// <summary>Model 4 — Chaparral (6 ft).</summary>
    Chaparral = 4,

    /// <summary>Model 5 — Brush (2 ft).</summary>
    Brush = 5,

    /// <summary>Model 6 — Dormant brush, hardwood slash.</summary>
    DormantBrush = 6,

    /// <summary>Model 7 — Southern rough.</summary>
    SouthernRough = 7,

    /// <summary>Model 8 — Closed timber litter.</summary>
    ClosedTimberLitter = 8,

    /// <summary>Model 9 — Hardwood litter.</summary>
    HardwoodLitter = 9,

    /// <summary>Model 10 — Timber (litter and understory).</summary>
    TimberLitterUnderstory = 10,

    /// <summary>Model 11 — Light logging slash.</summary>
    LightLoggingSlash = 11,

    /// <summary>Model 12 — Medium logging slash.</summary>
    MediumLoggingSlash = 12,

    /// <summary>Model 13 — Heavy logging slash.</summary>
    HeavyLoggingSlash = 13
}
