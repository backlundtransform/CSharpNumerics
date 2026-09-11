namespace CSharpNumerics.Physics.Environmental.Enums;

/// <summary>
/// Standard t-squared fire growth classes, defined by how long the fire takes to
/// reach 1055 kW (NFPA 72, NFPA 92, EN 1991-1-2).
/// </summary>
public enum FireGrowthRate
{
    /// <summary>Reaches 1055 kW in 600 s. Densely packed paper, some furnishings.</summary>
    Slow,

    /// <summary>Reaches 1055 kW in 300 s. Cabins, offices, shop displays.</summary>
    Medium,

    /// <summary>Reaches 1055 kW in 150 s. Upholstered furniture, stacked cartons.</summary>
    Fast,

    /// <summary>Reaches 1055 kW in 75 s. Pool fires, flammable liquids, foam plastics.</summary>
    UltraFast
}
