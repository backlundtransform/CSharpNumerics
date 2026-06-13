namespace CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;

/// <summary>
/// State of a single grid cell during a wildfire spread simulation.
/// </summary>
public enum CellBurnState
{
    /// <summary>Cell has not ignited.</summary>
    Unburned = 0,

    /// <summary>Cell is actively burning.</summary>
    Burning = 1,

    /// <summary>Cell has finished burning (fuel consumed).</summary>
    Burned = 2,

    /// <summary>Cell is a firebreak and cannot burn.</summary>
    Firebreak = 3
}
