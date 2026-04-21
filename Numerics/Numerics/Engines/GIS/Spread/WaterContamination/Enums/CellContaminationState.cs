namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination.Enums;

/// <summary>
/// State of a single grid cell during a water contamination spread simulation.
/// </summary>
public enum CellContaminationState
{
    /// <summary>Cell has no contamination (concentration = 0).</summary>
    Clean = 0,

    /// <summary>Cell currently has concentration above the toxicity threshold.</summary>
    Contaminated = 1,

    /// <summary>Cell was contaminated but concentration has decayed below threshold.</summary>
    Decayed = 2,

    /// <summary>Cell is an active contaminant injection source.</summary>
    Source = 3
}
