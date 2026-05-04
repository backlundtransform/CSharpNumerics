namespace CSharpNumerics.Engines.GIS.Spread.VolumetricContamination.Enums;

/// <summary>
/// Cell state for 3D volumetric contamination simulations.
/// </summary>
public enum ContaminationCellState3D
{
    /// <summary>No contamination detected.</summary>
    Clean = 0,

    /// <summary>Concentration above toxicity threshold.</summary>
    Contaminated = 1,

    /// <summary>Active source injector cell.</summary>
    Source = 2,

    /// <summary>Impermeable land/seabed barrier.</summary>
    Land = 3
}
