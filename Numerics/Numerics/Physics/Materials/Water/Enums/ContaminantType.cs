namespace CSharpNumerics.Physics.Materials.Water.Enums;

/// <summary>
/// Classification of aquatic contaminants by their primary hazard mechanism.
/// </summary>
public enum ContaminantType
{
    /// <summary>Radioactive isotopes (e.g. Cs-137, Sr-90, I-131). Decay produces ionising radiation.</summary>
    Radioactive = 0,

    /// <summary>Chemical pollutants (e.g. Benzene, Mercury, Cyanide). Toxicity by ingestion or contact.</summary>
    Chemical = 1,

    /// <summary>Biological agents (e.g. E. coli, Enterococcus). Pathogenic organisms.</summary>
    Biological = 2,

    /// <summary>Thermal pollution (heated effluent). Measured as temperature delta (°C).</summary>
    Thermal = 3
}
