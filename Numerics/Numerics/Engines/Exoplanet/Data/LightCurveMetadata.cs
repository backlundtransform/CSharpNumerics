using CSharpNumerics.Engines.Exoplanet.Enums;

namespace CSharpNumerics.Engines.Exoplanet.Data;

public class LightCurveMetadata
{
    public string TargetId { get; set; } = string.Empty;
    public string Mission { get; set; } = string.Empty;
    public double TimeOffset { get; set; }
    public CadenceType Cadence { get; set; }

    public LightCurveMetadata() { }

    public LightCurveMetadata(string targetId, string mission, double timeOffset, CadenceType cadence)
    {
        TargetId = targetId;
        Mission = mission;
        TimeOffset = timeOffset;
        Cadence = cadence;
    }
}
