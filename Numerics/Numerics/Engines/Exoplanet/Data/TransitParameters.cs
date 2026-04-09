namespace CSharpNumerics.Engines.Exoplanet.Data;

public class TransitParameters
{
    public double Period { get; set; }
    public double Epoch { get; set; }
    public double Depth { get; set; }
    public double Duration { get; set; }
    public double RadiusRatio { get; set; }
    public double ImpactParameter { get; set; }
    public double IngressDuration { get; set; }

    public TransitParameters() { }

    public TransitParameters(double period, double epoch, double depth, double duration,
        double radiusRatio, double impactParameter, double ingressDuration)
    {
        Period = period;
        Epoch = epoch;
        Depth = depth;
        Duration = duration;
        RadiusRatio = radiusRatio;
        ImpactParameter = impactParameter;
        IngressDuration = ingressDuration;
    }
}
