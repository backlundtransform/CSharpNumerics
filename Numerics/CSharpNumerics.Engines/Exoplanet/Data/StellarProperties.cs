using CSharpNumerics.Engines.Exoplanet.Enums;

namespace CSharpNumerics.Engines.Exoplanet.Data;

public class StellarProperties
{
    public double EffectiveTemp { get; set; }
    public double Radius { get; set; }
    public double Mass { get; set; }
    public double SurfaceGravity { get; set; }
    public double Metallicity { get; set; }
    public SpectralType Type { get; set; }

    public StellarProperties() { }

    public StellarProperties(double effectiveTemp, double radius, double mass,
        double surfaceGravity, double metallicity, SpectralType type)
    {
        EffectiveTemp = effectiveTemp;
        Radius = radius;
        Mass = mass;
        SurfaceGravity = surfaceGravity;
        Metallicity = metallicity;
        Type = type;
    }
}
