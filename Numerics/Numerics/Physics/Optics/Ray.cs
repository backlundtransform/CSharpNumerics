using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Represents a light ray with origin, direction, wavelength and intensity.
/// </summary>
public struct Ray
{
    /// <summary>Origin point of the ray in world space.</summary>
    public Vector Origin;

    /// <summary>Normalised direction of travel.</summary>
    public Vector Direction;

    /// <summary>Wavelength in nanometres (default 550 nm — green light).</summary>
    public double WavelengthNm;

    /// <summary>Relative intensity (0–1). Decreases with each interaction.</summary>
    public double Intensity;

    /// <summary>
    /// Creates a new ray.
    /// </summary>
    /// <param name="origin">Starting point.</param>
    /// <param name="direction">Direction of travel (will be normalised).</param>
    /// <param name="wavelengthNm">Wavelength in nm (default 550).</param>
    /// <param name="intensity">Initial intensity 0–1 (default 1).</param>
    public Ray(Vector origin, Vector direction, double wavelengthNm = 550.0, double intensity = 1.0)
    {
        Origin = origin;
        Direction = direction.GetUnitVector();
        WavelengthNm = wavelengthNm;
        Intensity = intensity;
    }

    /// <summary>Returns the point along the ray at parameter t: Origin + t * Direction.</summary>
    public Vector PointAt(double t) => Origin + t * Direction;

    public override string ToString() =>
        $"Ray(origin={Origin}, dir={Direction}, λ={WavelengthNm} nm, I={Intensity:F3})";
}
