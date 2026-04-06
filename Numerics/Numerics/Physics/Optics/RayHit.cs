using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Describes the result of a ray intersecting an optical surface.
/// </summary>
public readonly struct RayHit
{
    /// <summary>The world-space point where the ray hit the surface.</summary>
    public Vector Point { get; }

    /// <summary>Outward-facing surface normal at the hit point (unit length).</summary>
    public Vector Normal { get; }

    /// <summary>Distance from the ray origin to the hit point.</summary>
    public double Distance { get; }

    /// <summary>The optical medium on the far side of the surface (the medium the ray enters).</summary>
    public OpticalMedium MediumBeyond { get; }

    /// <summary>The surface that was hit.</summary>
    public IOpticalSurface Surface { get; }

    public RayHit(Vector point, Vector normal, double distance,
        OpticalMedium mediumBeyond, IOpticalSurface surface)
    {
        Point = point;
        Normal = normal;
        Distance = distance;
        MediumBeyond = mediumBeyond;
        Surface = surface;
    }
}
