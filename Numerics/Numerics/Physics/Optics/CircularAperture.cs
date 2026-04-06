using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// A circular aperture (iris) that blocks rays outside a given radius.
/// Rays that pass through continue unmodified; blocked rays return no hit.
/// Implements <see cref="IOpticalSurface"/> so the ray tracer can test it.
/// </summary>
public class CircularAperture : IOpticalSurface
{
    /// <summary>Centre of the aperture.</summary>
    public Vector Center { get; }

    /// <summary>Normal to the aperture plane (unit vector).</summary>
    public Vector Normal { get; }

    /// <summary>Opening radius.</summary>
    public double Radius { get; }

    public CircularAperture(Vector center, Vector normal, double radius)
    {
        Center = center;
        Normal = normal.GetUnitVector();
        Radius = radius;
    }

    /// <inheritdoc />
    /// <remarks>
    /// Returns a hit only for rays that pass through the opening.
    /// The caller should treat a non-null hit as "ray continues".
    /// </remarks>
    public RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        double denom = ray.Direction.Dot(Normal);
        if (Math.Abs(denom) < 1e-12) return null;

        double t = (Center - ray.Origin).Dot(Normal) / denom;
        if (t < tMin || t > tMax) return null;

        Vector hitPoint = ray.PointAt(t);
        double dist = (hitPoint - Center).GetMagnitude();
        if (dist > Radius) return null; // blocked

        Vector outNormal = denom < 0 ? Normal : (-1.0) * Normal;
        return new RayHit(hitPoint, outNormal, t,
            new OpticalMedium("Air", 1.000293), this);
    }
}
