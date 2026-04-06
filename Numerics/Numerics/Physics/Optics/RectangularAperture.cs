using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// A rectangular aperture (slit) that blocks rays outside a width × height opening.
/// </summary>
public class RectangularAperture : IOpticalSurface
{
    /// <summary>Centre of the aperture.</summary>
    public Vector Center { get; }

    /// <summary>Normal to the aperture plane (unit vector).</summary>
    public Vector Normal { get; }

    /// <summary>Local "right" direction on the aperture plane (unit vector).</summary>
    public Vector Right { get; }

    /// <summary>Local "up" direction on the aperture plane (unit vector).</summary>
    public Vector Up { get; }

    /// <summary>Half-width of the opening.</summary>
    public double HalfWidth { get; }

    /// <summary>Half-height of the opening.</summary>
    public double HalfHeight { get; }

    /// <summary>
    /// Creates a rectangular aperture.
    /// </summary>
    /// <param name="center">Centre position.</param>
    /// <param name="normal">Plane normal.</param>
    /// <param name="right">Local horizontal direction on the plane.</param>
    /// <param name="width">Full width of the opening.</param>
    /// <param name="height">Full height of the opening.</param>
    public RectangularAperture(Vector center, Vector normal, Vector right,
        double width, double height)
    {
        Center = center;
        Normal = normal.GetUnitVector();
        Right = right.GetUnitVector();
        Up = Normal.Cross(Right).GetUnitVector();
        HalfWidth = width / 2.0;
        HalfHeight = height / 2.0;
    }

    /// <inheritdoc />
    public RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        double denom = ray.Direction.Dot(Normal);
        if (Math.Abs(denom) < 1e-12) return null;

        double t = (Center - ray.Origin).Dot(Normal) / denom;
        if (t < tMin || t > tMax) return null;

        Vector hitPoint = ray.PointAt(t);
        Vector offset = hitPoint - Center;
        double u = Math.Abs(offset.Dot(Right));
        double v = Math.Abs(offset.Dot(Up));
        if (u > HalfWidth || v > HalfHeight) return null; // blocked

        Vector outNormal = denom < 0 ? Normal : (-1.0) * Normal;
        return new RayHit(hitPoint, outNormal, t,
            new OpticalMedium("Air", 1.000293), this);
    }
}
