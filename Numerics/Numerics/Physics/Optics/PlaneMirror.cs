using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// An infinite flat reflective surface (perfect mirror) defined by a point and a normal.
/// All incident light is reflected; no transmission occurs.
/// </summary>
public class PlaneMirror : IOpticalSurface
{
    /// <summary>A point on the mirror plane.</summary>
    public Vector Center { get; }

    /// <summary>Outward-facing normal of the mirror surface (unit vector).</summary>
    public Vector Normal { get; }

    /// <summary>Reflectivity (0–1). Default is 1.0 (perfect mirror).</summary>
    public double Reflectivity { get; }

    public PlaneMirror(Vector center, Vector normal, double reflectivity = 1.0)
    {
        Center = center;
        Normal = normal.GetUnitVector();
        Reflectivity = reflectivity;
    }

    /// <inheritdoc />
    public RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        double denom = ray.Direction.Dot(Normal);
        if (Math.Abs(denom) < 1e-12) return null; // parallel

        double t = (Center - ray.Origin).Dot(Normal) / denom;
        if (t < tMin || t > tMax) return null;

        Vector hitPoint = ray.PointAt(t);
        // Normal always faces toward the incoming ray
        Vector outNormal = denom < 0 ? Normal : (-1.0) * Normal;
        return new RayHit(hitPoint, outNormal, t,
            new OpticalMedium("Mirror", 0), this);
    }
}
