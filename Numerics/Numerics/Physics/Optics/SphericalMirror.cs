using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Optics.Enums;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// A spherical mirror (concave or convex) defined by centre of curvature,
/// pole position, radius of curvature, and aperture radius.
/// </summary>
public class SphericalMirror : IOpticalSurface
{
    /// <summary>Centre of curvature of the sphere.</summary>
    public Vector CentreOfCurvature { get; }

    /// <summary>Radius of curvature R (always positive).</summary>
    public double RadiusOfCurvature { get; }

    /// <summary>Focal length f = R / 2. Positive for concave, negative for convex.</summary>
    public double FocalLength { get; }

    /// <summary>Mirror type.</summary>
    public MirrorType Type { get; }

    /// <summary>Half-aperture radius that limits the reflective area.</summary>
    public double ApertureRadius { get; }

    /// <summary>Reflectivity (0–1).</summary>
    public double Reflectivity { get; }

    /// <summary>
    /// Creates a spherical mirror.
    /// </summary>
    /// <param name="centreOfCurvature">Centre of the sphere whose surface forms the mirror.</param>
    /// <param name="radiusOfCurvature">R &gt; 0.</param>
    /// <param name="type">Concave or convex.</param>
    /// <param name="apertureRadius">Maximum distance from the optical axis for reflections.</param>
    /// <param name="reflectivity">0–1.</param>
    public SphericalMirror(Vector centreOfCurvature, double radiusOfCurvature,
        MirrorType type, double apertureRadius = double.MaxValue,
        double reflectivity = 1.0)
    {
        CentreOfCurvature = centreOfCurvature;
        RadiusOfCurvature = radiusOfCurvature;
        Type = type;
        ApertureRadius = apertureRadius;
        Reflectivity = reflectivity;
        FocalLength = type == MirrorType.Concave
            ? radiusOfCurvature / 2.0
            : -radiusOfCurvature / 2.0;
    }

    /// <inheritdoc />
    public RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        // Solve |O + t·D − C|² = R² for t
        Vector oc = ray.Origin - CentreOfCurvature;
        double a = ray.Direction.Dot(ray.Direction);
        double b = 2.0 * oc.Dot(ray.Direction);
        double c = oc.Dot(oc) - RadiusOfCurvature * RadiusOfCurvature;
        double disc = b * b - 4.0 * a * c;
        if (disc < 0) return null;

        double sqrtDisc = Math.Sqrt(disc);
        double t1 = (-b - sqrtDisc) / (2.0 * a);
        double t2 = (-b + sqrtDisc) / (2.0 * a);

        // For concave mirror, the reflective surface is the inner side
        // For convex mirror, the reflective surface is the outer side
        double t = Type == MirrorType.Concave ? t2 : t1;
        if (t < tMin || t > tMax)
        {
            t = Type == MirrorType.Concave ? t1 : t2;
            if (t < tMin || t > tMax) return null;
        }

        Vector hitPoint = ray.PointAt(t);

        // Aperture check
        if (ApertureRadius < double.MaxValue)
        {
            Vector fromCentre = hitPoint - CentreOfCurvature;
            double distFromAxis = fromCentre.GetMagnitude();
            if (distFromAxis > ApertureRadius) return null;
        }

        // Normal points from centre of curvature toward hit point (outward)
        Vector normal = (hitPoint - CentreOfCurvature).GetUnitVector();
        // For concave mirror, flip normal so it faces the incoming ray
        if (Type == MirrorType.Concave)
            normal = (-1.0) * normal;

        return new RayHit(hitPoint, normal, t,
            new OpticalMedium("Mirror", 0), this);
    }

    /// <summary>
    /// Mirror equation: 1/f = 1/d_o + 1/d_i. Computes image distance.
    /// </summary>
    /// <param name="objectDistance">Object distance d_o (positive in front of mirror).</param>
    /// <returns>Image distance d_i (positive = real image, negative = virtual).</returns>
    public double ImageDistance(double objectDistance) =>
        1.0 / (1.0 / FocalLength - 1.0 / objectDistance);

    /// <summary>
    /// Lateral magnification: m = −d_i / d_o.
    /// </summary>
    public double Magnification(double objectDistance)
    {
        double di = ImageDistance(objectDistance);
        return -di / objectDistance;
    }
}
