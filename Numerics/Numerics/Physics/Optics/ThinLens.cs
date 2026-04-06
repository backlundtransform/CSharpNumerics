using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Optics.Enums;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Ideal thin lens modelled as a refracting plane that obeys the thin-lens equation
/// 1/f = 1/d_o + 1/d_i. Rays are refracted at the lens plane.
/// </summary>
public class ThinLens : IOpticalSurface
{
    /// <summary>Centre of the lens.</summary>
    public Vector Center { get; }

    /// <summary>Optical axis direction (unit vector, perpendicular to the lens plane).</summary>
    public Vector Axis { get; }

    /// <summary>Focal length. Positive for converging, negative for diverging.</summary>
    public double FocalLength { get; }

    /// <summary>Lens type.</summary>
    public LensType Type { get; }

    /// <summary>Aperture radius of the lens.</summary>
    public double ApertureRadius { get; }

    /// <summary>Refractive index of lens material (used for lensmaker's equation).</summary>
    public double RefractiveIndex { get; }

    /// <summary>
    /// Creates a thin lens.
    /// </summary>
    /// <param name="center">Centre position.</param>
    /// <param name="axis">Optical axis (perpendicular to lens plane). Will be normalised.</param>
    /// <param name="focalLength">Focal length |f|. Sign is determined by <paramref name="type"/>.</param>
    /// <param name="type">Converging or diverging.</param>
    /// <param name="apertureRadius">Maximum lens radius.</param>
    /// <param name="refractiveIndex">Refractive index of the glass.</param>
    public ThinLens(Vector center, Vector axis, double focalLength,
        LensType type = LensType.Converging, double apertureRadius = double.MaxValue,
        double refractiveIndex = 1.5168)
    {
        Center = center;
        Axis = axis.GetUnitVector();
        FocalLength = type == LensType.Converging ? Math.Abs(focalLength) : -Math.Abs(focalLength);
        Type = type;
        ApertureRadius = apertureRadius;
        RefractiveIndex = refractiveIndex;
    }

    /// <inheritdoc />
    public RayHit? Intersect(Ray ray, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        // Intersect with the lens plane
        double denom = ray.Direction.Dot(Axis);
        if (Math.Abs(denom) < 1e-12) return null; // parallel

        double t = (Center - ray.Origin).Dot(Axis) / denom;
        if (t < tMin || t > tMax) return null;

        Vector hitPoint = ray.PointAt(t);

        // Aperture check
        Vector offset = hitPoint - Center;
        double distFromAxis = (offset - offset.Dot(Axis) * Axis).GetMagnitude();
        if (distFromAxis > ApertureRadius) return null;

        Vector outNormal = denom < 0 ? Axis : (-1.0) * Axis;
        return new RayHit(hitPoint, outNormal, t,
            new OpticalMedium("Lens", RefractiveIndex), this);
    }

    /// <summary>
    /// Thin-lens equation: 1/f = 1/d_o + 1/d_i. Returns image distance.
    /// </summary>
    /// <param name="objectDistance">Object distance from lens (positive).</param>
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

    /// <summary>
    /// Lensmaker's equation for a thin lens in air:
    /// 1/f = (n − 1) · (1/R₁ − 1/R₂).
    /// </summary>
    /// <param name="n">Refractive index of lens material.</param>
    /// <param name="r1">Radius of curvature of the first surface (positive = convex toward light).</param>
    /// <param name="r2">Radius of curvature of the second surface (negative = convex toward light).</param>
    public static double LensmakerFocalLength(double n, double r1, double r2) =>
        1.0 / ((n - 1.0) * (1.0 / r1 - 1.0 / r2));

    /// <summary>
    /// Refract a ray through the thin lens using the paraxial approximation.
    /// A ray parallel to the axis is redirected through the focal point, etc.
    /// </summary>
    public Ray RefractRay(Ray incoming)
    {
        var hit = Intersect(incoming);
        if (!hit.HasValue) return incoming;

        Vector p = hit.Value.Point;
        // Decompose incoming direction into axial and transverse
        double axialComponent = incoming.Direction.Dot(Axis);
        Vector axialPart = axialComponent * Axis;
        Vector transversePart = incoming.Direction - axialPart;

        // Offset from optical axis on the lens
        Vector offset = p - Center;
        Vector heightOnLens = offset - offset.Dot(Axis) * Axis;

        // Thin lens formula: the lens adds a transverse kick = -h / f
        // (in the direction from the hit toward the axis)
        Vector deflection = (-1.0 / FocalLength) * heightOnLens;

        // Signum of axial component determines the propagation direction along the axis
        double sign = axialComponent >= 0 ? 1.0 : -1.0;
        Vector newDirection = (sign * Axis + transversePart + deflection).GetUnitVector();

        return new Ray(p, newDirection, incoming.WavelengthNm, incoming.Intensity);
    }
}
