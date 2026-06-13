using CSharpNumerics.Engines.Game.Objects;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Optics;
using System;

namespace CSharpNumerics.Engines.Game;

/// <summary>
/// Extension methods to cast optics <see cref="Ray"/>s against game-engine
/// bounding volumes (<see cref="AABB"/>, <see cref="BoundingSphere"/>).
/// Useful for game-engine light/visibility queries.
/// </summary>
public static class RaycastExtensions
{
    /// <summary>
    /// Tests whether a ray intersects an axis-aligned bounding box.
    /// Uses the slab method.
    /// </summary>
    /// <param name="ray">The optics ray to test.</param>
    /// <param name="aabb">The AABB to test against.</param>
    /// <param name="tMin">Minimum ray parameter.</param>
    /// <param name="tMax">Maximum ray parameter.</param>
    /// <returns>
    /// The distance along the ray to the nearest intersection,
    /// or null if no intersection occurs.
    /// </returns>
    public static double? IntersectAABB(this Ray ray, AABB aabb,
        double tMin = 1e-6, double tMax = double.MaxValue)
    {
        double tNear = tMin;
        double tFar = tMax;

        // X slab
        if (!SlabTest(ray.Origin.x, ray.Direction.x, aabb.Min.x, aabb.Max.x, ref tNear, ref tFar))
            return null;
        // Y slab
        if (!SlabTest(ray.Origin.y, ray.Direction.y, aabb.Min.y, aabb.Max.y, ref tNear, ref tFar))
            return null;
        // Z slab
        if (!SlabTest(ray.Origin.z, ray.Direction.z, aabb.Min.z, aabb.Max.z, ref tNear, ref tFar))
            return null;

        return tNear >= tMin ? tNear : tFar >= tMin ? tFar : (double?)null;
    }

    /// <summary>
    /// Tests whether a ray intersects a bounding sphere.
    /// </summary>
    /// <param name="ray">The optics ray to test.</param>
    /// <param name="sphere">The sphere to test against.</param>
    /// <param name="tMin">Minimum ray parameter.</param>
    /// <param name="tMax">Maximum ray parameter.</param>
    /// <returns>
    /// The distance along the ray to the nearest intersection,
    /// or null if no intersection occurs.
    /// </returns>
    public static double? IntersectSphere(this Ray ray, BoundingSphere sphere,
        double tMin = 1e-6, double tMax = double.MaxValue)
    {
        Vector oc = ray.Origin - sphere.Center;
        double a = ray.Direction.Dot(ray.Direction);
        double b = 2.0 * oc.Dot(ray.Direction);
        double c = oc.Dot(oc) - sphere.Radius * sphere.Radius;
        double disc = b * b - 4.0 * a * c;
        if (disc < 0) return null;

        double sqrtDisc = Math.Sqrt(disc);
        double t1 = (-b - sqrtDisc) / (2.0 * a);
        double t2 = (-b + sqrtDisc) / (2.0 * a);

        if (t1 >= tMin && t1 <= tMax) return t1;
        if (t2 >= tMin && t2 <= tMax) return t2;
        return null;
    }

    /// <summary>
    /// Returns the hit point and outward normal for a ray-sphere intersection.
    /// </summary>
    public static (Vector point, Vector normal)? IntersectSphereDetailed(this Ray ray,
        BoundingSphere sphere, double tMin = 1e-6, double tMax = double.MaxValue)
    {
        var t = ray.IntersectSphere(sphere, tMin, tMax);
        if (!t.HasValue) return null;

        Vector point = ray.PointAt(t.Value);
        Vector normal = (point - sphere.Center).GetUnitVector();
        return (point, normal);
    }

    private static bool SlabTest(double origin, double direction,
        double slabMin, double slabMax, ref double tNear, ref double tFar)
    {
        if (Math.Abs(direction) < 1e-15)
        {
            // Ray is parallel to slab — must be within the slab
            return origin >= slabMin && origin <= slabMax;
        }

        double invD = 1.0 / direction;
        double t0 = (slabMin - origin) * invD;
        double t1 = (slabMax - origin) * invD;
        if (t0 > t1) { double tmp = t0; t0 = t1; t1 = tmp; }

        tNear = t0 > tNear ? t0 : tNear;
        tFar = t1 < tFar ? t1 : tFar;

        return tNear <= tFar;
    }
}
