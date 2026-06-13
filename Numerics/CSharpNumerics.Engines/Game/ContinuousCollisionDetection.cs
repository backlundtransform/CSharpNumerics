using CSharpNumerics.Engines.Game.Objects;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game;

/// <summary>
/// Continuous Collision Detection (CCD) for fast-moving objects.
///
/// Prevents tunnelling by sweeping a sphere along its trajectory
/// and finding the earliest intersection with other geometry.
/// </summary>
public static class ContinuousCollisionDetection
{
    /// <summary>
    /// Result of a CCD sweep test.
    /// </summary>
    public struct SweepResult
    {
        /// <summary>Whether a collision was detected.</summary>
        public bool Hit;

        /// <summary>Parametric time of impact [0,1] along the sweep.</summary>
        public double TimeOfImpact;

        /// <summary>Position of the sphere center at impact.</summary>
        public Vector HitPosition;

        /// <summary>Contact normal at impact point.</summary>
        public Vector HitNormal;
    }

    /// <summary>
    /// Swept sphere vs. static sphere test.
    /// Tests if a sphere moving from startPos to endPos (with radius rA) 
    /// intersects a static sphere (center, rB).
    /// </summary>
    /// <param name="startPos">Start position of moving sphere.</param>
    /// <param name="endPos">End position of moving sphere.</param>
    /// <param name="radiusA">Radius of moving sphere.</param>
    /// <param name="center">Center of static sphere.</param>
    /// <param name="radiusB">Radius of static sphere.</param>
    /// <returns>Sweep result with time of impact.</returns>
    public static SweepResult SweptSphereVsSphere(
        Vector startPos, Vector endPos, double radiusA,
        Vector center, double radiusB)
    {
        var result = new SweepResult();
        double rSum = radiusA + radiusB;

        // Ray: P(t) = start + t * d,  t ∈ [0,1]
        var d = endPos - startPos;
        var m = startPos - center;

        double a = d.Dot(d);
        double b = 2.0 * m.Dot(d);
        double c = m.Dot(m) - rSum * rSum;

        // Already overlapping at start
        if (c <= 0)
        {
            result.Hit = true;
            result.TimeOfImpact = 0;
            result.HitPosition = startPos;
            double dist = m.GetMagnitude();
            result.HitNormal = dist > 1e-10 ? (1.0 / dist) * m : new Vector(0, 0, 1);
            return result;
        }

        double discriminant = b * b - 4.0 * a * c;
        if (discriminant < 0 || a < 1e-20)
            return result; // No intersection

        double sqrtD = Math.Sqrt(discriminant);
        double t = (-b - sqrtD) / (2.0 * a);

        if (t >= 0 && t <= 1.0)
        {
            result.Hit = true;
            result.TimeOfImpact = t;
            result.HitPosition = startPos + t * d;
            var normal = result.HitPosition - center;
            double mag = normal.GetMagnitude();
            result.HitNormal = mag > 1e-10 ? (1.0 / mag) * normal : new Vector(0, 0, 1);
        }

        return result;
    }

    /// <summary>
    /// Swept sphere vs. AABB test.
    /// Uses Minkowski sum: expand AABB by sphere radius and ray-cast against expanded box.
    /// </summary>
    public static SweepResult SweptSphereVsAABB(
        Vector startPos, Vector endPos, double radius, AABB box)
    {
        var result = new SweepResult();

        // Expand AABB by sphere radius
        var expandedMin = new Vector(box.Min.x - radius, box.Min.y - radius, box.Min.z - radius);
        var expandedMax = new Vector(box.Max.x + radius, box.Max.y + radius, box.Max.z + radius);

        var d = endPos - startPos;

        double tMin = 0, tMax = 1;
        int hitAxis = -1;
        bool hitLow = false;

        // Slab test on each axis
        if (!SlabTest(startPos.x, d.x, expandedMin.x, expandedMax.x, ref tMin, ref tMax, 0, ref hitAxis, ref hitLow))
            return result;
        if (!SlabTest(startPos.y, d.y, expandedMin.y, expandedMax.y, ref tMin, ref tMax, 1, ref hitAxis, ref hitLow))
            return result;
        if (!SlabTest(startPos.z, d.z, expandedMin.z, expandedMax.z, ref tMin, ref tMax, 2, ref hitAxis, ref hitLow))
            return result;

        if (tMin >= 0 && tMin <= 1)
        {
            result.Hit = true;
            result.TimeOfImpact = tMin;
            result.HitPosition = startPos + tMin * d;

            // Normal from the hit axis
            result.HitNormal = hitAxis switch
            {
                0 => new Vector(hitLow ? -1 : 1, 0, 0),
                1 => new Vector(0, hitLow ? -1 : 1, 0),
                _ => new Vector(0, 0, hitLow ? -1 : 1),
            };
        }

        return result;
    }

    /// <summary>
    /// Swept sphere vs. ground plane at z = height.
    /// </summary>
    public static SweepResult SweptSphereVsPlane(
        Vector startPos, Vector endPos, double radius, double planeHeight = 0)
    {
        var result = new SweepResult();

        double startBottom = startPos.z - radius;
        double endBottom = endPos.z - radius;

        if (startBottom >= planeHeight && endBottom >= planeHeight)
            return result; // No collision

        if (startBottom < planeHeight)
        {
            // Already penetrating
            result.Hit = true;
            result.TimeOfImpact = 0;
            result.HitPosition = startPos;
            result.HitNormal = new Vector(0, 0, 1);
            return result;
        }

        // Find t where z - radius = planeHeight
        double dz = endPos.z - startPos.z;
        if (Math.Abs(dz) < 1e-15) return result;

        double t = (startBottom - planeHeight) / (startBottom - endBottom);
        if (t >= 0 && t <= 1)
        {
            var d = endPos - startPos;
            result.Hit = true;
            result.TimeOfImpact = t;
            result.HitPosition = startPos + t * d;
            result.HitNormal = new Vector(0, 0, 1);
        }

        return result;
    }

    private static bool SlabTest(double origin, double dir, double lo, double hi,
        ref double tMin, ref double tMax, int axis, ref int hitAxis, ref bool hitLow)
    {
        if (Math.Abs(dir) < 1e-15)
        {
            return origin >= lo && origin <= hi;
        }

        double invD = 1.0 / dir;
        double t1 = (lo - origin) * invD;
        double t2 = (hi - origin) * invD;
        bool enteredLow = true;
        if (t1 > t2)
        {
            (t1, t2) = (t2, t1);
            enteredLow = false;
        }

        if (t1 > tMin)
        {
            tMin = t1;
            hitAxis = axis;
            hitLow = enteredLow;
        }
        tMax = Math.Min(tMax, t2);

        return tMax >= tMin;
    }
}
