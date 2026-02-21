using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Applied.Objects
{
    /// <summary>
    /// Axis-Aligned Bounding Box defined by minimum and maximum corners.
    /// Used for fast broad-phase overlap tests.
    /// </summary>
    public struct AABB
    {
        /// <summary>Minimum corner (lowest x, y, z).</summary>
        public Vector Min;

        /// <summary>Maximum corner (highest x, y, z).</summary>
        public Vector Max;

        public AABB(Vector min, Vector max)
        {
            Min = min;
            Max = max;
        }

        /// <summary>Center of the box: (Min + Max) / 2.</summary>
        public Vector Center => 0.5 * (Min + Max);

        /// <summary>Half-extents from center to each face.</summary>
        public Vector HalfExtents => 0.5 * (Max - Min);

        /// <summary>Full size along each axis.</summary>
        public Vector Size => Max - Min;

        /// <summary>Volume of the box.</summary>
        public double Volume
        {
            get
            {
                var s = Size;
                return s.x * s.y * s.z;
            }
        }

        /// <summary>
        /// Creates an AABB from a center point and half-extents.
        /// </summary>
        public static AABB FromCenterExtents(Vector center, Vector halfExtents)
        {
            return new AABB(center - halfExtents, center + halfExtents);
        }

        /// <summary>
        /// Returns true if the point is inside or on the surface of this box.
        /// </summary>
        public bool Contains(Vector point)
        {
            return point.x >= Min.x && point.x <= Max.x
                && point.y >= Min.y && point.y <= Max.y
                && point.z >= Min.z && point.z <= Max.z;
        }

        /// <summary>
        /// Returns the closest point on or inside this AABB to the given point.
        /// </summary>
        public Vector ClosestPoint(Vector point)
        {
            return new Vector(
                Math.Clamp(point.x, Min.x, Max.x),
                Math.Clamp(point.y, Min.y, Max.y),
                Math.Clamp(point.z, Min.z, Max.z));
        }

        /// <summary>
        /// Returns the smallest AABB that contains both this and another AABB.
        /// </summary>
        public AABB Merge(AABB other)
        {
            return new AABB(
                new Vector(Math.Min(Min.x, other.Min.x), Math.Min(Min.y, other.Min.y), Math.Min(Min.z, other.Min.z)),
                new Vector(Math.Max(Max.x, other.Max.x), Math.Max(Max.y, other.Max.y), Math.Max(Max.z, other.Max.z)));
        }

        /// <summary>
        /// Expands this AABB uniformly by margin in all directions.
        /// </summary>
        public AABB Expand(double margin)
        {
            var m = new Vector(margin, margin, margin);
            return new AABB(Min - m, Max + m);
        }

        public override string ToString() => $"AABB(Min={Min}, Max={Max})";
    }
}
