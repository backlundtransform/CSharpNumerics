using CSharpNumerics.Physics.Applied.Objects;
using Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Applied
{
    /// <summary>
    /// Provides overlap tests (broad phase) and contact generation (narrow phase)
    /// for bounding volumes and simple shapes.
    /// </summary>
    public static class CollisionDetection
    {
        #region Broad Phase — Overlap Tests

        /// <summary>
        /// Tests whether two AABBs overlap (inclusive of touching).
        /// </summary>
        public static bool Intersects(this AABB a, AABB b)
        {
            return a.Min.x <= b.Max.x && a.Max.x >= b.Min.x
                && a.Min.y <= b.Max.y && a.Max.y >= b.Min.y
                && a.Min.z <= b.Max.z && a.Max.z >= b.Min.z;
        }

        /// <summary>
        /// Tests whether two bounding spheres overlap.
        /// </summary>
        public static bool Intersects(this BoundingSphere a, BoundingSphere b)
        {
            var d = b.Center - a.Center;
            double distSq = d.Dot(d);
            double radiusSum = a.Radius + b.Radius;
            return distSq <= radiusSum * radiusSum;
        }

        /// <summary>
        /// Tests whether an AABB and a bounding sphere overlap.
        /// </summary>
        public static bool Intersects(this AABB box, BoundingSphere sphere)
        {
            var closest = box.ClosestPoint(sphere.Center);
            var d = sphere.Center - closest;
            return d.Dot(d) <= sphere.Radius * sphere.Radius;
        }

        /// <summary>
        /// Tests whether a bounding sphere and an AABB overlap.
        /// </summary>
        public static bool Intersects(this BoundingSphere sphere, AABB box)
        {
            return box.Intersects(sphere);
        }

        #endregion

        #region Narrow Phase — Contact Generation

        /// <summary>
        /// Computes the contact point between two intersecting spheres.
        /// Returns null if they are not overlapping.
        /// Normal points from A toward B.
        /// </summary>
        public static ContactPoint? SphereSphereContact(this BoundingSphere a, BoundingSphere b)
        {
            var delta = b.Center - a.Center;
            double distSq = delta.Dot(delta);
            double radiusSum = a.Radius + b.Radius;

            if (distSq > radiusSum * radiusSum)
                return null;

            double dist = Math.Sqrt(distSq);
            Vector normal;

            if (dist < 1e-15)
            {
                // Coincident centers — pick arbitrary separation axis
                normal = new Vector(1, 0, 0);
                dist = 0;
            }
            else
            {
                normal = (1.0 / dist) * delta;
            }

            double penetration = radiusSum - dist;
            var position = a.Center + a.Radius * normal;

            return new ContactPoint(position, normal, penetration);
        }

        /// <summary>
        /// Computes the contact point between a sphere and an AABB.
        /// Returns null if they are not overlapping.
        /// Normal points from box toward sphere.
        /// </summary>
        public static ContactPoint? SphereAABBContact(this BoundingSphere sphere, AABB box)
        {
            var closest = box.ClosestPoint(sphere.Center);
            var delta = sphere.Center - closest;
            double distSq = delta.Dot(delta);

            if (distSq > sphere.Radius * sphere.Radius)
                return null;

            double dist = Math.Sqrt(distSq);
            Vector normal;

            if (dist < 1e-15)
            {
                // Sphere center is inside the box — find the shallowest exit axis
                normal = ShallowExitNormal(sphere.Center, box);
                double penetration = sphere.Radius + ShallowExitDepth(sphere.Center, box);
                return new ContactPoint(closest, normal, penetration);
            }
            else
            {
                normal = (1.0 / dist) * delta;
            }

            double pen = sphere.Radius - dist;
            return new ContactPoint(closest, normal, pen);
        }

        /// <summary>
        /// Computes the contact point between an AABB and a sphere.
        /// Normal points from sphere toward box (reversed from SphereAABBContact).
        /// </summary>
        public static ContactPoint? AABBSphereContact(this AABB box, BoundingSphere sphere)
        {
            var contact = sphere.SphereAABBContact(box);
            if (contact is null) return null;

            var c = contact.Value;
            return new ContactPoint(c.Position, -1.0 * c.Normal, c.PenetrationDepth);
        }

        #endregion

        #region Helpers

        private static Vector ShallowExitNormal(Vector point, AABB box)
        {
            double dx = Math.Min(point.x - box.Min.x, box.Max.x - point.x);
            double dy = Math.Min(point.y - box.Min.y, box.Max.y - point.y);
            double dz = Math.Min(point.z - box.Min.z, box.Max.z - point.z);

            if (dx <= dy && dx <= dz)
                return new Vector(point.x - box.Min.x < box.Max.x - point.x ? -1 : 1, 0, 0);
            if (dy <= dz)
                return new Vector(0, point.y - box.Min.y < box.Max.y - point.y ? -1 : 1, 0);
            return new Vector(0, 0, point.z - box.Min.z < box.Max.z - point.z ? -1 : 1);
        }

        private static double ShallowExitDepth(Vector point, AABB box)
        {
            double dx = Math.Min(point.x - box.Min.x, box.Max.x - point.x);
            double dy = Math.Min(point.y - box.Min.y, box.Max.y - point.y);
            double dz = Math.Min(point.z - box.Min.z, box.Max.z - point.z);
            return Math.Min(dx, Math.Min(dy, dz));
        }

        #endregion
    }
}
