using Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Applied.Objects
{
    /// <summary>
    /// Bounding sphere defined by center and radius.
    /// Used for fast broad-phase overlap tests.
    /// </summary>
    public struct BoundingSphere
    {
        /// <summary>Center of the sphere in world space.</summary>
        public Vector Center;

        /// <summary>Radius of the sphere.</summary>
        public double Radius;

        public BoundingSphere(Vector center, double radius)
        {
            Center = center;
            Radius = Math.Max(0, radius);
        }

        /// <summary>Volume: 4/3·π·r³.</summary>
        public double Volume => (4.0 / 3.0) * Math.PI * Radius * Radius * Radius;

        /// <summary>Surface area: 4·π·r².</summary>
        public double SurfaceArea => 4.0 * Math.PI * Radius * Radius;

        /// <summary>
        /// Returns true if the point is inside or on the surface of this sphere.
        /// </summary>
        public bool Contains(Vector point)
        {
            var d = point - Center;
            return d.Dot(d) <= Radius * Radius;
        }

        /// <summary>
        /// Returns the closest point on the surface of this sphere to the given point.
        /// If the point is at the center, returns a point along +X.
        /// </summary>
        public Vector ClosestPoint(Vector point)
        {
            var d = point - Center;
            double dist = d.GetMagnitude();
            if (dist < 1e-15)
                return Center + new Vector(Radius, 0, 0);
            return Center + (Radius / dist) * d;
        }

        /// <summary>
        /// Returns the smallest sphere that contains both this and another sphere.
        /// </summary>
        public BoundingSphere Merge(BoundingSphere other)
        {
            var delta = other.Center - Center;
            double dist = delta.GetMagnitude();

            // One sphere contains the other
            if (dist + other.Radius <= Radius) return this;
            if (dist + Radius <= other.Radius) return other;

            double newRadius = (dist + Radius + other.Radius) * 0.5;
            var direction = dist > 1e-15 ? (1.0 / dist) * delta : new Vector(0, 0, 0);
            var newCenter = Center + (newRadius - Radius) * direction;
            return new BoundingSphere(newCenter, newRadius);
        }

        public override string ToString() => $"Sphere(Center={Center}, R={Radius:F3})";
    }
}
