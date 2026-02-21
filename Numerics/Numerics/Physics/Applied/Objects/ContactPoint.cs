using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Applied.Objects
{
    /// <summary>
    /// Describes a contact point between two colliding bodies.
    /// Used by the collision response system to compute impulses.
    /// </summary>
    public struct ContactPoint
    {
        /// <summary>Contact position in world space.</summary>
        public Vector Position;

        /// <summary>
        /// Contact normal pointing from body A toward body B.
        /// Must be a unit vector.
        /// </summary>
        public Vector Normal;

        /// <summary>
        /// Overlap depth (positive = penetrating, zero = touching).
        /// </summary>
        public double PenetrationDepth;

        public ContactPoint(Vector position, Vector normal, double penetrationDepth)
        {
            Position = position;
            Normal = normal;
            PenetrationDepth = penetrationDepth;
        }

        public override string ToString()
            => $"Contact(Pos={Position}, Normal={Normal}, Depth={PenetrationDepth:F4})";
    }
}
