using CSharpNumerics.Physics.Applied.Objects;
using CSharpNumerics.Physics.Objects;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Applied
{
    /// <summary>
    /// Impulse-based collision response for rigid bodies.
    /// Resolves velocity (impulse) and position (penetration correction) at contact points.
    /// </summary>
    public static class CollisionResponse
    {
        /// <summary>
        /// Resolves a collision between two rigid bodies at a contact point using impulse-based response.
        /// Updates velocities and angular velocities of both bodies.
        /// Applies both normal impulse (bounce) and tangential friction impulse.
        /// </summary>
        /// <param name="a">First rigid body (normal points away from A).</param>
        /// <param name="b">Second rigid body (normal points toward B).</param>
        /// <param name="contact">Contact information (normal from A toward B).</param>
        /// <param name="restitution">Coefficient of restitution: 0 = perfectly inelastic, 1 = elastic.</param>
        /// <param name="friction">Coefficient of friction for tangential impulse (default 0).</param>
        public static void ResolveCollision(
            ref RigidBody a,
            ref RigidBody b,
            ContactPoint contact,
            double restitution,
            double friction = 0.0)
        {
            if (a.IsStatic && b.IsStatic) return;

            var n = contact.Normal;

            // Vectors from center of mass to contact point
            var rA = contact.Position - a.Position;
            var rB = contact.Position - b.Position;

            // Relative velocity at contact point
            var vA = a.Velocity + a.AngularVelocity.Cross(rA);
            var vB = b.Velocity + b.AngularVelocity.Cross(rB);
            var vRel = vB - vA;

            // Relative velocity along normal
            double vn = vRel.Dot(n);

            // Don't resolve if bodies are separating
            if (vn > 0) return;

            // Effective inverse mass at contact
            double invMassSum = a.InverseMass + b.InverseMass;

            var rAxN = rA.Cross(n);
            var rBxN = rB.Cross(n);
            double angularTermA = (a.InverseInertiaTensorWorld * rAxN).Dot(rAxN);
            double angularTermB = (b.InverseInertiaTensorWorld * rBxN).Dot(rBxN);

            double denominator = invMassSum + angularTermA + angularTermB;
            if (denominator < 1e-30) return;

            // Normal impulse magnitude
            double jn = -(1.0 + restitution) * vn / denominator;

            // Apply normal impulse
            var impulse = jn * n;

            a.Velocity = a.Velocity - a.InverseMass * impulse;
            b.Velocity = b.Velocity + b.InverseMass * impulse;
            a.AngularVelocity = a.AngularVelocity - a.InverseInertiaTensorWorld * rA.Cross(impulse);
            b.AngularVelocity = b.AngularVelocity + b.InverseInertiaTensorWorld * rB.Cross(impulse);

            // Friction impulse (Coulomb model)
            if (friction > 0)
            {
                // Recompute relative velocity after normal impulse
                vA = a.Velocity + a.AngularVelocity.Cross(rA);
                vB = b.Velocity + b.AngularVelocity.Cross(rB);
                vRel = vB - vA;

                // Tangent direction
                var tangent = vRel - vRel.Dot(n) * n;
                double tangentSpeed = tangent.GetMagnitude();
                if (tangentSpeed > 1e-15)
                {
                    tangent = (1.0 / tangentSpeed) * tangent;

                    // Friction denominator
                    var rAxT = rA.Cross(tangent);
                    var rBxT = rB.Cross(tangent);
                    double angTanA = (a.InverseInertiaTensorWorld * rAxT).Dot(rAxT);
                    double angTanB = (b.InverseInertiaTensorWorld * rBxT).Dot(rBxT);
                    double denomT = invMassSum + angTanA + angTanB;

                    double jt = -tangentSpeed / denomT;

                    // Clamp to Coulomb cone: |jt| ≤ μ·jn
                    jt = Math.Clamp(jt, -friction * jn, friction * jn);

                    var frictionImpulse = jt * tangent;

                    a.Velocity = a.Velocity - a.InverseMass * frictionImpulse;
                    b.Velocity = b.Velocity + b.InverseMass * frictionImpulse;
                    a.AngularVelocity = a.AngularVelocity - a.InverseInertiaTensorWorld * rA.Cross(frictionImpulse);
                    b.AngularVelocity = b.AngularVelocity + b.InverseInertiaTensorWorld * rB.Cross(frictionImpulse);
                }
            }
        }

        /// <summary>
        /// Applies positional correction to prevent bodies from sinking into each other (Baumgarte stabilization).
        /// Call after ResolveCollision each frame.
        /// </summary>
        /// <param name="a">First rigid body.</param>
        /// <param name="b">Second rigid body.</param>
        /// <param name="contact">Contact information.</param>
        /// <param name="correctionFraction">How much penetration to correct per step (0.2–0.8 typical). Default 0.4.</param>
        /// <param name="slop">Penetration depth below which no correction is applied (avoids jitter). Default 0.01.</param>
        public static void CorrectPositions(
            ref RigidBody a,
            ref RigidBody b,
            ContactPoint contact,
            double correctionFraction = 0.4,
            double slop = 0.01)
        {
            if (a.IsStatic && b.IsStatic) return;

            double penetration = contact.PenetrationDepth;
            if (penetration <= slop) return;

            double invMassSum = a.InverseMass + b.InverseMass;
            if (invMassSum < 1e-30) return;

            var correction = (correctionFraction * (penetration - slop) / invMassSum) * contact.Normal;

            a.Position = a.Position - a.InverseMass * correction;
            b.Position = b.Position + b.InverseMass * correction;
        }
    }
}
