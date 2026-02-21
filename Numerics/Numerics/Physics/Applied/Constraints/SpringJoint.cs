using CSharpNumerics.Physics.Objects;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Physics.Applied.Constraints
{
    /// <summary>
    /// Soft spring-damper constraint between two anchor points.
    /// Applies F = -k·(|d| - L₀) - c·v_rel as an impulse once per frame.
    /// Unlike hard constraints, this does not enforce an exact relationship —
    /// it applies a restoring force, so the bodies oscillate around equilibrium.
    /// </summary>
    public class SpringJoint : IConstraint
    {
        /// <summary>Index of body A in the bodies array.</summary>
        public int BodyA;

        /// <summary>Index of body B in the bodies array.</summary>
        public int BodyB;

        /// <summary>Anchor point in body A's local frame.</summary>
        public Vector LocalAnchorA;

        /// <summary>Anchor point in body B's local frame.</summary>
        public Vector LocalAnchorB;

        /// <summary>Spring stiffness in N/m.</summary>
        public double Stiffness;

        /// <summary>Damping coefficient in N·s/m.</summary>
        public double Damping;

        /// <summary>Rest length of the spring in meters.</summary>
        public double RestLength;

        private double _dt;
        private bool _applied;

        public SpringJoint(int bodyA, int bodyB, Vector localAnchorA, Vector localAnchorB,
            double stiffness, double damping, double restLength)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = localAnchorA;
            LocalAnchorB = localAnchorB;
            Stiffness = stiffness;
            Damping = damping;
            RestLength = restLength;
        }

        public void PreStep(RigidBody[] bodies, double dt)
        {
            _dt = dt;
            _applied = false;
        }

        public void Solve(RigidBody[] bodies)
        {
            // Soft constraint — apply once, not iterated
            if (_applied) return;
            _applied = true;

            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            var rA = a.Orientation * LocalAnchorA;
            var rB = b.Orientation * LocalAnchorB;
            var worldA = a.Position + rA;
            var worldB = b.Position + rB;

            var d = worldB - worldA;
            double dist = d.GetMagnitude();
            if (dist < 1e-15) return;

            var n = (1.0 / dist) * d;
            double stretch = dist - RestLength;

            // Relative velocity along spring axis
            var vA = a.Velocity + a.AngularVelocity.Cross(rA);
            var vB = b.Velocity + b.AngularVelocity.Cross(rB);
            double vRel = n.Dot(vB - vA);

            // Spring-damper impulse: J = (k·stretch + c·v_rel) · dt
            double j = (Stiffness * stretch + Damping * vRel) * _dt;

            // Impulse pulls A toward B (+n) and B toward A (-n)
            var impulse = j * n;
            a.Velocity = a.Velocity + a.InverseMass * impulse;
            a.AngularVelocity = a.AngularVelocity + a.InverseInertiaTensorWorld * rA.Cross(impulse);
            b.Velocity = b.Velocity - b.InverseMass * impulse;
            b.AngularVelocity = b.AngularVelocity - b.InverseInertiaTensorWorld * rB.Cross(impulse);
        }
    }
}
