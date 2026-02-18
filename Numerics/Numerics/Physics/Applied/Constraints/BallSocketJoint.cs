using CSharpNumerics.Physics.Objects;
using Numerics.Objects;

namespace CSharpNumerics.Physics.Applied.Constraints
{
    /// <summary>
    /// Constrains two anchor points to be coincident (zero distance).
    /// C = worldAnchorB - worldAnchorA = 0 (3 equations).
    /// Allows free rotation — the bodies can spin freely around the joint point.
    /// Used for shoulders, ragdoll joints, pendulums, and chain links.
    /// </summary>
    public class BallSocketJoint : IConstraint
    {
        /// <summary>Index of body A in the bodies array.</summary>
        public int BodyA;

        /// <summary>Index of body B in the bodies array.</summary>
        public int BodyB;

        /// <summary>Anchor point in body A's local frame.</summary>
        public Vector LocalAnchorA;

        /// <summary>Anchor point in body B's local frame.</summary>
        public Vector LocalAnchorB;

        /// <summary>Baumgarte stabilization factor. Default 0.2.</summary>
        public double BaumgarteFactor = 0.2;

        private Vector _rA, _rB, _bias;
        private double _effMassX, _effMassY, _effMassZ;

        public BallSocketJoint(int bodyA, int bodyB, Vector localAnchorA, Vector localAnchorB)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = localAnchorA;
            LocalAnchorB = localAnchorB;
        }

        /// <summary>
        /// Creates a ball-socket joint at a world-space pivot point.
        /// Local anchors are computed from the bodies' current positions.
        /// </summary>
        public static BallSocketJoint FromWorldPivot(int bodyA, int bodyB, Vector worldPivot, RigidBody[] bodies)
        {
            var localA = worldPivot - bodies[bodyA].Position;
            var localB = worldPivot - bodies[bodyB].Position;
            return new BallSocketJoint(bodyA, bodyB, localA, localB);
        }

        public void PreStep(RigidBody[] bodies, double dt)
        {
            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            _rA = a.Orientation * LocalAnchorA;
            _rB = b.Orientation * LocalAnchorB;

            var worldA = a.Position + _rA;
            var worldB = b.Position + _rB;
            _bias = (BaumgarteFactor / dt) * (worldB - worldA);

            _effMassX = ComputeEffectiveMass(a, b, _rA, _rB, new Vector(1, 0, 0));
            _effMassY = ComputeEffectiveMass(a, b, _rA, _rB, new Vector(0, 1, 0));
            _effMassZ = ComputeEffectiveMass(a, b, _rA, _rB, new Vector(0, 0, 1));
        }

        public void Solve(RigidBody[] bodies)
        {
            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            var vA = a.Velocity + a.AngularVelocity.Cross(_rA);
            var vB = b.Velocity + b.AngularVelocity.Cross(_rB);
            var vRel = vB - vA;

            var impulse = new Vector(
                -_effMassX * (vRel.x + _bias.x),
                -_effMassY * (vRel.y + _bias.y),
                -_effMassZ * (vRel.z + _bias.z));

            a.Velocity = a.Velocity - a.InverseMass * impulse;
            a.AngularVelocity = a.AngularVelocity - a.InverseInertiaTensorWorld * _rA.Cross(impulse);
            b.Velocity = b.Velocity + b.InverseMass * impulse;
            b.AngularVelocity = b.AngularVelocity + b.InverseInertiaTensorWorld * _rB.Cross(impulse);
        }

        private static double ComputeEffectiveMass(RigidBody a, RigidBody b, Vector rA, Vector rB, Vector axis)
        {
            double K = a.InverseMass + b.InverseMass;
            var rAxN = rA.Cross(axis);
            var rBxN = rB.Cross(axis);
            K += (a.InverseInertiaTensorWorld * rAxN).Dot(rAxN);
            K += (b.InverseInertiaTensorWorld * rBxN).Dot(rBxN);
            return K > 1e-30 ? 1.0 / K : 0;
        }
    }
}
