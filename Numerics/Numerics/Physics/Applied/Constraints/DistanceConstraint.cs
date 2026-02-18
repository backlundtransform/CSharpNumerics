using CSharpNumerics.Physics.Objects;
using Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Applied.Constraints
{
    /// <summary>
    /// Maintains a fixed distance between two anchor points on two bodies.
    /// C = |worldB - worldA| - L = 0.
    /// Equivalent to a rigid rod connecting the two points.
    /// </summary>
    public class DistanceConstraint : IConstraint
    {
        /// <summary>Index of body A in the bodies array.</summary>
        public int BodyA;

        /// <summary>Index of body B in the bodies array.</summary>
        public int BodyB;

        /// <summary>Anchor point in body A's local frame (offset from center of mass).</summary>
        public Vector LocalAnchorA;

        /// <summary>Anchor point in body B's local frame (offset from center of mass).</summary>
        public Vector LocalAnchorB;

        /// <summary>Target distance between the two anchor points in meters.</summary>
        public double Distance;

        /// <summary>Baumgarte stabilization factor. Higher = faster drift correction but stiffer. Default 0.2.</summary>
        public double BaumgarteFactor = 0.2;

        private Vector _n, _rA, _rB;
        private double _effectiveMass, _bias;

        public DistanceConstraint(int bodyA, int bodyB, Vector localAnchorA, Vector localAnchorB, double distance)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = localAnchorA;
            LocalAnchorB = localAnchorB;
            Distance = distance;
        }

        /// <summary>
        /// Creates a distance constraint using the current positions of two bodies.
        /// The distance is computed automatically from their current separation.
        /// </summary>
        public static DistanceConstraint FromCurrentPositions(int bodyA, int bodyB, RigidBody[] bodies)
        {
            var d = bodies[bodyB].Position - bodies[bodyA].Position;
            return new DistanceConstraint(bodyA, bodyB, new Vector(0, 0, 0), new Vector(0, 0, 0), d.GetMagnitude());
        }

        public void PreStep(RigidBody[] bodies, double dt)
        {
            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            _rA = a.Orientation * LocalAnchorA;
            _rB = b.Orientation * LocalAnchorB;

            var worldA = a.Position + _rA;
            var worldB = b.Position + _rB;
            var d = worldB - worldA;
            double dist = d.GetMagnitude();

            _n = dist > 1e-15 ? (1.0 / dist) * d : new Vector(1, 0, 0);

            // Effective mass: K = 1/(invMa + invMb + angular terms)
            double K = a.InverseMass + b.InverseMass;
            var rAxN = _rA.Cross(_n);
            var rBxN = _rB.Cross(_n);
            K += (a.InverseInertiaTensorWorld * rAxN).Dot(rAxN);
            K += (b.InverseInertiaTensorWorld * rBxN).Dot(rBxN);
            _effectiveMass = K > 1e-30 ? 1.0 / K : 0;

            // Baumgarte bias: push velocity toward resolving positional error
            _bias = BaumgarteFactor / dt * (dist - Distance);
        }

        public void Solve(RigidBody[] bodies)
        {
            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            var vA = a.Velocity + a.AngularVelocity.Cross(_rA);
            var vB = b.Velocity + b.AngularVelocity.Cross(_rB);
            double Jv = _n.Dot(vB - vA);

            double lambda = -_effectiveMass * (Jv + _bias);

            var impulse = lambda * _n;
            a.Velocity = a.Velocity - a.InverseMass * impulse;
            a.AngularVelocity = a.AngularVelocity - a.InverseInertiaTensorWorld * _rA.Cross(impulse);
            b.Velocity = b.Velocity + b.InverseMass * impulse;
            b.AngularVelocity = b.AngularVelocity + b.InverseInertiaTensorWorld * _rB.Cross(impulse);
        }
    }
}
