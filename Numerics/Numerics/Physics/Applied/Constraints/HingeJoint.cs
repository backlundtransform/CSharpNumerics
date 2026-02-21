using CSharpNumerics.Physics.Objects;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Applied.Constraints
{
    /// <summary>
    /// Constrains two bodies to share a pivot point and only rotate about a single axis.
    /// Combines a 3-DOF position constraint (ball-socket) with a 2-DOF angular constraint
    /// that locks the two rotation axes perpendicular to the hinge.
    /// Used for doors, elbows, wheels, and any single-axis revolute joint.
    /// </summary>
    public class HingeJoint : IConstraint
    {
        /// <summary>Index of body A in the bodies array.</summary>
        public int BodyA;

        /// <summary>Index of body B in the bodies array.</summary>
        public int BodyB;

        /// <summary>Anchor point in body A's local frame.</summary>
        public Vector LocalAnchorA;

        /// <summary>Anchor point in body B's local frame.</summary>
        public Vector LocalAnchorB;

        /// <summary>Hinge rotation axis in world space (unit vector).</summary>
        public Vector HingeAxis;

        /// <summary>Baumgarte stabilization factor for the position constraint. Default 0.2.</summary>
        public double BaumgarteFactor = 0.2;

        // Position constraint state (same as BallSocket)
        private Vector _rA, _rB, _positionBias;
        private double _effMassX, _effMassY, _effMassZ;

        // Angular constraint state (2 locked axes)
        private Vector _perp1, _perp2;
        private double _angEffMass1, _angEffMass2;

        public HingeJoint(int bodyA, int bodyB, Vector localAnchorA, Vector localAnchorB, Vector hingeAxis)
        {
            BodyA = bodyA;
            BodyB = bodyB;
            LocalAnchorA = localAnchorA;
            LocalAnchorB = localAnchorB;
            HingeAxis = hingeAxis.GetUnitVector();

            // Build perpendicular basis for the locked axes
            var temp = Math.Abs(HingeAxis.x) < 0.9 ? new Vector(1, 0, 0) : new Vector(0, 1, 0);
            _perp1 = HingeAxis.Cross(temp);
            _perp1 = (1.0 / _perp1.GetMagnitude()) * _perp1;
            _perp2 = HingeAxis.Cross(_perp1);
        }

        /// <summary>
        /// Creates a hinge joint at a world-space pivot with a given rotation axis.
        /// Local anchors are computed from the bodies' current positions.
        /// </summary>
        public static HingeJoint FromWorldPivot(int bodyA, int bodyB, Vector worldPivot, Vector hingeAxis, RigidBody[] bodies)
        {
            var localA = worldPivot - bodies[bodyA].Position;
            var localB = worldPivot - bodies[bodyB].Position;
            return new HingeJoint(bodyA, bodyB, localA, localB, hingeAxis);
        }

        public void PreStep(RigidBody[] bodies, double dt)
        {
            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            _rA = a.Orientation * LocalAnchorA;
            _rB = b.Orientation * LocalAnchorB;

            // --- Position constraint ---
            var worldA = a.Position + _rA;
            var worldB = b.Position + _rB;
            _positionBias = (BaumgarteFactor / dt) * (worldB - worldA);

            _effMassX = ComputeLinearEffectiveMass(a, b, _rA, _rB, new Vector(1, 0, 0));
            _effMassY = ComputeLinearEffectiveMass(a, b, _rA, _rB, new Vector(0, 1, 0));
            _effMassZ = ComputeLinearEffectiveMass(a, b, _rA, _rB, new Vector(0, 0, 1));

            // --- Angular constraint effective masses ---
            _angEffMass1 = ComputeAngularEffectiveMass(a, b, _perp1);
            _angEffMass2 = ComputeAngularEffectiveMass(a, b, _perp2);
        }

        public void Solve(RigidBody[] bodies)
        {
            ref var a = ref bodies[BodyA];
            ref var b = ref bodies[BodyB];

            // --- Position constraint (same as BallSocket) ---
            var vA = a.Velocity + a.AngularVelocity.Cross(_rA);
            var vB = b.Velocity + b.AngularVelocity.Cross(_rB);
            var vRel = vB - vA;

            var posImpulse = new Vector(
                -_effMassX * (vRel.x + _positionBias.x),
                -_effMassY * (vRel.y + _positionBias.y),
                -_effMassZ * (vRel.z + _positionBias.z));

            a.Velocity = a.Velocity - a.InverseMass * posImpulse;
            a.AngularVelocity = a.AngularVelocity - a.InverseInertiaTensorWorld * _rA.Cross(posImpulse);
            b.Velocity = b.Velocity + b.InverseMass * posImpulse;
            b.AngularVelocity = b.AngularVelocity + b.InverseInertiaTensorWorld * _rB.Cross(posImpulse);

            // --- Angular constraint: kill relative ω on the two locked axes ---
            var omegaRel = b.AngularVelocity - a.AngularVelocity;

            double angLambda1 = -_angEffMass1 * _perp1.Dot(omegaRel);
            double angLambda2 = -_angEffMass2 * _perp2.Dot(omegaRel);

            var angImpulse = angLambda1 * _perp1 + angLambda2 * _perp2;
            a.AngularVelocity = a.AngularVelocity - a.InverseInertiaTensorWorld * angImpulse;
            b.AngularVelocity = b.AngularVelocity + b.InverseInertiaTensorWorld * angImpulse;
        }

        private static double ComputeLinearEffectiveMass(RigidBody a, RigidBody b, Vector rA, Vector rB, Vector axis)
        {
            double K = a.InverseMass + b.InverseMass;
            var rAxN = rA.Cross(axis);
            var rBxN = rB.Cross(axis);
            K += (a.InverseInertiaTensorWorld * rAxN).Dot(rAxN);
            K += (b.InverseInertiaTensorWorld * rBxN).Dot(rBxN);
            return K > 1e-30 ? 1.0 / K : 0;
        }

        private static double ComputeAngularEffectiveMass(RigidBody a, RigidBody b, Vector axis)
        {
            double K = (a.InverseInertiaTensorWorld * axis).Dot(axis)
                     + (b.InverseInertiaTensorWorld * axis).Dot(axis);
            return K > 1e-30 ? 1.0 / K : 0;
        }
    }
}
