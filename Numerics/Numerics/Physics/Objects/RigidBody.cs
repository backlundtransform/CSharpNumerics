using Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Objects
{
    /// <summary>
    /// Represents a rigid body with mass, inertia, position, velocity, and orientation.
    /// </summary>
    public struct RigidBody
    {
        // Material properties
        public double Mass;
        public double InverseMass;
        public Matrix InertiaTensorBody;
        public Matrix InverseInertiaTensorBody;

        // Dynamic state
        public Vector Position;
        public Vector Velocity;
        public Vector AngularVelocity;
        public Matrix Orientation;

        // Force/torque accumulators
        public Vector AccumulatedForce;
        public Vector AccumulatedTorque;

        private static readonly Matrix Identity3 = new Matrix(new double[,]
        {
            { 1, 0, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 }
        });

        /// <summary>
        /// Creates a rigid body with given mass and inertia tensor (in body frame).
        /// </summary>
        /// <param name="mass">Mass in kg. Use 0 for a static (immovable) body.</param>
        /// <param name="inertiaTensor">3×3 inertia tensor in body frame.</param>
        public RigidBody(double mass, Matrix inertiaTensor)
        {
            Mass = mass;
            InverseMass = mass > 0 ? 1.0 / mass : 0;
            InertiaTensorBody = inertiaTensor;
            InverseInertiaTensorBody = mass > 0 ? inertiaTensor.Inverse() : new Matrix(3, 3);

            Position = new Vector(0, 0, 0);
            Velocity = new Vector(0, 0, 0);
            AngularVelocity = new Vector(0, 0, 0);
            Orientation = new Matrix(Identity3.values);

            AccumulatedForce = new Vector(0, 0, 0);
            AccumulatedTorque = new Vector(0, 0, 0);
        }

        /// <summary>
        /// Creates a rigid body with given mass, inertia tensor, initial position and velocity.
        /// </summary>
        public RigidBody(double mass, Matrix inertiaTensor, Vector position, Vector velocity)
            : this(mass, inertiaTensor)
        {
            Position = position;
            Velocity = velocity;
        }

        /// <summary>
        /// Returns true if this body is static (infinite mass, does not move).
        /// </summary>
        public bool IsStatic => Mass <= 0;

        /// <summary>
        /// Computes the inverse inertia tensor in world frame: R · I⁻¹_body · Rᵀ.
        /// </summary>
        public Matrix InverseInertiaTensorWorld
        {
            get
            {
                var Rt = Orientation.Transpose();
                return Orientation * InverseInertiaTensorBody * Rt;
            }
        }

        /// <summary>
        /// Linear acceleration from accumulated forces: a = F/m.
        /// </summary>
        public Vector LinearAcceleration => InverseMass * AccumulatedForce;

        /// <summary>
        /// Angular acceleration from accumulated torques: α = I⁻¹_world · τ.
        /// </summary>
        public Vector AngularAcceleration => InverseInertiaTensorWorld * AccumulatedTorque;

        /// <summary>
        /// Linear momentum: p = mv.
        /// </summary>
        public Vector LinearMomentum => Mass * Velocity;

        /// <summary>
        /// Angular momentum in world frame: L = I_world · ω.
        /// </summary>
        public Vector AngularMomentum
        {
            get
            {
                var Rt = Orientation.Transpose();
                var Iworld = Orientation * InertiaTensorBody * Rt;
                return Iworld * AngularVelocity;
            }
        }

        /// <summary>
        /// Total kinetic energy (translational + rotational):
        /// KE = ½mv² + ½ωᵀIω.
        /// </summary>
        public double KineticEnergy
        {
            get
            {
                double translational = 0.5 * Mass * Velocity.Dot(Velocity);
                var Iw = InertiaTensorBody * AngularVelocity;
                double rotational = 0.5 * AngularVelocity.Dot(Iw);
                return translational + rotational;
            }
        }

        /// <summary>
        /// Applies a force at the center of mass (no torque generated).
        /// </summary>
        public void ApplyForce(Vector force)
        {
            AccumulatedForce = AccumulatedForce + force;
        }

        /// <summary>
        /// Applies a force at a world-space point, generating both force and torque.
        /// τ = (worldPoint - position) × force.
        /// </summary>
        public void ApplyForceAtPoint(Vector force, Vector worldPoint)
        {
            AccumulatedForce = AccumulatedForce + force;
            var r = worldPoint - Position;
            AccumulatedTorque = AccumulatedTorque + r.Cross(force);
        }

        /// <summary>
        /// Applies a pure torque (no linear force).
        /// </summary>
        public void ApplyTorque(Vector torque)
        {
            AccumulatedTorque = AccumulatedTorque + torque;
        }

        /// <summary>
        /// Resets accumulated force and torque to zero (call after each integration step).
        /// </summary>
        public void ClearForces()
        {
            AccumulatedForce = new Vector(0, 0, 0);
            AccumulatedTorque = new Vector(0, 0, 0);
        }

        #region Static Factory Methods

        /// <summary>
        /// Creates a static (immovable) rigid body at the given position.
        /// </summary>
        public static RigidBody CreateStatic(Vector position)
        {
            return new RigidBody(0, new Matrix(3, 3)) { Position = position };
        }

        /// <summary>
        /// Creates a solid sphere rigid body.
        /// I = diag(2/5·mr², 2/5·mr², 2/5·mr²).
        /// </summary>
        public static RigidBody CreateSolidSphere(double mass, double radius)
        {
            double I = 2.0 / 5.0 * mass * radius * radius;
            return new RigidBody(mass, DiagonalInertia(I, I, I));
        }

        /// <summary>
        /// Creates a hollow sphere (shell) rigid body.
        /// I = diag(2/3·mr², 2/3·mr², 2/3·mr²).
        /// </summary>
        public static RigidBody CreateHollowSphere(double mass, double radius)
        {
            double I = 2.0 / 3.0 * mass * radius * radius;
            return new RigidBody(mass, DiagonalInertia(I, I, I));
        }

        /// <summary>
        /// Creates a solid box (cuboid) rigid body.
        /// Ix = m/12·(h²+d²), Iy = m/12·(w²+d²), Iz = m/12·(w²+h²).
        /// </summary>
        public static RigidBody CreateSolidBox(double mass, double width, double height, double depth)
        {
            double w2 = width * width, h2 = height * height, d2 = depth * depth;
            double k = mass / 12.0;
            return new RigidBody(mass, DiagonalInertia(k * (h2 + d2), k * (w2 + d2), k * (w2 + h2)));
        }

        /// <summary>
        /// Creates a solid cylinder rigid body (symmetry axis along Z).
        /// Ix = Iy = m/12·(3r²+h²), Iz = ½mr².
        /// </summary>
        public static RigidBody CreateSolidCylinder(double mass, double radius, double height)
        {
            double r2 = radius * radius;
            double Idiameter = mass / 12.0 * (3 * r2 + height * height);
            double Iaxis = 0.5 * mass * r2;
            return new RigidBody(mass, DiagonalInertia(Idiameter, Idiameter, Iaxis));
        }

        /// <summary>
        /// Creates a hollow cylinder (tube) rigid body (symmetry axis along Z).
        /// Iz = ½m(r₁²+r₂²), Ix = Iy = m/12·(3(r₁²+r₂²)+h²).
        /// </summary>
        public static RigidBody CreateHollowCylinder(double mass, double innerRadius, double outerRadius, double height)
        {
            double r2sum = innerRadius * innerRadius + outerRadius * outerRadius;
            double Iaxis = 0.5 * mass * r2sum;
            double Idiameter = mass / 12.0 * (3 * r2sum + height * height);
            return new RigidBody(mass, DiagonalInertia(Idiameter, Idiameter, Iaxis));
        }

        /// <summary>
        /// Creates a thin rod rigid body (along Z-axis, about center).
        /// Ix = Iy = m·L²/12, Iz ≈ 0.
        /// </summary>
        public static RigidBody CreateThinRod(double mass, double length)
        {
            double I = mass * length * length / 12.0;
            return new RigidBody(mass, DiagonalInertia(I, I, 0));
        }

        private static Matrix DiagonalInertia(double Ix, double Iy, double Iz)
        {
            return new Matrix(new double[,]
            {
                { Ix, 0,  0  },
                { 0,  Iy, 0  },
                { 0,  0,  Iz }
            });
        }

        #endregion
    }
}
