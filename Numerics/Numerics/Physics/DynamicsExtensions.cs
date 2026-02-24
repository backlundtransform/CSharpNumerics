using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Objects;
using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;


namespace CSharpNumerics.Physics
{
    /// <summary>
    /// Provides extension methods for particle dynamics — forces, momentum, energy, work, and power.
    /// </summary>
    public static class DynamicsExtensions
    {
        #region Newton's Laws

        /// <summary>
        /// Computes acceleration from a net force and mass: a = F/m (Newton's 2nd law).
        /// </summary>
        /// <param name="force">Net force vector in Newtons.</param>
        /// <param name="mass">Mass in kg.</param>
        public static Vector Acceleration(this Vector force, double mass)
        {
            if (mass <= 0) throw new ArgumentException("Mass must be greater than zero.");
            return force / mass;
        }

        /// <summary>
        /// Computes the force required to produce a given acceleration: F = ma.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="acceleration">Acceleration vector in m/s².</param>
        public static Vector Force(this double mass, Vector acceleration)
        {
            return mass * acceleration;
        }

        /// <summary>
        /// Computes the net force from multiple force vectors: F_net = ΣF.
        /// </summary>
        /// <param name="first">First force vector.</param>
        /// <param name="rest">Additional force vectors.</param>
        public static Vector NetForce(this Vector first, params Vector[] rest)
        {
            var sum = first;
            for (int i = 0; i < rest.Length; i++)
                sum = sum + rest[i];
            return sum;
        }

        /// <summary>
        /// Computes the weight force of a mass near Earth's surface: W = m·g (downward along -Z).
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="direction">Direction of gravity (default: -Z).</param>
        public static Vector Weight(this double mass, Vector? direction = null)
        {
            var dir = direction ?? new Vector(0, 0, -1);
            return mass * PhysicsConstants.GravitationalAcceleration * dir.GetUnitVector();
        }

        #endregion

        #region Momentum & Impulse

        /// <summary>
        /// Computes linear momentum: p = mv.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="velocity">Velocity vector in m/s.</param>
        public static Vector Momentum(this double mass, Vector velocity)
        {
            return mass * velocity;
        }

        /// <summary>
        /// Computes scalar momentum: p = mv.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="velocity">Velocity in m/s.</param>
        public static double Momentum(this double mass, double velocity)
        {
            return mass * velocity;
        }

        /// <summary>
        /// Computes impulse from a constant force over time: J = F·Δt.
        /// </summary>
        /// <param name="force">Constant force vector in Newtons.</param>
        /// <param name="duration">Time duration in seconds.</param>
        public static Vector Impulse(this Vector force, double duration)
        {
            return duration * force;
        }

        /// <summary>
        /// Computes impulse as change in momentum: J = m·Δv.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="velocityBefore">Velocity before impact.</param>
        /// <param name="velocityAfter">Velocity after impact.</param>
        public static Vector ImpulseFromVelocityChange(this double mass, Vector velocityBefore, Vector velocityAfter)
        {
            return mass * (velocityAfter - velocityBefore);
        }

        /// <summary>
        /// Computes the velocity after applying an impulse: v' = v + J/m.
        /// </summary>
        /// <param name="impulse">Impulse vector in N·s.</param>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="currentVelocity">Current velocity vector.</param>
        public static Vector ApplyImpulse(this Vector impulse, double mass, Vector currentVelocity)
        {
            if (mass <= 0) throw new ArgumentException("Mass must be greater than zero.");
            return currentVelocity + impulse / mass;
        }

        #endregion

        #region Energy

        /// <summary>
        /// Computes kinetic energy: KE = ½mv².
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="velocity">Velocity vector in m/s.</param>
        public static double KineticEnergy(this double mass, Vector velocity)
        {
            return 0.5 * mass * velocity.Dot(velocity);
        }

        /// <summary>
        /// Computes kinetic energy from scalar speed: KE = ½mv².
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="speed">Speed in m/s.</param>
        public static double KineticEnergy(this double mass, double speed)
        {
            return 0.5 * mass * speed * speed;
        }

        /// <summary>
        /// Computes gravitational potential energy near Earth's surface: PE = mgh.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="height">Height above reference in meters.</param>
        public static double PotentialEnergy(this double mass, double height)
        {
            return mass * PhysicsConstants.GravitationalAcceleration * height;
        }

        /// <summary>
        /// Computes gravitational potential energy between two masses: U = -Gm₁m₂/r.
        /// </summary>
        /// <param name="mass1">First mass in kg.</param>
        /// <param name="mass2">Second mass in kg.</param>
        /// <param name="distance">Distance between centers in meters.</param>
        public static double GravitationalPotentialEnergy(this double mass1, double mass2, double distance)
        {
            if (distance <= 0) throw new ArgumentException("Distance must be greater than zero.");
            return -PhysicsConstants.GravitationalConstant * mass1 * mass2 / distance;
        }

        /// <summary>
        /// Computes the total mechanical energy: E = KE + PE = ½mv² + mgh.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="velocity">Velocity vector in m/s.</param>
        /// <param name="height">Height above reference in meters.</param>
        public static double MechanicalEnergy(this double mass, Vector velocity, double height)
        {
            return mass.KineticEnergy(velocity) + mass.PotentialEnergy(height);
        }

        /// <summary>
        /// Computes the speed required from energy conservation: v = √(2·KE/m).
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="kineticEnergy">Kinetic energy in Joules.</param>
        public static double SpeedFromKineticEnergy(this double mass, double kineticEnergy)
        {
            if (mass <= 0) throw new ArgumentException("Mass must be greater than zero.");
            if (kineticEnergy < 0) throw new ArgumentException("Kinetic energy cannot be negative.");
            return Math.Sqrt(2 * kineticEnergy / mass);
        }

        #endregion

        #region Work & Power

        /// <summary>
        /// Computes work done by a force over a displacement: W = F·d.
        /// </summary>
        /// <param name="force">Force vector in Newtons.</param>
        /// <param name="displacement">Displacement vector in meters.</param>
        public static double Work(this Vector force, Vector displacement)
        {
            return force.Dot(displacement);
        }

        /// <summary>
        /// Computes work done by a scalar force along a distance: W = F·d·cos(θ).
        /// </summary>
        /// <param name="force">Force magnitude in Newtons.</param>
        /// <param name="distance">Distance in meters.</param>
        /// <param name="angleRadians">Angle between force and displacement in radians (default 0 = parallel).</param>
        public static double Work(this double force, double distance, double angleRadians = 0.0)
        {
            return force * distance * Math.Cos(angleRadians);
        }

        /// <summary>
        /// Computes instantaneous power: P = F·v.
        /// </summary>
        /// <param name="force">Force vector in Newtons.</param>
        /// <param name="velocity">Velocity vector in m/s.</param>
        public static double Power(this Vector force, Vector velocity)
        {
            return force.Dot(velocity);
        }

        /// <summary>
        /// Computes average power: P = W/Δt.
        /// </summary>
        /// <param name="work">Work done in Joules.</param>
        /// <param name="duration">Time duration in seconds.</param>
        public static double AveragePower(this double work, double duration)
        {
            if (duration <= 0) throw new ArgumentException("Duration must be greater than zero.");
            return work / duration;
        }

        /// <summary>
        /// Computes the work–energy theorem result: the change in kinetic energy.
        /// ΔKE = ½m(v₂² - v₁²).
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="velocityBefore">Velocity before.</param>
        /// <param name="velocityAfter">Velocity after.</param>
        public static double WorkEnergyTheorem(this double mass, Vector velocityBefore, Vector velocityAfter)
        {
            return 0.5 * mass * (velocityAfter.Dot(velocityAfter) - velocityBefore.Dot(velocityBefore));
        }

        #endregion

        #region Elastic & Inelastic Collisions (1D)

        /// <summary>
        /// Computes the final velocity of object 1 in a 1D elastic collision.
        /// v₁' = ((m₁-m₂)v₁ + 2m₂v₂) / (m₁+m₂).
        /// </summary>
        /// <param name="m1">Mass of object 1 in kg.</param>
        /// <param name="v1">Initial velocity of object 1.</param>
        /// <param name="m2">Mass of object 2 in kg.</param>
        /// <param name="v2">Initial velocity of object 2.</param>
        public static double ElasticCollisionVelocity(this double m1, double v1, double m2, double v2)
        {
            return ((m1 - m2) * v1 + 2 * m2 * v2) / (m1 + m2);
        }

        /// <summary>
        /// Computes the final velocities of both objects in a 1D elastic collision.
        /// </summary>
        /// <param name="m1">Mass of object 1 in kg.</param>
        /// <param name="v1">Initial velocity of object 1.</param>
        /// <param name="m2">Mass of object 2 in kg.</param>
        /// <param name="v2">Initial velocity of object 2.</param>
        /// <returns>Tuple of (v1', v2') final velocities.</returns>
        public static (double v1Final, double v2Final) ElasticCollision(this double m1, double v1, double m2, double v2)
        {
            double totalMass = m1 + m2;
            double v1f = ((m1 - m2) * v1 + 2 * m2 * v2) / totalMass;
            double v2f = ((m2 - m1) * v2 + 2 * m1 * v1) / totalMass;
            return (v1f, v2f);
        }

        /// <summary>
        /// Computes the final velocity of a perfectly inelastic (sticky) collision: v' = (m₁v₁+m₂v₂)/(m₁+m₂).
        /// </summary>
        /// <param name="m1">Mass of object 1 in kg.</param>
        /// <param name="v1">Velocity of object 1 before collision.</param>
        /// <param name="m2">Mass of object 2 in kg.</param>
        /// <param name="v2">Velocity of object 2 before collision.</param>
        public static double InelasticCollisionVelocity(this double m1, double v1, double m2, double v2)
        {
            return (m1 * v1 + m2 * v2) / (m1 + m2);
        }

        /// <summary>
        /// Computes the final velocity vector of a perfectly inelastic collision in 3D.
        /// v' = (m₁v₁ + m₂v₂) / (m₁+m₂).
        /// </summary>
        /// <param name="m1">Mass of object 1 in kg.</param>
        /// <param name="v1">Velocity vector of object 1.</param>
        /// <param name="m2">Mass of object 2 in kg.</param>
        /// <param name="v2">Velocity vector of object 2.</param>
        public static Vector InelasticCollisionVelocity(this double m1, Vector v1, double m2, Vector v2)
        {
            return (1.0 / (m1 + m2)) * (m1 * v1 + m2 * v2);
        }

        /// <summary>
        /// Computes energy lost in a perfectly inelastic collision.
        /// </summary>
        /// <param name="m1">Mass of object 1 in kg.</param>
        /// <param name="v1">Velocity of object 1 before collision.</param>
        /// <param name="m2">Mass of object 2 in kg.</param>
        /// <param name="v2">Velocity of object 2 before collision.</param>
        public static double InelasticCollisionEnergyLoss(this double m1, double v1, double m2, double v2)
        {
            double keBefore = 0.5 * m1 * v1 * v1 + 0.5 * m2 * v2 * v2;
            double vFinal = m1.InelasticCollisionVelocity(v1, m2, v2);
            double keAfter = 0.5 * (m1 + m2) * vFinal * vFinal;
            return keBefore - keAfter;
        }

        /// <summary>
        /// Computes the coefficient of restitution from before/after velocities.
        /// e = -(v₁' - v₂') / (v₁ - v₂).
        /// </summary>
        /// <param name="v1Before">Velocity of object 1 before collision.</param>
        /// <param name="v2Before">Velocity of object 2 before collision.</param>
        /// <param name="v1After">Velocity of object 1 after collision.</param>
        /// <param name="v2After">Velocity of object 2 after collision.</param>
        public static double CoefficientOfRestitution(
            this double v1Before,
            double v2Before,
            double v1After,
            double v2After)
        {
            double approach = v1Before - v2Before;
            if (Math.Abs(approach) < 1e-15) throw new ArgumentException("Objects have the same initial velocity.");
            return -(v1After - v2After) / approach;
        }

        #endregion

        #region Moment of Inertia (Scalar)

        /// <summary>
        /// Moment of inertia of a solid sphere: I = 2/5·mr².
        /// </summary>
        public static double MomentOfInertiaSolidSphere(this double mass, double radius)
            => 2.0 / 5.0 * mass * radius * radius;

        /// <summary>
        /// Moment of inertia of a hollow sphere (thin shell): I = 2/3·mr².
        /// </summary>
        public static double MomentOfInertiaHollowSphere(this double mass, double radius)
            => 2.0 / 3.0 * mass * radius * radius;

        /// <summary>
        /// Moment of inertia of a solid cylinder about its symmetry axis: I = ½mr².
        /// </summary>
        public static double MomentOfInertiaSolidCylinder(this double mass, double radius)
            => 0.5 * mass * radius * radius;

        /// <summary>
        /// Moment of inertia of a hollow cylinder about its symmetry axis: I = ½m(r₁²+r₂²).
        /// </summary>
        public static double MomentOfInertiaHollowCylinder(this double mass, double innerRadius, double outerRadius)
            => 0.5 * mass * (innerRadius * innerRadius + outerRadius * outerRadius);

        /// <summary>
        /// Moment of inertia of a thin rod about its center: I = mL²/12.
        /// </summary>
        public static double MomentOfInertiaThinRod(this double mass, double length)
            => mass * length * length / 12.0;

        /// <summary>
        /// Moment of inertia of a thin rod about one end: I = mL²/3.
        /// </summary>
        public static double MomentOfInertiaThinRodEnd(this double mass, double length)
            => mass * length * length / 3.0;

        /// <summary>
        /// Moment of inertia of a solid box about an axis through the center perpendicular to face (a×b):
        /// I = m/12·(a²+b²).
        /// </summary>
        public static double MomentOfInertiaSolidBox(this double mass, double sideA, double sideB)
            => mass / 12.0 * (sideA * sideA + sideB * sideB);

        /// <summary>
        /// Parallel axis theorem (scalar): I_new = I_cm + m·d².
        /// </summary>
        /// <param name="momentOfInertia">Moment of inertia about center of mass.</param>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="distance">Distance from CM to new axis in meters.</param>
        public static double ParallelAxis(this double momentOfInertia, double mass, double distance)
            => momentOfInertia + mass * distance * distance;

        #endregion

        #region Inertia Tensor (3×3 Matrix)

        /// <summary>
        /// 3×3 inertia tensor of a solid sphere: I = diag(2/5·mr²).
        /// </summary>
        public static Matrix InertiaTensorSolidSphere(this double mass, double radius)
        {
            double I = 2.0 / 5.0 * mass * radius * radius;
            return new Matrix(new double[,] { { I, 0, 0 }, { 0, I, 0 }, { 0, 0, I } });
        }

        /// <summary>
        /// 3×3 inertia tensor of a hollow sphere: I = diag(2/3·mr²).
        /// </summary>
        public static Matrix InertiaTensorHollowSphere(this double mass, double radius)
        {
            double I = 2.0 / 3.0 * mass * radius * radius;
            return new Matrix(new double[,] { { I, 0, 0 }, { 0, I, 0 }, { 0, 0, I } });
        }

        /// <summary>
        /// 3×3 inertia tensor of a solid box (cuboid) about its center.
        /// </summary>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="width">Width along X.</param>
        /// <param name="height">Height along Y.</param>
        /// <param name="depth">Depth along Z.</param>
        public static Matrix InertiaTensorSolidBox(this double mass, double width, double height, double depth)
        {
            double w2 = width * width, h2 = height * height, d2 = depth * depth;
            double k = mass / 12.0;
            return new Matrix(new double[,]
            {
                { k * (h2 + d2), 0, 0 },
                { 0, k * (w2 + d2), 0 },
                { 0, 0, k * (w2 + h2) }
            });
        }

        /// <summary>
        /// 3×3 inertia tensor of a solid cylinder (symmetry axis along Z).
        /// </summary>
        public static Matrix InertiaTensorSolidCylinder(this double mass, double radius, double height)
        {
            double r2 = radius * radius;
            double Idiameter = mass / 12.0 * (3 * r2 + height * height);
            double Iaxis = 0.5 * mass * r2;
            return new Matrix(new double[,]
            {
                { Idiameter, 0, 0 },
                { 0, Idiameter, 0 },
                { 0, 0, Iaxis }
            });
        }

        /// <summary>
        /// Parallel axis theorem for a 3×3 inertia tensor: I_new = I_cm + m·(d²·E - d⊗d).
        /// </summary>
        /// <param name="inertiaTensor">Inertia tensor about center of mass.</param>
        /// <param name="mass">Mass in kg.</param>
        /// <param name="offset">Displacement vector from CM to new point.</param>
        public static Matrix ParallelAxis(this Matrix inertiaTensor, double mass, Vector offset)
        {
            double d2 = offset.Dot(offset);
            double[] d = [offset.x, offset.y, offset.z];
            var result = new double[3, 3];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    result[i, j] = inertiaTensor.values[i, j]
                        + mass * (d2 * (i == j ? 1.0 : 0.0) - d[i] * d[j]);
                }
            }
            return new Matrix(result);
        }

        #endregion

        #region Torque & Rotational Dynamics

        /// <summary>
        /// Computes torque from a moment arm and force: τ = r × F.
        /// </summary>
        /// <param name="momentArm">Vector from pivot to point of force application.</param>
        /// <param name="force">Force vector in Newtons.</param>
        public static Vector Torque(this Vector momentArm, Vector force)
        {
            return momentArm.Cross(force);
        }

        /// <summary>
        /// Computes angular momentum: L = I·ω.
        /// </summary>
        /// <param name="inertiaTensor">3×3 inertia tensor.</param>
        /// <param name="angularVelocity">Angular velocity vector in rad/s.</param>
        public static Vector AngularMomentum(this Matrix inertiaTensor, Vector angularVelocity)
        {
            return inertiaTensor * angularVelocity;
        }

        /// <summary>
        /// Computes angular acceleration: α = I⁻¹·τ.
        /// </summary>
        /// <param name="inverseInertiaTensor">Inverse of the 3×3 inertia tensor.</param>
        /// <param name="torque">Torque vector in N·m.</param>
        public static Vector AngularAcceleration(this Matrix inverseInertiaTensor, Vector torque)
        {
            return inverseInertiaTensor * torque;
        }

        /// <summary>
        /// Computes rotational kinetic energy: KE_rot = ½ωᵀIω.
        /// </summary>
        /// <param name="inertiaTensor">3×3 inertia tensor.</param>
        /// <param name="angularVelocity">Angular velocity vector in rad/s.</param>
        public static double RotationalKineticEnergy(this Matrix inertiaTensor, Vector angularVelocity)
        {
            var Iw = inertiaTensor * angularVelocity;
            return 0.5 * angularVelocity.Dot(Iw);
        }

        /// <summary>
        /// Computes the scalar torque magnitude for rotation about a single axis: τ = r·F·sin(θ).
        /// </summary>
        /// <param name="force">Force magnitude in Newtons.</param>
        /// <param name="momentArm">Distance from pivot in meters.</param>
        /// <param name="angleRadians">Angle between force and moment arm in radians.</param>
        public static double Torque(this double force, double momentArm, double angleRadians = Math.PI / 2)
        {
            return force * momentArm * Math.Sin(angleRadians);
        }

        /// <summary>
        /// Computes scalar angular momentum: L = I·ω.
        /// </summary>
        public static double AngularMomentum(this double momentOfInertia, double angularVelocity)
        {
            return momentOfInertia * angularVelocity;
        }

        /// <summary>
        /// Computes scalar rotational kinetic energy: KE = ½Iω².
        /// </summary>
        public static double RotationalKineticEnergy(this double momentOfInertia, double angularVelocity)
        {
            return 0.5 * momentOfInertia * angularVelocity * angularVelocity;
        }

        #endregion

        #region RigidBody Integration

        /// <summary>
        /// Integrates a rigid body one time step using semi-implicit (symplectic) Euler.
        /// Updates velocity before position for better stability in oscillatory systems.
        /// Forces and torques must be applied before calling this method.
        /// Accumulators are cleared after integration.
        /// Delegates to <see cref="DifferentialEquationExtensions.SemiImplicitEuler"/> for the linear state.
        /// </summary>
        /// <param name="body">The rigid body to integrate.</param>
        /// <param name="dt">Time step in seconds.</param>
        public static void IntegrateSemiImplicitEuler(this ref RigidBody body, double dt)
        {
            if (body.IsStatic) { body.ClearForces(); return; }

            var linearAccel = body.LinearAcceleration;
            var angularAccel = body.AngularAcceleration;

            // Pack linear state: [px, py, pz, vx, vy, vz] (half = 3)
            var y0 = new VectorN([
                body.Position.x, body.Position.y, body.Position.z,
                body.Velocity.x, body.Velocity.y, body.Velocity.z
            ]);

            Func<(double t, VectorN y), VectorN> derivative = state =>
                new VectorN([
                    state.y[3], state.y[4], state.y[5],
                    linearAccel.x, linearAccel.y, linearAccel.z
                ]);

            var result = derivative.SemiImplicitEuler(0, dt, dt, y0);

            body.Position = new Vector(result[0], result[1], result[2]);
            body.Velocity = new Vector(result[3], result[4], result[5]);
            body.AngularVelocity = body.AngularVelocity + dt * angularAccel;

            body.ClearForces();
        }

        /// <summary>
        /// Integrates a rigid body one time step using velocity Verlet.
        /// Provides O(dt²) accuracy and excellent energy conservation.
        /// The forceFunc computes all forces/torques at each evaluation — no need to pre-apply forces.
        /// Delegates to <see cref="DifferentialEquationExtensions.VelocityVerlet"/> for the state update.
        /// </summary>
        /// <param name="body">The rigid body to integrate.</param>
        /// <param name="forceFunc">Function that computes (force, torque) given the body's current state.</param>
        /// <param name="dt">Time step in seconds.</param>
        public static void IntegrateVelocityVerlet(
            this ref RigidBody body,
            Func<RigidBody, (Vector force, Vector torque)> forceFunc,
            double dt)
        {
            if (body.IsStatic) { body.ClearForces(); return; }

            var invMass = body.InverseMass;
            var invI = body.InverseInertiaTensorWorld;
            var mass = body.Mass;
            var inertiaTensor = body.InertiaTensorBody;

            // Pack state: [px, py, pz, θx, θy, θz, vx, vy, vz, ωx, ωy, ωz] (half = 6)
            // θ slots are dummy angular positions (not tracked by RigidBody, but needed for Verlet layout)
            var y0 = new VectorN([
                body.Position.x, body.Position.y, body.Position.z,
                0, 0, 0,
                body.Velocity.x, body.Velocity.y, body.Velocity.z,
                body.AngularVelocity.x, body.AngularVelocity.y, body.AngularVelocity.z
            ]);

            Func<(double t, VectorN y), VectorN> accelFunc = state =>
            {
                var tempBody = new RigidBody(mass, inertiaTensor)
                {
                    Position = new Vector(state.y[0], state.y[1], state.y[2]),
                    Velocity = new Vector(state.y[6], state.y[7], state.y[8]),
                    AngularVelocity = new Vector(state.y[9], state.y[10], state.y[11])
                };
                var (f, torque) = forceFunc(tempBody);
                var linA = invMass * f;
                var angA = invI * torque;
                return new VectorN([linA.x, linA.y, linA.z, angA.x, angA.y, angA.z]);
            };

            var result = accelFunc.VelocityVerlet(0, dt, dt, y0);

            body.Position = new Vector(result[0], result[1], result[2]);
            body.Velocity = new Vector(result[6], result[7], result[8]);
            body.AngularVelocity = new Vector(result[9], result[10], result[11]);

            body.ClearForces();
        }

        /// <summary>
        /// Integrates a rigid body one time step using explicit Euler.
        /// Simple but least stable — use semi-implicit Euler or Verlet for better results.
        /// Delegates to <see cref="DifferentialEquationExtensions.EulerMethod"/> for the linear state.
        /// </summary>
        /// <param name="body">The rigid body to integrate.</param>
        /// <param name="dt">Time step in seconds.</param>
        public static void IntegrateEuler(this ref RigidBody body, double dt)
        {
            if (body.IsStatic) { body.ClearForces(); return; }

            var linearAccel = body.LinearAcceleration;
            var angularAccel = body.AngularAcceleration;

            // Pack linear state: [px, py, pz, vx, vy, vz]
            var y0 = new VectorN([
                body.Position.x, body.Position.y, body.Position.z,
                body.Velocity.x, body.Velocity.y, body.Velocity.z
            ]);

            Func<(double t, VectorN y), VectorN> derivative = state =>
                new VectorN([
                    state.y[3], state.y[4], state.y[5],
                    linearAccel.x, linearAccel.y, linearAccel.z
                ]);

            var result = derivative.EulerMethod(0, dt, dt, y0);

            body.Position = new Vector(result[0], result[1], result[2]);
            body.Velocity = new Vector(result[3], result[4], result[5]);
            body.AngularVelocity = body.AngularVelocity + dt * angularAccel;

            body.ClearForces();
        }

        #endregion

        #region Common Force Models

        /// <summary>
        /// Computes a spring force between two points using Hooke's law: F = -k·(|Δr| - L₀)·r̂.
        /// The returned force acts on the object at <paramref name="position"/> pulling it toward <paramref name="anchor"/>.
        /// </summary>
        /// <param name="k">Spring stiffness in N/m.</param>
        /// <param name="restLength">Natural (rest) length of the spring in meters.</param>
        /// <param name="position">Position of the attached object.</param>
        /// <param name="anchor">Position of the anchor point.</param>
        public static Vector SpringForce(this double k, double restLength, Vector position, Vector anchor)
        {
            var delta = position - anchor;
            double distance = delta.GetMagnitude();
            if (distance < 1e-15) return new Vector(0, 0, 0);
            var direction = (1.0 / distance) * delta;
            return -k * (distance - restLength) * direction;
        }

        /// <summary>
        /// Computes a linear damping (viscous) force: F = -c·v.
        /// Opposes the direction of motion proportionally to speed.
        /// </summary>
        /// <param name="c">Damping coefficient in N·s/m.</param>
        /// <param name="velocity">Velocity vector of the object in m/s.</param>
        public static Vector DampingForce(this double c, Vector velocity)
        {
            return -c * velocity;
        }

        /// <summary>
        /// Computes aerodynamic drag: F = -½·Cd·ρ·A·|v|·v.
        /// The force is proportional to v² and opposes the direction of motion.
        /// </summary>
        /// <param name="dragCoefficient">Drag coefficient Cd (dimensionless).</param>
        /// <param name="fluidDensity">Fluid density ρ in kg/m³ (air ≈ 1.225).</param>
        /// <param name="crossSectionArea">Cross-sectional area A in m².</param>
        /// <param name="velocity">Velocity vector of the object in m/s.</param>
        public static Vector DragForce(this double dragCoefficient, double fluidDensity, double crossSectionArea, Vector velocity)
        {
            double speed = velocity.GetMagnitude();
            if (speed < 1e-15) return new Vector(0, 0, 0);
            return -0.5 * dragCoefficient * fluidDensity * crossSectionArea * speed * velocity;
        }

        /// <summary>
        /// Computes a friction force given a friction coefficient and normal force magnitude.
        /// Static friction is returned when speed is below <paramref name="staticThreshold"/>;
        /// kinetic friction is returned otherwise.
        /// 
        /// Kinetic: F = -μk·|N|·v̂
        /// Static:  F = -min(|applied|, μs·|N|)·applied̂  (opposes applied tangential force)
        /// </summary>
        /// <param name="mu">Friction coefficient (μk for kinetic, or μs for static).</param>
        /// <param name="normalForceMagnitude">Magnitude of the normal force |N| in Newtons.</param>
        /// <param name="velocity">Current velocity of the object in m/s.</param>
        /// <param name="appliedTangentialForce">
        /// The tangential component of the applied force (used for static friction direction).
        /// If null, only kinetic friction is computed.
        /// </param>
        /// <param name="staticThreshold">Speed below which static friction applies (default 1e-6 m/s).</param>
        public static Vector FrictionForce(
            this double mu,
            double normalForceMagnitude,
            Vector velocity,
            Vector? appliedTangentialForce = null,
            double staticThreshold = 1e-6)
        {
            if (normalForceMagnitude < 0) normalForceMagnitude = 0;
            double speed = velocity.GetMagnitude();

            if (speed > staticThreshold)
            {
                // Kinetic friction: opposes velocity direction
                return -mu * normalForceMagnitude * (1.0 / speed) * velocity;
            }

            // Static friction: opposes the applied tangential force
            if (appliedTangentialForce is Vector applied)
            {
                double appliedMag = applied.GetMagnitude();
                if (appliedMag < 1e-15) return new Vector(0, 0, 0);
                double maxStaticForce = mu * normalForceMagnitude;
                double frictionMag = Math.Min(appliedMag, maxStaticForce);
                return -frictionMag * (1.0 / appliedMag) * applied;
            }

            return new Vector(0, 0, 0);
        }

        #endregion
    }
}
