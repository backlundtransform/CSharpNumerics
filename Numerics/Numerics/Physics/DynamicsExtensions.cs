using CSharpNumerics.Physics.Constants;
using System;
using Numerics.Objects;

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
    }
}
