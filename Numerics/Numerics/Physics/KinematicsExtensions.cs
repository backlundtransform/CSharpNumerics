using CSharpNumerics.Physics.Constants;
using System;
using Numerics.Objects;

namespace CSharpNumerics.Physics
{
    /// <summary>
    /// Provides extension methods for kinematic calculations involving scalars and vectors.
    /// </summary>
    public static class KinematicsExtensions
    {
        #region Free Fall

        /// <summary>
        /// Computes velocity from free fall from height h: v = sqrt(2 * g * h).
        /// </summary>
        /// <param name="height">Height in meters.</param>
        public static double FreeFallVelocity(this double height)
        {
            return Math.Sqrt(2 * PhysicsConstants.GravitationalAcceleration * height);
        }
        /// <summary>
        /// Computes free fall velocity vector (assuming downward along a given direction).
        /// </summary>
        /// <param name="height">Height in meters (scalar magnitude).</param>
        /// <param name="direction">Unit vector indicating direction of fall (default: downward -Z).</param>
        public static Vector FreeFallVelocity(this double height, Vector? direction = null)
        {
            var dir = direction ?? new Vector(0, 0, -1); 
            return  Math.Sqrt(2 * PhysicsConstants.GravitationalAcceleration * height) * dir.GetUnitVector();
        }

        /// <summary>
        /// Computes the time it takes to fall from a certain height: t = sqrt(2h / g).
        /// </summary>
        /// <param name="height">Height in meters.</param>
        public static double FreeFallTime(this double height)
        {
            return Math.Sqrt((2 * height) / PhysicsConstants.GravitationalAcceleration);
        }

        #endregion

        #region Constant Velocity

        /// <summary>
        /// Computes position assuming constant velocity: x = x0 + v * t.
        /// </summary>
        public static double PositionFromConstantVelocity(
            this double velocity,
            double time,
            double initialPosition = 0.0)
        {
            return initialPosition + velocity * time;
        }

        /// <summary>
        /// Computes position vector assuming constant velocity: r = r0 + v * t.
        /// </summary>
        public static Vector PositionFromConstantVelocity(
            this Vector velocity,
            double time,
            Vector initialPosition)
        {
            return initialPosition + time * velocity;
        }

        #endregion

        #region Constant Acceleration

        /// <summary>
        /// Computes velocity assuming constant acceleration: v = v0 + a * t.
        /// </summary>
        public static double VelocityFromConstantAcceleration(
            this double acceleration,
            double time,
            double initialVelocity = 0.0)
        {
            return initialVelocity + acceleration * time;
        }

        /// <summary>
        /// Computes velocity vector assuming constant acceleration: v = v0 + a * t.
        /// </summary>
        public static Vector VelocityFromConstantAcceleration(
            this Vector acceleration,
            double time,
            Vector initialVelocity)
        {
            return initialVelocity + time * acceleration;
        }

        /// <summary>
        /// Computes position assuming constant acceleration: s = s0 + v0*t + 0.5*a*t².
        /// </summary>
        public static double PositionFromConstantAcceleration(
            this double acceleration,
            double time,
            double initialVelocity = 0.0,
            double initialPosition = 0.0)
        {
            return initialPosition
                 + initialVelocity * time
                 + 0.5 * acceleration * time * time;
        }

        /// <summary>
        /// Computes position vector assuming constant acceleration: r = r0 + v0*t + 0.5*a*t².
        /// </summary>
        public static Vector PositionFromConstantAcceleration(
            this Vector acceleration,
            double time,
            Vector initialVelocity,
            Vector initialPosition)
        {
            return initialPosition
                 + time * initialVelocity
                 + 0.5 * time * time * acceleration;
        }

        #endregion

        #region Time-Independent & Helper Equations (SUVAT)

        /// <summary>
        /// Computes final velocity when time is unknown (Torricelli's Law): v = sqrt(v0² + 2*a*Δs).
        /// </summary>
        public static double VelocityFromDisplacement(
            this double acceleration,
            double displacement,
            double initialVelocity = 0.0)
        {
            return Math.Sqrt(Math.Pow(initialVelocity, 2) + 2 * acceleration * displacement);
        }

        /// <summary>
        /// Computes final velocity vector when time is unknown: v = sqrt(v0² + 2*a*Δs) in each component.
        /// </summary>
        public static Vector VelocityFromDisplacement(this Vector acceleration, Vector displacement, Vector? initialVelocity = null)
        {
            var v0 = initialVelocity ?? new Vector(0, 0, 0);
            return new Vector(
                Math.Sqrt(v0.x * v0.x + 2 * acceleration.x * displacement.x),
                Math.Sqrt(v0.y * v0.y + 2 * acceleration.y * displacement.y),
                Math.Sqrt(v0.z * v0.z + 2 * acceleration.z * displacement.z)
            );
        }

        /// <summary>
        /// Computes displacement when time is unknown: Δs = (v² - v0²) / (2*a).
        /// </summary>
        public static double DisplacementFromVelocities(
            this double acceleration,
            double finalVelocity,
            double initialVelocity = 0.0)
        {
            if (Math.Abs(acceleration) < 1e-10) throw new DivideByZeroException("Acceleration is too close to zero.");
            return (Math.Pow(finalVelocity, 2) - Math.Pow(initialVelocity, 2)) / (2 * acceleration);
        }

        /// <summary>
        /// Computes displacement vector from velocities: Δs = (v² - v0²) / (2*a) component-wise.
        /// </summary>
        public static Vector DisplacementFromVelocities(this Vector acceleration, Vector finalVelocity, Vector? initialVelocity = null)
        {
            var v0 = initialVelocity ?? new Vector(0, 0, 0);
            return new Vector(
                acceleration.x.Abs() < 1e-10 ? 0 : (finalVelocity.x * finalVelocity.x - v0.x * v0.x) / (2 * acceleration.x),
                acceleration.y.Abs() < 1e-10 ? 0 : (finalVelocity.y * finalVelocity.y - v0.y * v0.y) / (2 * acceleration.y),
                acceleration.z.Abs() < 1e-10 ? 0 : (finalVelocity.z * finalVelocity.z - v0.z * v0.z) / (2 * acceleration.z)
            );
        }

        /// <summary>
        /// Computes time required to reach a specific velocity: t = (v - v0) / a.
        /// </summary>
        public static double TimeToReachVelocity(
            this double acceleration,
            double finalVelocity,
            double initialVelocity = 0.0)
        {
            if (Math.Abs(acceleration) < 1e-10) throw new DivideByZeroException("Acceleration is too close to zero.");
            return (finalVelocity - initialVelocity) / acceleration;
        }
        /// <summary>
        /// Computes time to reach a final velocity vector component-wise: t = (v - v0)/a.
        /// </summary>
        public static Vector TimeToReachVelocity(this Vector acceleration, Vector finalVelocity, Vector? initialVelocity = null)
        {
            var v0 = initialVelocity ?? new Vector(0, 0, 0);
            return new Vector(
                acceleration.x.Abs() < 1e-10 ? 0 : (finalVelocity.x - v0.x) / acceleration.x,
                acceleration.y.Abs() < 1e-10 ? 0 : (finalVelocity.y - v0.y) / acceleration.y,
                acceleration.z.Abs() < 1e-10 ? 0 : (finalVelocity.z - v0.z) / acceleration.z
            );
        }

        /// <summary>
        /// Computes displacement using average velocity: Δs = t * (v0 + v) / 2.
        /// </summary>
        public static double DisplacementFromAverageVelocity(
            this double time,
            double initialVelocity,
            double finalVelocity)
        {
            return time * (initialVelocity + finalVelocity) / 2.0;
        }

        /// <summary>
        /// Computes displacement using average velocity vector: Δs = t * (v0 + v)/2 component-wise.
        /// </summary>
        public static Vector DisplacementFromAverageVelocity(this double time, Vector initialVelocity, Vector finalVelocity)
        {
            return 0.5 * time * (initialVelocity + finalVelocity);
        }

        #endregion

        #region Circular Motion

        /// <summary>
        /// Computes centripetal acceleration: a = v² / r.
        /// </summary>
        /// <param name="velocity">Tangential velocity.</param>
        /// <param name="radius">Radius of the circular path.</param>
        public static double CentripetalAcceleration(this double velocity, double radius)
        {
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return Math.Pow(velocity, 2) / radius;
        }
        /// <summary>
        /// Computes centripetal acceleration vector: a = v² / r, direction toward center.
        /// </summary>
        /// <param name="velocityVector">Tangential velocity vector</param>
        /// <param name="radiusVector">Vector from center to object</param>
        public static Vector CentripetalAcceleration(this Vector velocityVector, Vector radiusVector)
        {
            var r = radiusVector.GetMagnitude();
            if (r <= 0) throw new ArgumentException("Radius must be greater than zero.");
            var vMag = velocityVector.GetMagnitude();
            return -(vMag * vMag / r) * radiusVector.GetUnitVector(); // acceleration toward center
        }
        #endregion
    }
}