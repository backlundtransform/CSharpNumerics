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

        #region Projectile Motion

        /// <summary>
        /// Creates an initial velocity vector from launch speed and angle in the XZ plane.
        /// v₀ = (v·cos(θ), 0, v·sin(θ)).
        /// </summary>
        /// <param name="speed">Launch speed in m/s.</param>
        /// <param name="launchAngleRadians">Launch angle from horizontal in radians.</param>
        public static Vector ProjectileVelocityFromAngle(this double speed, double launchAngleRadians)
        {
            return new Vector(
                speed * Math.Cos(launchAngleRadians),
                0,
                speed * Math.Sin(launchAngleRadians)
            );
        }

        /// <summary>
        /// Computes the position vector of a projectile at time t.
        /// Assumes gravity acts along -Z: r(t) = v₀·t + ½g·t² with g = (0, 0, -g).
        /// </summary>
        /// <param name="initialVelocity">Initial velocity vector of the projectile.</param>
        /// <param name="time">Time elapsed since launch in seconds.</param>
        /// <param name="initialHeight">Launch height above ground in meters (default 0).</param>
        public static Vector ProjectilePosition(this Vector initialVelocity, double time, double initialHeight = 0.0)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            return new Vector(
                initialVelocity.x * time,
                initialVelocity.y * time,
                initialHeight + initialVelocity.z * time - 0.5 * g * time * time
            );
        }

        /// <summary>
        /// Computes the velocity vector of a projectile at time t.
        /// v(t) = v₀ + g·t where g = (0, 0, -g).
        /// </summary>
        /// <param name="initialVelocity">Initial velocity vector of the projectile.</param>
        /// <param name="time">Time elapsed since launch in seconds.</param>
        public static Vector ProjectileVelocity(this Vector initialVelocity, double time)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            return new Vector(
                initialVelocity.x,
                initialVelocity.y,
                initialVelocity.z - g * time
            );
        }

        /// <summary>
        /// Computes the total time of flight for a projectile until it returns to ground (z = 0).
        /// T = (v₀z + √(v₀z² + 2g·h₀)) / g.
        /// </summary>
        /// <param name="initialVelocity">Initial velocity vector of the projectile.</param>
        /// <param name="initialHeight">Launch height above ground in meters (default 0).</param>
        public static double ProjectileTimeOfFlight(this Vector initialVelocity, double initialHeight = 0.0)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            double vz = initialVelocity.z;
            return (vz + Math.Sqrt(vz * vz + 2 * g * initialHeight)) / g;
        }

        /// <summary>
        /// Computes the maximum height reached by a projectile: H = h₀ + v₀z² / (2g).
        /// </summary>
        /// <param name="initialVelocity">Initial velocity vector of the projectile.</param>
        /// <param name="initialHeight">Launch height above ground in meters (default 0).</param>
        public static double ProjectileMaxHeight(this Vector initialVelocity, double initialHeight = 0.0)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            double vz = initialVelocity.z;
            return initialHeight + (vz * vz) / (2 * g);
        }

        /// <summary>
        /// Computes the horizontal range of a projectile (distance in XY plane when z returns to 0).
        /// Range = √(v₀x² + v₀y²) · T.
        /// </summary>
        /// <param name="initialVelocity">Initial velocity vector of the projectile.</param>
        /// <param name="initialHeight">Launch height above ground in meters (default 0).</param>
        public static double ProjectileRange(this Vector initialVelocity, double initialHeight = 0.0)
        {
            double T = initialVelocity.ProjectileTimeOfFlight(initialHeight);
            double horizontalSpeed = Math.Sqrt(
                initialVelocity.x * initialVelocity.x + initialVelocity.y * initialVelocity.y);
            return horizontalSpeed * T;
        }

        /// <summary>
        /// Computes time of flight for a projectile from speed and angle: T = 2v₀sin(θ)/g.
        /// Assumes launch and landing at the same height.
        /// </summary>
        /// <param name="speed">Launch speed in m/s.</param>
        /// <param name="launchAngleRadians">Launch angle from horizontal in radians.</param>
        public static double ProjectileTimeOfFlight(this double speed, double launchAngleRadians)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            return 2 * speed * Math.Sin(launchAngleRadians) / g;
        }

        /// <summary>
        /// Computes the maximum height of a projectile: H = v₀²sin²(θ) / (2g).
        /// </summary>
        /// <param name="speed">Launch speed in m/s.</param>
        /// <param name="launchAngleRadians">Launch angle from horizontal in radians.</param>
        public static double ProjectileMaxHeight(this double speed, double launchAngleRadians)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            double vz = speed * Math.Sin(launchAngleRadians);
            return (vz * vz) / (2 * g);
        }

        /// <summary>
        /// Computes the horizontal range of a projectile: R = v₀²sin(2θ) / g.
        /// Assumes launch and landing at the same height.
        /// </summary>
        /// <param name="speed">Launch speed in m/s.</param>
        /// <param name="launchAngleRadians">Launch angle from horizontal in radians.</param>
        public static double ProjectileRange(this double speed, double launchAngleRadians)
        {
            double g = PhysicsConstants.GravitationalAcceleration;
            return (speed * speed * Math.Sin(2 * launchAngleRadians)) / g;
        }

        #endregion
    }
}