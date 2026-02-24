using CSharpNumerics.Physics.Constants;
using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;

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

        /// <summary>
        /// Computes angular speed from tangential speed and radius: ω = v / r.
        /// </summary>
        /// <param name="speed">Tangential speed in m/s.</param>
        /// <param name="radius">Radius of the circular path in meters.</param>
        public static double AngularSpeed(this double speed, double radius)
        {
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return speed / radius;
        }

        /// <summary>
        /// Computes the angular velocity vector: ω = (r × v) / |r|².
        /// The result points along the axis of rotation (right-hand rule).
        /// </summary>
        /// <param name="tangentialVelocity">Tangential velocity vector.</param>
        /// <param name="radiusVector">Vector from center to object.</param>
        public static Vector AngularVelocity(this Vector tangentialVelocity, Vector radiusVector)
        {
            var r2 = radiusVector.Dot(radiusVector);
            if (r2 <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return (1.0 / r2) * radiusVector.Cross(tangentialVelocity);
        }

        /// <summary>
        /// Computes the period of uniform circular motion: T = 2πr / v.
        /// </summary>
        /// <param name="speed">Tangential speed in m/s.</param>
        /// <param name="radius">Radius of the circular path in meters.</param>
        public static double Period(this double speed, double radius)
        {
            if (speed <= 0) throw new ArgumentException("Speed must be greater than zero.");
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return 2 * Math.PI * radius / speed;
        }

        /// <summary>
        /// Computes the frequency of uniform circular motion: f = v / (2πr).
        /// </summary>
        /// <param name="speed">Tangential speed in m/s.</param>
        /// <param name="radius">Radius of the circular path in meters.</param>
        public static double Frequency(this double speed, double radius)
        {
            if (speed <= 0) throw new ArgumentException("Speed must be greater than zero.");
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return speed / (2 * Math.PI * radius);
        }

        /// <summary>
        /// Computes the tangential velocity vector from angular velocity and radius: v = ω × r.
        /// </summary>
        /// <param name="angularVelocity">Angular velocity vector (along rotation axis).</param>
        /// <param name="radiusVector">Vector from center to object.</param>
        public static Vector TangentialVelocity(this Vector angularVelocity, Vector radiusVector)
        {
            return angularVelocity.Cross(radiusVector);
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

        #region Orbital Mechanics

        /// <summary>
        /// Computes gravitational field strength at a distance from a mass: g = GM/r².
        /// </summary>
        /// <param name="centralMass">Mass of the central body in kg.</param>
        /// <param name="distance">Distance from the center of mass in meters.</param>
        public static double GravitationalFieldStrength(this double centralMass, double distance)
        {
            if (distance <= 0) throw new ArgumentException("Distance must be greater than zero.");
            return PhysicsConstants.GravitationalConstant * centralMass / (distance * distance);
        }

        /// <summary>
        /// Computes gravitational force between two masses: F = G·m₁·m₂/r².
        /// </summary>
        /// <param name="mass1">First mass in kg.</param>
        /// <param name="mass2">Second mass in kg.</param>
        /// <param name="distance">Distance between centers in meters.</param>
        public static double GravitationalForce(this double mass1, double mass2, double distance)
        {
            if (distance <= 0) throw new ArgumentException("Distance must be greater than zero.");
            return PhysicsConstants.GravitationalConstant * mass1 * mass2 / (distance * distance);
        }

        /// <summary>
        /// Computes circular orbital speed: v = √(GM/r).
        /// </summary>
        /// <param name="centralMass">Mass of the central body in kg.</param>
        /// <param name="radius">Orbital radius in meters.</param>
        public static double OrbitalSpeed(this double centralMass, double radius)
        {
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return Math.Sqrt(PhysicsConstants.GravitationalConstant * centralMass / radius);
        }

        /// <summary>
        /// Computes the orbital period for a circular orbit: T = 2π√(r³/(GM)).
        /// </summary>
        /// <param name="centralMass">Mass of the central body in kg.</param>
        /// <param name="radius">Orbital radius in meters.</param>
        public static double OrbitalPeriod(this double centralMass, double radius)
        {
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return 2 * Math.PI * Math.Sqrt(radius * radius * radius / (PhysicsConstants.GravitationalConstant * centralMass));
        }

        /// <summary>
        /// Computes escape velocity from a body: v = √(2GM/r).
        /// </summary>
        /// <param name="centralMass">Mass of the body in kg.</param>
        /// <param name="radius">Distance from center of mass in meters.</param>
        public static double EscapeVelocity(this double centralMass, double radius)
        {
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            return Math.Sqrt(2 * PhysicsConstants.GravitationalConstant * centralMass / radius);
        }

        /// <summary>
        /// Computes position on a circular orbit at time t in the XY plane.
        /// r(t) = R·(cos(ωt), sin(ωt), 0) where ω = √(GM/R³).
        /// </summary>
        /// <param name="centralMass">Mass of the central body in kg.</param>
        /// <param name="radius">Orbital radius in meters.</param>
        /// <param name="time">Time elapsed in seconds.</param>
        public static Vector OrbitalPosition(this double centralMass, double radius, double time)
        {
            double omega = Math.Sqrt(PhysicsConstants.GravitationalConstant * centralMass / (radius * radius * radius));
            double angle = omega * time;
            return new Vector(
                radius * Math.Cos(angle),
                radius * Math.Sin(angle),
                0
            );
        }

        /// <summary>
        /// Computes velocity on a circular orbit at time t in the XY plane.
        /// v(t) = Rω·(-sin(ωt), cos(ωt), 0) where ω = √(GM/R³).
        /// </summary>
        /// <param name="centralMass">Mass of the central body in kg.</param>
        /// <param name="radius">Orbital radius in meters.</param>
        /// <param name="time">Time elapsed in seconds.</param>
        public static Vector OrbitalVelocity(this double centralMass, double radius, double time)
        {
            double omega = Math.Sqrt(PhysicsConstants.GravitationalConstant * centralMass / (radius * radius * radius));
            double v = radius * omega;
            double angle = omega * time;
            return new Vector(
                -v * Math.Sin(angle),
                v * Math.Cos(angle),
                0
            );
        }

        /// <summary>
        /// Computes centripetal acceleration on a circular orbit at time t.
        /// a(t) = -ω²R·(cos(ωt), sin(ωt), 0) = -(GM/R²)·r̂, directed toward center.
        /// </summary>
        /// <param name="centralMass">Mass of the central body in kg.</param>
        /// <param name="radius">Orbital radius in meters.</param>
        /// <param name="time">Time elapsed in seconds.</param>
        public static Vector OrbitalAcceleration(this double centralMass, double radius, double time)
        {
            double omega = Math.Sqrt(PhysicsConstants.GravitationalConstant * centralMass / (radius * radius * radius));
            double a = omega * omega * radius; // = GM/R²
            double angle = omega * time;
            return new Vector(
                -a * Math.Cos(angle),
                -a * Math.Sin(angle),
                0
            );
        }

        #endregion

        #region Relative Motion

        /// <summary>
        /// Computes the relative velocity of an object with respect to a reference frame.
        /// v_rel = v_object - v_reference.
        /// </summary>
        /// <param name="objectVelocity">Velocity of the object.</param>
        /// <param name="referenceVelocity">Velocity of the reference frame.</param>
        public static Vector RelativeVelocity(this Vector objectVelocity, Vector referenceVelocity)
        {
            return objectVelocity - referenceVelocity;
        }

        /// <summary>
        /// Computes the relative position of an object with respect to a reference point.
        /// r_rel = r_object - r_reference.
        /// </summary>
        /// <param name="objectPosition">Position of the object.</param>
        /// <param name="referencePosition">Position of the reference point.</param>
        public static Vector RelativePosition(this Vector objectPosition, Vector referencePosition)
        {
            return objectPosition - referencePosition;
        }

        /// <summary>
        /// Computes the relative acceleration of an object with respect to a reference frame.
        /// a_rel = a_object - a_reference.
        /// </summary>
        /// <param name="objectAcceleration">Acceleration of the object.</param>
        /// <param name="referenceAcceleration">Acceleration of the reference frame.</param>
        public static Vector RelativeAcceleration(this Vector objectAcceleration, Vector referenceAcceleration)
        {
            return objectAcceleration - referenceAcceleration;
        }

        /// <summary>
        /// Computes the closing speed (approach/separation rate) between two objects.
        /// Positive value means the objects are approaching; negative means separating.
        /// </summary>
        /// <param name="objectVelocity">Velocity of the object.</param>
        /// <param name="referenceVelocity">Velocity of the reference.</param>
        /// <param name="objectPosition">Position of the object.</param>
        /// <param name="referencePosition">Position of the reference.</param>
        public static double ClosingSpeed(
            this Vector objectVelocity,
            Vector referenceVelocity,
            Vector objectPosition,
            Vector referencePosition)
        {
            var relVel = objectVelocity - referenceVelocity;
            var relPos = objectPosition - referencePosition;
            var distance = relPos.GetMagnitude();
            if (distance < 1e-15) return relVel.GetMagnitude();
            // Project relative velocity onto the line connecting the two objects
            // Negative dot means approaching
            return -(relVel.Dot(relPos)) / distance;
        }

        /// <summary>
        /// Computes the position of an object in a reference frame moving at constant velocity.
        /// r'(t) = (r_object + v_object·t) - (r_reference + v_reference·t).
        /// </summary>
        /// <param name="objectVelocity">Velocity of the object.</param>
        /// <param name="referenceVelocity">Velocity of the reference frame.</param>
        /// <param name="time">Time elapsed in seconds.</param>
        /// <param name="objectInitialPosition">Initial position of the object (default origin).</param>
        /// <param name="referenceInitialPosition">Initial position of the reference frame (default origin).</param>
        public static Vector RelativePositionAtTime(
            this Vector objectVelocity,
            Vector referenceVelocity,
            double time,
            Vector? objectInitialPosition = null,
            Vector? referenceInitialPosition = null)
        {
            var r0Obj = objectInitialPosition ?? new Vector(0, 0, 0);
            var r0Ref = referenceInitialPosition ?? new Vector(0, 0, 0);
            return (r0Obj + time * objectVelocity) - (r0Ref + time * referenceVelocity);
        }

        /// <summary>
        /// Computes the time of closest approach between two objects moving at constant velocity.
        /// Minimizes |r_rel(t)|² = |Δr₀ + Δv·t|².
        /// Returns the time at which distance is minimized (can be negative if closest approach was in the past).
        /// </summary>
        /// <param name="objectVelocity">Velocity of the object.</param>
        /// <param name="referenceVelocity">Velocity of the reference.</param>
        /// <param name="objectPosition">Initial position of the object.</param>
        /// <param name="referencePosition">Initial position of the reference.</param>
        public static double TimeOfClosestApproach(
            this Vector objectVelocity,
            Vector referenceVelocity,
            Vector objectPosition,
            Vector referencePosition)
        {
            var deltaR = objectPosition - referencePosition;
            var deltaV = objectVelocity - referenceVelocity;
            double vv = deltaV.Dot(deltaV);
            if (vv < 1e-15) return 0; // no relative motion
            return -deltaR.Dot(deltaV) / vv;
        }

        /// <summary>
        /// Computes the minimum distance between two objects moving at constant velocity.
        /// </summary>
        /// <param name="objectVelocity">Velocity of the object.</param>
        /// <param name="referenceVelocity">Velocity of the reference.</param>
        /// <param name="objectPosition">Initial position of the object.</param>
        /// <param name="referencePosition">Initial position of the reference.</param>
        public static double MinimumDistance(
            this Vector objectVelocity,
            Vector referenceVelocity,
            Vector objectPosition,
            Vector referencePosition)
        {
            double t = objectVelocity.TimeOfClosestApproach(referenceVelocity, objectPosition, referencePosition);
            var relPos = objectVelocity.RelativePositionAtTime(referenceVelocity, t, objectPosition, referencePosition);
            return relPos.GetMagnitude();
        }

        #endregion
    }
}