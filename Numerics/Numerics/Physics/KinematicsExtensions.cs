using CSharpNumerics.Physics.Constants;
using System;
using Numerics.Objects;

namespace CSharpNumerics.Physics;
public static class KinematicsExtensions
{
    /// <summary>
    /// Computes velocity from free fall from height h.
    /// v = sqrt(2 g h)
    /// Height in meters.
    public static double FreeFallVelocity(this double height)
    {
        return Math.Sqrt(2 * PhysicsConstants.GravitationalAcceleration * height);
    }
    /// <summary>
    /// Computes velocity assuming constant acceleration:
    /// v = v0 + a t.
    /// Acceleration in m/s², time in seconds.
    /// </summary>
    public static double VelocityFromConstantAcceleration(
        this double acceleration,
        double time,
        double initialVelocity=0)
    {
        return initialVelocity + acceleration* time;
    }

    public static Vector VelocityFromConstantAcceleration(this 
    Vector acceleration,
    double time,
    Vector initialVelocity
    )
    {
        return initialVelocity + time *acceleration;
    }

    /// <summary>
    /// Computes position assuming constant velocity:
    /// x = x0 + v t.
    /// Velocity in m/s, time in seconds.
    /// </summary>
    public static double PositionFromConstantVelocity(
        this double velocity,
        double time,
        double initialPosition = 0.0)
    {
        return initialPosition + velocity * time;
    }
    public static Vector PositionFromConstantVelocity(
        this Vector velocity,
        double time,
        Vector initialPosition)
    {
        return initialPosition + time * velocity;
    }

    /// <summary>
    /// Computes position assuming constant acceleration:
    /// s = s0 + v0 t + ½ a t².
    /// Acceleration in m/s², time in seconds.
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
}

