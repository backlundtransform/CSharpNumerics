using CSharpNumerics.Constants;
using Numerics;

namespace System;
public static class PhysicsExtensions
{
    /// <summary>
    /// Computes velocity from free fall from height h.
    /// v = sqrt(2 g h)
    /// Height in meters.
    public static double FreeFallVelocity(double height)
    {
        return Math.Sqrt(2 * PhysicsConstants.GravitationalAcceleration * height);
    }
}

