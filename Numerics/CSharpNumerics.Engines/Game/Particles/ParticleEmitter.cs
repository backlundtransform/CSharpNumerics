using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Particles;

/// <summary>
/// Defines emission parameters for a particle emitter.
/// Controls spawn rate, direction, spread, and initial particle properties.
/// </summary>
public class ParticleEmitter
{
    /// <summary>World-space position of the emitter.</summary>
    public Vector Position { get; set; } = new Vector(0, 0, 0);

    /// <summary>Emission direction (unit vector).</summary>
    public Vector Direction { get; set; } = new Vector(0, 0, 1);

    /// <summary>Number of particles emitted per second.</summary>
    public double EmissionRate { get; set; } = 50;

    /// <summary>Half-angle of the emission cone in radians.</summary>
    public double ConeAngle { get; set; } = 0.3;

    /// <summary>Initial speed of emitted particles (m/s).</summary>
    public double InitialSpeed { get; set; } = 5.0;

    /// <summary>Random variation in speed [0,1].</summary>
    public double SpeedVariance { get; set; } = 0.2;

    /// <summary>Lifetime of emitted particles in seconds.</summary>
    public double ParticleLifetime { get; set; } = 3.0;

    /// <summary>Random variation in lifetime [0,1].</summary>
    public double LifetimeVariance { get; set; } = 0.1;

    /// <summary>Mass of each emitted particle.</summary>
    public double ParticleMass { get; set; } = 0.01;

    /// <summary>Whether the emitter is currently active.</summary>
    public bool IsActive { get; set; } = true;

    private double _emitAccumulator;

    /// <summary>
    /// Computes how many particles to emit this frame and generates them.
    /// </summary>
    /// <returns>Number of particles to emit.</returns>
    internal int ComputeEmitCount(double dt)
    {
        if (!IsActive) return 0;
        _emitAccumulator += EmissionRate * dt;
        int count = (int)_emitAccumulator;
        _emitAccumulator -= count;
        return count;
    }

    /// <summary>
    /// Generate initial velocity for one particle with randomization.
    /// </summary>
    internal Vector GenerateVelocity(Random rng)
    {
        // Random direction within cone
        double theta = ConeAngle * Math.Sqrt(rng.NextDouble());
        double phi = 2.0 * Math.PI * rng.NextDouble();

        // Create a local frame from Direction
        var dir = Direction;
        double mag = dir.GetMagnitude();
        if (mag < 1e-10) dir = new Vector(0, 0, 1);
        else dir = (1.0 / mag) * dir;

        // Find perpendicular vectors
        var up = Math.Abs(dir.z) < 0.9 ? new Vector(0, 0, 1) : new Vector(1, 0, 0);
        var right = dir.Cross(up);
        double rightMag = right.GetMagnitude();
        if (rightMag > 1e-10) right = (1.0 / rightMag) * right;
        up = right.Cross(dir);

        double sinT = Math.Sin(theta);
        var localDir = Math.Cos(theta) * dir + sinT * Math.Cos(phi) * right + sinT * Math.Sin(phi) * up;

        double speed = InitialSpeed * (1.0 + SpeedVariance * (2.0 * rng.NextDouble() - 1.0));
        return speed * localDir;
    }

    /// <summary>
    /// Generate lifetime for one particle.
    /// </summary>
    internal double GenerateLifetime(Random rng)
    {
        return ParticleLifetime * (1.0 + LifetimeVariance * (2.0 * rng.NextDouble() - 1.0));
    }
}
