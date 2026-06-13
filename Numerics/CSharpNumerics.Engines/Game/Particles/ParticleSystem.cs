using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.Particles;

/// <summary>
/// Game particle system supporting emission, lifetime, physics forces,
/// and collision with the world.
///
/// Each particle has position, velocity, age, and lifetime. Forces include
/// gravity, drag, and optional wind from a fluid solver.
///
/// Particles are stored in a contiguous array with swap-remove for O(1) kill.
/// </summary>
public class ParticleSystem
{
    /// <summary>Per-particle data.</summary>
    public struct ParticleData
    {
        public Vector Position;
        public Vector Velocity;
        public double Age;
        public double Lifetime;
        public double Mass;
    }

    private ParticleData[] _particles;
    private int _count;
    private readonly List<ParticleEmitter> _emitters = new();
    private readonly Random _rng;

    /// <summary>Global gravity acceleration.</summary>
    public Vector Gravity { get; set; } = new Vector(0, 0, -9.81);

    /// <summary>Linear drag coefficient. Force = -drag * velocity.</summary>
    public double DragCoefficient { get; set; } = 0.1;

    /// <summary>Global wind velocity (added to drag computation).</summary>
    public Vector Wind { get; set; } = new Vector(0, 0, 0);

    /// <summary>Maximum number of particles.</summary>
    public int MaxParticles { get; }

    /// <summary>Current alive particle count.</summary>
    public int AliveCount => _count;

    /// <summary>Read-only access to alive particles.</summary>
    public ReadOnlySpan<ParticleData> Particles => _particles.AsSpan(0, _count);

    /// <summary>Ground plane height (particles bounce off this Z). Set to NaN to disable.</summary>
    public double GroundHeight { get; set; } = double.NaN;

    /// <summary>Coefficient of restitution for ground bounce.</summary>
    public double GroundRestitution { get; set; } = 0.3;

    /// <summary>
    /// Creates a particle system.
    /// </summary>
    /// <param name="maxParticles">Maximum number of particles.</param>
    /// <param name="seed">Random seed for emission.</param>
    public ParticleSystem(int maxParticles = 10000, int seed = 42)
    {
        MaxParticles = maxParticles;
        _particles = new ParticleData[maxParticles];
        _count = 0;
        _rng = new Random(seed);
    }

    /// <summary>Add an emitter to the system.</summary>
    public void AddEmitter(ParticleEmitter emitter) => _emitters.Add(emitter);

    /// <summary>Remove an emitter.</summary>
    public bool RemoveEmitter(ParticleEmitter emitter) => _emitters.Remove(emitter);

    /// <summary>Number of emitters.</summary>
    public int EmitterCount => _emitters.Count;

    /// <summary>
    /// Spawn a single particle manually.
    /// </summary>
    public void Emit(Vector position, Vector velocity, double lifetime, double mass = 0.01)
    {
        if (_count >= MaxParticles) return;
        _particles[_count] = new ParticleData
        {
            Position = position,
            Velocity = velocity,
            Age = 0,
            Lifetime = lifetime,
            Mass = mass
        };
        _count++;
    }

    /// <summary>
    /// Advance the particle system by dt.
    /// Emits new particles, integrates physics, kills expired particles.
    /// </summary>
    public void Update(double dt)
    {
        // 1. Emit from all emitters
        for (int e = 0; e < _emitters.Count; e++)
        {
            var emitter = _emitters[e];
            int toEmit = emitter.ComputeEmitCount(dt);
            for (int k = 0; k < toEmit; k++)
            {
                if (_count >= MaxParticles) break;
                _particles[_count] = new ParticleData
                {
                    Position = emitter.Position,
                    Velocity = emitter.GenerateVelocity(_rng),
                    Age = 0,
                    Lifetime = emitter.GenerateLifetime(_rng),
                    Mass = emitter.ParticleMass
                };
                _count++;
            }
        }

        // 2. Integrate + age + kill
        int i = 0;
        while (i < _count)
        {
            ref var p = ref _particles[i];
            p.Age += dt;

            if (p.Age >= p.Lifetime)
            {
                // Swap-remove
                _particles[i] = _particles[_count - 1];
                _count--;
                continue;
            }

            // Forces
            var relVel = p.Velocity - Wind;
            var force = p.Mass * Gravity - DragCoefficient * relVel;
            var accel = (1.0 / p.Mass) * force;

            p.Velocity = p.Velocity + dt * accel;
            p.Position = p.Position + dt * p.Velocity;

            // Ground collision
            if (!double.IsNaN(GroundHeight) && p.Position.z < GroundHeight)
            {
                p.Position = new Vector(p.Position.x, p.Position.y, GroundHeight);
                if (p.Velocity.z < 0)
                    p.Velocity = new Vector(p.Velocity.x, p.Velocity.y, -p.Velocity.z * GroundRestitution);
            }

            i++;
        }
    }

    /// <summary>
    /// Apply a wind field from a 2D fluid solver to all particles.
    /// The wind velocity at each particle's XY position is sampled from the grid.
    /// </summary>
    /// <param name="windX">X-velocity grid (flat array).</param>
    /// <param name="windY">Y-velocity grid (flat array).</param>
    /// <param name="gridNx">Grid X size.</param>
    /// <param name="gridNy">Grid Y size.</param>
    /// <param name="gridOriginX">World X of grid origin.</param>
    /// <param name="gridOriginY">World Y of grid origin.</param>
    /// <param name="cellSize">Grid cell size in world units.</param>
    /// <param name="strength">Wind force multiplier.</param>
    /// <param name="dt">Timestep.</param>
    public void ApplyWindField(ReadOnlySpan<double> windX, ReadOnlySpan<double> windY,
        int gridNx, int gridNy, double gridOriginX, double gridOriginY,
        double cellSize, double strength, double dt)
    {
        for (int i = 0; i < _count; i++)
        {
            ref var p = ref _particles[i];

            double gx = (p.Position.x - gridOriginX) / cellSize;
            double gy = (p.Position.y - gridOriginY) / cellSize;

            int ix = (int)gx;
            int iy = (int)gy;
            if (ix < 0 || ix >= gridNx - 1 || iy < 0 || iy >= gridNy - 1) continue;

            double fx = gx - ix;
            double fy = gy - iy;

            int idx00 = ix + iy * gridNx;
            int idx10 = idx00 + 1;
            int idx01 = idx00 + gridNx;
            int idx11 = idx01 + 1;

            double wx = (1 - fx) * (1 - fy) * windX[idx00] + fx * (1 - fy) * windX[idx10]
                       + (1 - fx) * fy * windX[idx01] + fx * fy * windX[idx11];
            double wy = (1 - fx) * (1 - fy) * windY[idx00] + fx * (1 - fy) * windY[idx10]
                       + (1 - fx) * fy * windY[idx01] + fx * fy * windY[idx11];

            var windForce = strength * new Vector(wx, wy, 0);
            p.Velocity = p.Velocity + dt * (1.0 / p.Mass) * windForce;
        }
    }

    /// <summary>
    /// Kill all particles.
    /// </summary>
    public void Clear() => _count = 0;
}
