using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.FluidDynamics.SPH;

/// <summary>
/// Smoothed Particle Hydrodynamics solver for liquid simulation.
///
/// Particles carry mass, density, pressure, and velocity. Inter-particle
/// forces are computed using smoothing kernels:
///   - Poly6 kernel for density estimation
///   - Spiky kernel gradient for pressure force
///   - Viscosity kernel Laplacian for viscous diffusion
///
/// Based on Müller et al. 2003 "Particle-Based Fluid Simulation for Interactive Applications".
/// </summary>
public class SPHSolver
{
    /// <summary>Particle data.</summary>
    public struct Particle
    {
        public Vector Position;
        public Vector Velocity;
        public Vector Force;
        public double Density;
        public double Pressure;
    }

    private Particle[] _particles;
    private int _count;

    /// <summary>Smoothing radius (support of the kernel).</summary>
    public double SmoothingRadius { get; set; } = 0.1;

    /// <summary>Particle mass.</summary>
    public double ParticleMass { get; set; } = 0.02;

    /// <summary>Rest density of the fluid (kg/m³).</summary>
    public double RestDensity { get; set; } = 1000.0;

    /// <summary>Gas constant for pressure computation (stiffness).</summary>
    public double GasConstant { get; set; } = 2000.0;

    /// <summary>Dynamic viscosity coefficient.</summary>
    public double Viscosity { get; set; } = 1.0;

    /// <summary>External gravity.</summary>
    public Vector Gravity { get; set; } = new Vector(0, 0, -9.81);

    /// <summary>Timestep.</summary>
    public double Dt { get; set; } = 0.001;

    /// <summary>Domain bounds (particles bounce off walls).</summary>
    public Vector BoundsMin { get; set; } = new Vector(-1, -1, -1);

    /// <summary>Domain bounds.</summary>
    public Vector BoundsMax { get; set; } = new Vector(1, 1, 1);

    /// <summary>Coefficient of restitution for boundary collisions.</summary>
    public double BoundaryRestitution { get; set; } = 0.3;

    /// <summary>Number of particles.</summary>
    public int ParticleCount => _count;

    /// <summary>Read-only access to particles.</summary>
    public ReadOnlySpan<Particle> Particles => _particles.AsSpan(0, _count);

    // Precomputed kernel constants
    private double _poly6Coeff;
    private double _spikyGradCoeff;
    private double _viscLapCoeff;

    /// <summary>
    /// Creates an SPH solver with initial particle positions.
    /// </summary>
    public SPHSolver(IReadOnlyList<Vector> initialPositions)
    {
        _count = initialPositions.Count;
        _particles = new Particle[_count];
        for (int i = 0; i < _count; i++)
        {
            _particles[i].Position = initialPositions[i];
            _particles[i].Velocity = new Vector(0, 0, 0);
        }
        UpdateKernelConstants();
    }

    /// <summary>
    /// Add a particle at runtime.
    /// </summary>
    public void AddParticle(Vector position, Vector velocity = default)
    {
        if (_count >= _particles.Length)
            Array.Resize(ref _particles, _particles.Length * 2);

        _particles[_count].Position = position;
        _particles[_count].Velocity = velocity;
        _count++;
    }

    /// <summary>
    /// Advance the simulation by one timestep.
    /// </summary>
    public void Step()
    {
        UpdateKernelConstants();
        ComputeDensityPressure();
        ComputeForces();
        Integrate();
        EnforceBoundary();
    }

    /// <summary>
    /// Advance by multiple steps.
    /// </summary>
    public void Step(int numSteps)
    {
        for (int i = 0; i < numSteps; i++)
            Step();
    }

    private void UpdateKernelConstants()
    {
        double h = SmoothingRadius;
        double h2 = h * h;
        double h9 = h2 * h2 * h2 * h2 * h; // h^9

        // Poly6: W(r,h) = 315/(64πh⁹) * (h²-r²)³
        _poly6Coeff = 315.0 / (64.0 * Math.PI * h9);

        // Spiky gradient: ∇W = -45/(πh⁶) * (h-r)² * r̂
        double h6 = h2 * h2 * h2;
        _spikyGradCoeff = -45.0 / (Math.PI * h6);

        // Viscosity Laplacian: ∇²W = 45/(πh⁶) * (h-r)
        _viscLapCoeff = 45.0 / (Math.PI * h6);
    }

    private void ComputeDensityPressure()
    {
        double h2 = SmoothingRadius * SmoothingRadius;
        double m = ParticleMass;

        for (int i = 0; i < _count; i++)
        {
            double density = 0;
            for (int j = 0; j < _count; j++)
            {
                var diff = _particles[j].Position - _particles[i].Position;
                double r2 = diff.Dot(diff);
                if (r2 < h2)
                {
                    double w = h2 - r2;
                    density += m * _poly6Coeff * w * w * w;
                }
            }

            _particles[i].Density = Math.Max(density, 1e-6);
            _particles[i].Pressure = GasConstant * (_particles[i].Density - RestDensity);
        }
    }

    private void ComputeForces()
    {
        double h = SmoothingRadius;
        double m = ParticleMass;

        for (int i = 0; i < _count; i++)
        {
            var fPressure = new Vector(0, 0, 0);
            var fViscosity = new Vector(0, 0, 0);

            for (int j = 0; j < _count; j++)
            {
                if (i == j) continue;

                var diff = _particles[j].Position - _particles[i].Position;
                double r = diff.GetMagnitude();
                if (r < h && r > 1e-10)
                {
                    var rHat = (1.0 / r) * diff;

                    // Pressure force (spiky kernel gradient)
                    double pAvg = (_particles[i].Pressure + _particles[j].Pressure) / 2.0;
                    double hMinusR = h - r;
                    fPressure = fPressure + (m * pAvg / _particles[j].Density *
                        _spikyGradCoeff * hMinusR * hMinusR) * rHat;

                    // Viscosity force (viscosity kernel Laplacian)
                    var velDiff = _particles[j].Velocity - _particles[i].Velocity;
                    fViscosity = fViscosity + (m / _particles[j].Density *
                        _viscLapCoeff * (h - r)) * velDiff;
                }
            }

            fViscosity = Viscosity * fViscosity;

            // Gravity
            var fGravity = _particles[i].Density * Gravity;

            _particles[i].Force = fPressure + fViscosity + fGravity;
        }
    }

    private void Integrate()
    {
        for (int i = 0; i < _count; i++)
        {
            var accel = (1.0 / _particles[i].Density) * _particles[i].Force;
            _particles[i].Velocity = _particles[i].Velocity + Dt * accel;
            _particles[i].Position = _particles[i].Position + Dt * _particles[i].Velocity;
        }
    }

    private void EnforceBoundary()
    {
        double e = BoundaryRestitution;
        for (int i = 0; i < _count; i++)
        {
            ref var p = ref _particles[i];
            BounceAxis(ref p, 0, BoundsMin.x, BoundsMax.x, e);
            BounceAxis(ref p, 1, BoundsMin.y, BoundsMax.y, e);
            BounceAxis(ref p, 2, BoundsMin.z, BoundsMax.z, e);
        }
    }

    private static void BounceAxis(ref Particle p, int axis, double lo, double hi, double e)
    {
        double pos = axis == 0 ? p.Position.x : axis == 1 ? p.Position.y : p.Position.z;
        double vel = axis == 0 ? p.Velocity.x : axis == 1 ? p.Velocity.y : p.Velocity.z;

        if (pos < lo)
        {
            SetAxis(ref p.Position, axis, lo);
            SetAxis(ref p.Velocity, axis, -vel * e);
        }
        else if (pos > hi)
        {
            SetAxis(ref p.Position, axis, hi);
            SetAxis(ref p.Velocity, axis, -vel * e);
        }
    }

    private static void SetAxis(ref Vector v, int axis, double val)
    {
        v = axis switch
        {
            0 => new Vector(val, v.y, v.z),
            1 => new Vector(v.x, val, v.z),
            _ => new Vector(v.x, v.y, val)
        };
    }
}
