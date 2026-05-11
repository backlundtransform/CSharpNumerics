using CSharpNumerics.Engines.Game.Performance;
using CSharpNumerics.Physics.FluidDynamics.Buoyancy;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Real-time 3D Navier-Stokes fluid solver optimized for games.
/// Extension of <see cref="GameFluidSolver2D"/> to three dimensions using the
/// same Stable Fluids method (unconditionally stable semi-Lagrangian advection).
/// 
/// Typical grid sizes: 32³ (fast) to 128³ (detailed).
/// Fields stored as flat arrays indexed by [i + j*nx + k*nx*ny].
/// </summary>
public class GameFluidSolver3D
{
    private readonly int _nx, _ny, _nz, _nxy, _size;
    private readonly FluidConfig _config;

    private double[] _u, _v, _w;
    private double[] _uPrev, _vPrev, _wPrev;
    private double[] _density, _densityPrev;
    private double[] _temperature, _tempPrev;

    // Cached temporary arrays for Project() — avoids per-step allocation
    private double[] _projDiv;
    private double[] _projP;

    private readonly List<FluidEmitter> _emitters = new();
    private readonly List<FluidObstacle> _obstacles = new();

    /// <summary>Read-only density field.</summary>
    public ReadOnlySpan<double> Density => _density;
    /// <summary>Read-only X-velocity.</summary>
    public ReadOnlySpan<double> VelocityX => _u;
    /// <summary>Read-only Y-velocity.</summary>
    public ReadOnlySpan<double> VelocityY => _v;
    /// <summary>Read-only Z-velocity.</summary>
    public ReadOnlySpan<double> VelocityZ => _w;

    /// <summary>Grid X size (including ghost cells).</summary>
    public int NX => _nx;
    /// <summary>Grid Y size.</summary>
    public int NY => _ny;
    /// <summary>Grid Z size.</summary>
    public int NZ => _nz;

    /// <summary>
    /// Creates a 3D game fluid solver.
    /// </summary>
    public GameFluidSolver3D(FluidConfig config)
    {
        _config = config;
        _nx = config.GridX + 2;
        _ny = config.GridY + 2;
        _nz = config.GridZ + 2;
        _nxy = _nx * _ny;
        _size = _nxy * _nz;

        _u = new double[_size];
        _v = new double[_size];
        _w = new double[_size];
        _uPrev = new double[_size];
        _vPrev = new double[_size];
        _wPrev = new double[_size];
        _density = new double[_size];
        _densityPrev = new double[_size];
        _temperature = new double[_size];
        _tempPrev = new double[_size];

        Array.Fill(_temperature, config.AmbientTemperature);
        Array.Fill(_tempPrev, config.AmbientTemperature);

        _projDiv = new double[_size];
        _projP = new double[_size];
    }

    /// <summary>Adds an emitter.</summary>
    public void AddEmitter(FluidEmitter emitter) => _emitters.Add(emitter);
    /// <summary>Removes an emitter.</summary>
    public bool RemoveEmitter(FluidEmitter emitter) => _emitters.Remove(emitter);
    /// <summary>Adds an obstacle.</summary>
    public void AddObstacle(FluidObstacle obstacle) => _obstacles.Add(obstacle);
    /// <summary>Removes an obstacle.</summary>
    public bool RemoveObstacle(FluidObstacle obstacle) => _obstacles.Remove(obstacle);

    /// <summary>Gets density at a grid cell.</summary>
    public double GetDensity(int i, int j, int k) => _density[i + j * _nx + k * _nxy];

    /// <summary>Gets velocity at a grid cell.</summary>
    public (double u, double v, double w) GetVelocity(int i, int j, int k)
    {
        int idx = i + j * _nx + k * _nxy;
        return (_u[idx], _v[idx], _w[idx]);
    }

    /// <summary>
    /// Samples velocity at a world-space position using trilinear interpolation.
    /// </summary>
    public (double u, double v, double w) SampleVelocity(double worldX, double worldY, double worldZ)
    {
        double gx = worldX / _config.CellSize + 0.5;
        double gy = worldY / _config.CellSize + 0.5;
        double gz = worldZ / _config.CellSize + 0.5;
        return (Interpolate(_u, gx, gy, gz), Interpolate(_v, gx, gy, gz), Interpolate(_w, gx, gy, gz));
    }

    /// <summary>
    /// Advances the 3D fluid simulation by one time step.
    /// </summary>
    public void Step(double dt)
    {
        // Emitters
        foreach (var e in _emitters)
        {
            e.Apply3D(_density, _u, _v, _w, _nx, _ny, _nz, dt);
            if (_config.EnableBuoyancy)
                ApplyEmitterTemperature(e, dt);
        }

        // Buoyancy (Z-axis = vertical in 3D)
        if (_config.EnableBuoyancy)
        {
            for (int k = 1; k < _nz - 1; k++)
                for (int j = 1; j < _ny - 1; j++)
                    for (int i = 1; i < _nx - 1; i++)
                    {
                        int idx = i + j * _nx + k * _nxy;
                        _w[idx] = BuoyancyForce.ApplyBuoyancy(
                            _w[idx], _temperature[idx], _config.AmbientTemperature,
                            dt, _config.BuoyancyBeta);
                    }
        }

        // Velocity step
        Swap(ref _u, ref _uPrev); Swap(ref _v, ref _vPrev); Swap(ref _w, ref _wPrev);
        Diffuse(1, _u, _uPrev, _config.Viscosity, dt);
        Diffuse(2, _v, _vPrev, _config.Viscosity, dt);
        Diffuse(3, _w, _wPrev, _config.Viscosity, dt);
        Project(_u, _v, _w);

        Swap(ref _u, ref _uPrev); Swap(ref _v, ref _vPrev); Swap(ref _w, ref _wPrev);
        Advect(1, _u, _uPrev, _uPrev, _vPrev, _wPrev, dt);
        Advect(2, _v, _vPrev, _uPrev, _vPrev, _wPrev, dt);
        Advect(3, _w, _wPrev, _uPrev, _vPrev, _wPrev, dt);
        Project(_u, _v, _w);

        // Vorticity confinement
        if (_config.VorticityConfinementStrength > 0)
            VorticityConfinement.Apply3D(_u, _v, _w, _nx, _ny, _nz,
                _config.VorticityConfinementStrength, dt);

        // Density step
        Swap(ref _density, ref _densityPrev);
        Diffuse(0, _density, _densityPrev, _config.DensityDiffusion, dt);
        Swap(ref _density, ref _densityPrev);
        Advect(0, _density, _densityPrev, _u, _v, _w, dt);

        // Temperature step
        if (_config.EnableBuoyancy)
        {
            Swap(ref _temperature, ref _tempPrev);
            Diffuse(0, _temperature, _tempPrev, _config.DensityDiffusion, dt);
            Swap(ref _temperature, ref _tempPrev);
            Advect(0, _temperature, _tempPrev, _u, _v, _w, dt);
        }

        // Obstacles
        foreach (var obs in _obstacles)
            obs.Apply3D(_u, _v, _w, _nx, _ny, _nz);

        // Density decay (SIMD-accelerated when available)
        SimdMath.Decay(_density, 0.999);
    }

    /// <summary>Resets all fields.</summary>
    public void Reset()
    {
        Array.Clear(_u, 0, _size);
        Array.Clear(_v, 0, _size);
        Array.Clear(_w, 0, _size);
        Array.Clear(_uPrev, 0, _size);
        Array.Clear(_vPrev, 0, _size);
        Array.Clear(_wPrev, 0, _size);
        Array.Clear(_density, 0, _size);
        Array.Clear(_densityPrev, 0, _size);
        Array.Fill(_temperature, _config.AmbientTemperature);
        Array.Fill(_tempPrev, _config.AmbientTemperature);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Core routines (3D extensions)
    // ═══════════════════════════════════════════════════════════════

    private void Diffuse(int b, double[] x, double[] x0, double diff, double dt)
    {
        double a = dt * diff * (_nx - 2) * (_ny - 2);
        int iter = _config.DiffusionIterations;

        for (int l = 0; l < iter; l++)
        {
            for (int k = 1; k < _nz - 1; k++)
                for (int j = 1; j < _ny - 1; j++)
                    for (int i = 1; i < _nx - 1; i++)
                    {
                        int idx = i + j * _nx + k * _nxy;
                        x[idx] = (x0[idx] + a * (
                            x[idx - 1] + x[idx + 1] +
                            x[idx - _nx] + x[idx + _nx] +
                            x[idx - _nxy] + x[idx + _nxy]
                        )) / (1 + 6 * a);
                    }
            SetBoundary(b, x);
        }
    }

    private void Advect(int b, double[] d, double[] d0, double[] u, double[] v, double[] w, double dt)
    {
        double dt0 = dt * (_nx - 2);

        for (int k = 1; k < _nz - 1; k++)
            for (int j = 1; j < _ny - 1; j++)
                for (int i = 1; i < _nx - 1; i++)
                {
                    int idx = i + j * _nx + k * _nxy;
                    double x = i - dt0 * u[idx];
                    double y = j - dt0 * v[idx];
                    double z = k - dt0 * w[idx];

                    x = Math.Clamp(x, 0.5, _nx - 1.5);
                    y = Math.Clamp(y, 0.5, _ny - 1.5);
                    z = Math.Clamp(z, 0.5, _nz - 1.5);

                    d[idx] = Interpolate(d0, x, y, z);
                }
        SetBoundary(b, d);
    }

    private void Project(double[] u, double[] v, double[] w)
    {
        Array.Clear(_projDiv, 0, _size);
        Array.Clear(_projP, 0, _size);
        var div = _projDiv;
        var p = _projP;
        double h = 1.0 / Math.Max(_nx - 2, Math.Max(_ny - 2, _nz - 2));

        for (int k = 1; k < _nz - 1; k++)
            for (int j = 1; j < _ny - 1; j++)
                for (int i = 1; i < _nx - 1; i++)
                {
                    int idx = i + j * _nx + k * _nxy;
                    div[idx] = -0.5 * h * (
                        u[idx + 1] - u[idx - 1] +
                        v[idx + _nx] - v[idx - _nx] +
                        w[idx + _nxy] - w[idx - _nxy]);
                }
        SetBoundary(0, div);
        SetBoundary(0, p);

        for (int l = 0; l < _config.PoissonIterations; l++)
        {
            for (int k = 1; k < _nz - 1; k++)
                for (int j = 1; j < _ny - 1; j++)
                    for (int i = 1; i < _nx - 1; i++)
                    {
                        int idx = i + j * _nx + k * _nxy;
                        p[idx] = (div[idx] +
                            p[idx - 1] + p[idx + 1] +
                            p[idx - _nx] + p[idx + _nx] +
                            p[idx - _nxy] + p[idx + _nxy]) / 6.0;
                    }
            SetBoundary(0, p);
        }

        for (int k = 1; k < _nz - 1; k++)
            for (int j = 1; j < _ny - 1; j++)
                for (int i = 1; i < _nx - 1; i++)
                {
                    int idx = i + j * _nx + k * _nxy;
                    u[idx] -= 0.5 * (p[idx + 1] - p[idx - 1]) / h;
                    v[idx] -= 0.5 * (p[idx + _nx] - p[idx - _nx]) / h;
                    w[idx] -= 0.5 * (p[idx + _nxy] - p[idx - _nxy]) / h;
                }
        SetBoundary(1, u);
        SetBoundary(2, v);
        SetBoundary(3, w);
    }

    private void SetBoundary(int b, double[] x)
    {
        // X-faces
        for (int k = 1; k < _nz - 1; k++)
            for (int j = 1; j < _ny - 1; j++)
            {
                int l = j * _nx + k * _nxy;
                x[l] = b == 1 ? -x[1 + l] : x[1 + l];
                x[(_nx - 1) + l] = b == 1 ? -x[(_nx - 2) + l] : x[(_nx - 2) + l];
            }
        // Y-faces
        for (int k = 1; k < _nz - 1; k++)
            for (int i = 1; i < _nx - 1; i++)
            {
                int l = i + k * _nxy;
                x[l] = b == 2 ? -x[l + _nx] : x[l + _nx];
                x[l + (_ny - 1) * _nx] = b == 2 ? -x[l + (_ny - 2) * _nx] : x[l + (_ny - 2) * _nx];
            }
        // Z-faces
        for (int j = 1; j < _ny - 1; j++)
            for (int i = 1; i < _nx - 1; i++)
            {
                int l = i + j * _nx;
                x[l] = b == 3 ? -x[l + _nxy] : x[l + _nxy];
                x[l + (_nz - 1) * _nxy] = b == 3 ? -x[l + (_nz - 2) * _nxy] : x[l + (_nz - 2) * _nxy];
            }
    }

    private double Interpolate(double[] f, double x, double y, double z)
    {
        int i0 = (int)x, j0 = (int)y, k0 = (int)z;
        int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;
        double sx = x - i0, sy = y - j0, sz = z - k0;

        i0 = Math.Clamp(i0, 0, _nx - 1); i1 = Math.Clamp(i1, 0, _nx - 1);
        j0 = Math.Clamp(j0, 0, _ny - 1); j1 = Math.Clamp(j1, 0, _ny - 1);
        k0 = Math.Clamp(k0, 0, _nz - 1); k1 = Math.Clamp(k1, 0, _nz - 1);

        double c00 = f[i0 + j0 * _nx + k0 * _nxy] * (1 - sx) + f[i1 + j0 * _nx + k0 * _nxy] * sx;
        double c10 = f[i0 + j1 * _nx + k0 * _nxy] * (1 - sx) + f[i1 + j1 * _nx + k0 * _nxy] * sx;
        double c01 = f[i0 + j0 * _nx + k1 * _nxy] * (1 - sx) + f[i1 + j0 * _nx + k1 * _nxy] * sx;
        double c11 = f[i0 + j1 * _nx + k1 * _nxy] * (1 - sx) + f[i1 + j1 * _nx + k1 * _nxy] * sx;

        double c0 = c00 * (1 - sy) + c10 * sy;
        double c1 = c01 * (1 - sy) + c11 * sy;

        return c0 * (1 - sz) + c1 * sz;
    }

    private void ApplyEmitterTemperature(FluidEmitter e, double dt)
    {
        int cx = (int)e.Position.x, cy = (int)e.Position.y, cz = (int)e.Position.z;
        int r = (int)Math.Ceiling(e.Radius);

        for (int k = Math.Max(1, cz - r); k <= Math.Min(_nz - 2, cz + r); k++)
            for (int j = Math.Max(1, cy - r); j <= Math.Min(_ny - 2, cy + r); j++)
                for (int i = Math.Max(1, cx - r); i <= Math.Min(_nx - 2, cx + r); i++)
                {
                    double dx = i - e.Position.x;
                    double dy = j - e.Position.y;
                    double dz = k - e.Position.z;
                    double dist = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                    if (dist > e.Radius) continue;

                    double falloff = 1.0 - dist / e.Radius;
                    int idx = i + j * _nx + k * _nxy;
                    _temperature[idx] += (e.Temperature - _config.AmbientTemperature) * falloff * dt;
                }
    }

    private static void Swap(ref double[] a, ref double[] b)
    {
        var tmp = a; a = b; b = tmp;
    }
}
