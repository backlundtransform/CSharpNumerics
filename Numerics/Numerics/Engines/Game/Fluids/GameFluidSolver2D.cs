using CSharpNumerics.Engines.Game.Performance;
using CSharpNumerics.Physics.FluidDynamics.Buoyancy;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.Fluids;

/// <summary>
/// Real-time 2D Navier-Stokes fluid solver optimized for games.
/// Based on Jos Stam's "Stable Fluids" method — unconditionally stable at any time step.
/// 
/// Pipeline per step:
///   1. Add forces (emitters, buoyancy)
///   2. Diffuse velocity
///   3. Project (pressure solve → divergence-free)
///   4. Advect velocity (semi-Lagrangian)
///   5. Project again
///   6. Diffuse density
///   7. Advect density
///   8. Vorticity confinement (optional)
///   9. Apply obstacles
/// 
/// Fields are stored as flat arrays indexed by [i + j * nx].
/// Boundary cells (i=0, i=nx-1, j=0, j=ny-1) are ghost cells.
/// </summary>
public class GameFluidSolver2D
{
    private readonly int _nx, _ny, _size;
    private readonly FluidConfig _config;

    // Velocity fields
    private double[] _u, _v;
    private double[] _uPrev, _vPrev;

    // Scalar fields
    private double[] _density;
    private double[] _densityPrev;
    private double[] _temperature;
    private double[] _tempPrev;

    // Cached temporary arrays for Project() — avoids per-step allocation
    private double[] _projDiv;
    private double[] _projP;

    // Emitters and obstacles
    private readonly List<FluidEmitter> _emitters = new();
    private readonly List<FluidObstacle> _obstacles = new();

    /// <summary>Read-only access to the density field.</summary>
    public ReadOnlySpan<double> Density => _density;

    /// <summary>Read-only access to the X-velocity field.</summary>
    public ReadOnlySpan<double> VelocityX => _u;

    /// <summary>Read-only access to the Y-velocity field.</summary>
    public ReadOnlySpan<double> VelocityY => _v;

    /// <summary>Read-only access to the temperature field.</summary>
    public ReadOnlySpan<double> Temperature => _temperature;

    /// <summary>Grid width (including ghost cells).</summary>
    public int NX => _nx;

    /// <summary>Grid height (including ghost cells).</summary>
    public int NY => _ny;

    /// <summary>
    /// Creates a 2D game fluid solver.
    /// </summary>
    /// <param name="config">Fluid configuration.</param>
    public GameFluidSolver2D(FluidConfig config)
    {
        _config = config;
        _nx = config.GridX + 2; // +2 for ghost cells
        _ny = config.GridY + 2;
        _size = _nx * _ny;

        _u = new double[_size];
        _v = new double[_size];
        _uPrev = new double[_size];
        _vPrev = new double[_size];
        _density = new double[_size];
        _densityPrev = new double[_size];
        _temperature = new double[_size];
        _tempPrev = new double[_size];

        // Initialize temperature to ambient
        Array.Fill(_temperature, config.AmbientTemperature);
        Array.Fill(_tempPrev, config.AmbientTemperature);

        _projDiv = new double[_size];
        _projP = new double[_size];
    }

    /// <summary>Adds an emitter to the simulation.</summary>
    public void AddEmitter(FluidEmitter emitter) => _emitters.Add(emitter);

    /// <summary>Removes an emitter.</summary>
    public bool RemoveEmitter(FluidEmitter emitter) => _emitters.Remove(emitter);

    /// <summary>Removes all emitters.</summary>
    public void ClearEmitters() => _emitters.Clear();

    /// <summary>Adds an obstacle to the simulation.</summary>
    public void AddObstacle(FluidObstacle obstacle) => _obstacles.Add(obstacle);

    /// <summary>Removes an obstacle.</summary>
    public bool RemoveObstacle(FluidObstacle obstacle) => _obstacles.Remove(obstacle);

    /// <summary>Removes all obstacles.</summary>
    public void ClearObstacles() => _obstacles.Clear();

    /// <summary>
    /// Gets the density at a grid cell.
    /// </summary>
    public double GetDensity(int i, int j) => _density[i + j * _nx];

    /// <summary>
    /// Gets velocity at a grid cell.
    /// </summary>
    public (double u, double v) GetVelocity(int i, int j)
    {
        int idx = i + j * _nx;
        return (_u[idx], _v[idx]);
    }

    /// <summary>
    /// Samples the velocity at a world-space position using bilinear interpolation.
    /// </summary>
    /// <param name="worldX">X position in metres.</param>
    /// <param name="worldY">Y position in metres.</param>
    /// <returns>Interpolated velocity (u, v) in m/s.</returns>
    public (double u, double v) SampleVelocity(double worldX, double worldY)
    {
        double gx = worldX / _config.CellSize + 0.5;
        double gy = worldY / _config.CellSize + 0.5;
        return (Interpolate(_u, gx, gy), Interpolate(_v, gx, gy));
    }

    /// <summary>
    /// Advances the fluid simulation by one time step.
    /// </summary>
    /// <param name="dt">Time step in seconds.</param>
    public void Step(double dt)
    {
        // Apply emitters
        foreach (var e in _emitters)
        {
            e.Apply2D(_density, _u, _v, _nx, _ny, dt);
            if (_config.EnableBuoyancy)
            {
                // Emitters also inject temperature
                ApplyEmitterTemperature(e, dt);
            }
        }

        // Buoyancy: hot gas rises (applied to v-field as vertical acceleration)
        if (_config.EnableBuoyancy)
        {
            for (int j = 1; j < _ny - 1; j++)
            {
                for (int i = 1; i < _nx - 1; i++)
                {
                    int idx = i + j * _nx;
                    _v[idx] = BuoyancyForce.ApplyBuoyancy(
                        _v[idx], _temperature[idx], _config.AmbientTemperature,
                        dt, _config.BuoyancyBeta);
                }
            }
        }

        // --- Velocity step ---
        // Diffuse
        Swap(ref _u, ref _uPrev);
        Swap(ref _v, ref _vPrev);
        Diffuse(1, _u, _uPrev, _config.Viscosity, dt);
        Diffuse(2, _v, _vPrev, _config.Viscosity, dt);

        // Project to make divergence-free
        Project(_u, _v);

        // Advect
        Swap(ref _u, ref _uPrev);
        Swap(ref _v, ref _vPrev);
        Advect(1, _u, _uPrev, _uPrev, _vPrev, dt);
        Advect(2, _v, _vPrev, _uPrev, _vPrev, dt);

        // Project again
        Project(_u, _v);

        // Vorticity confinement
        if (_config.VorticityConfinementStrength > 0)
            VorticityConfinement.Apply2D(_u, _v, _nx, _ny, _config.VorticityConfinementStrength, dt);

        // --- Density step ---
        Swap(ref _density, ref _densityPrev);
        Diffuse(0, _density, _densityPrev, _config.DensityDiffusion, dt);
        Swap(ref _density, ref _densityPrev);
        Advect(0, _density, _densityPrev, _u, _v, dt);

        // --- Temperature step (if buoyancy enabled) ---
        if (_config.EnableBuoyancy)
        {
            Swap(ref _temperature, ref _tempPrev);
            Diffuse(0, _temperature, _tempPrev, _config.DensityDiffusion, dt);
            Swap(ref _temperature, ref _tempPrev);
            Advect(0, _temperature, _tempPrev, _u, _v, dt);
        }

        // Apply obstacles
        foreach (var obs in _obstacles)
            obs.Apply2D(_u, _v, _nx, _ny);

        // Density decay (SIMD-accelerated when available)
        SimdMath.Decay(_density, 0.999);
    }

    /// <summary>Resets all fields to zero.</summary>
    public void Reset()
    {
        Array.Clear(_u, 0, _size);
        Array.Clear(_v, 0, _size);
        Array.Clear(_uPrev, 0, _size);
        Array.Clear(_vPrev, 0, _size);
        Array.Clear(_density, 0, _size);
        Array.Clear(_densityPrev, 0, _size);
        Array.Fill(_temperature, _config.AmbientTemperature);
        Array.Fill(_tempPrev, _config.AmbientTemperature);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Core Stable Fluids routines
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Gauss-Seidel relaxation for implicit diffusion: (I − dt·ν·∇²)x = x0.
    /// </summary>
    private void Diffuse(int b, double[] x, double[] x0, double diff, double dt)
    {
        double a = dt * diff * (_nx - 2) * (_ny - 2);
        int iter = _config.DiffusionIterations;

        for (int k = 0; k < iter; k++)
        {
            for (int j = 1; j < _ny - 1; j++)
            {
                for (int i = 1; i < _nx - 1; i++)
                {
                    int idx = i + j * _nx;
                    x[idx] = (x0[idx] + a * (
                        x[idx - 1] + x[idx + 1] +
                        x[idx - _nx] + x[idx + _nx]
                    )) / (1 + 4 * a);
                }
            }
            SetBoundary(b, x);
        }
    }

    /// <summary>
    /// Semi-Lagrangian advection: trace particle back in time, interpolate.
    /// </summary>
    private void Advect(int b, double[] d, double[] d0, double[] u, double[] v, double dt)
    {
        double dt0x = dt * (_nx - 2);
        double dt0y = dt * (_ny - 2);

        for (int j = 1; j < _ny - 1; j++)
        {
            for (int i = 1; i < _nx - 1; i++)
            {
                int idx = i + j * _nx;
                double x = i - dt0x * u[idx];
                double y = j - dt0y * v[idx];

                x = Math.Clamp(x, 0.5, _nx - 1.5);
                y = Math.Clamp(y, 0.5, _ny - 1.5);

                d[idx] = Interpolate(d0, x, y);
            }
        }
        SetBoundary(b, d);
    }

    /// <summary>
    /// Helmholtz-Hodge decomposition: subtract pressure gradient to enforce ∇·v = 0.
    /// </summary>
    private void Project(double[] u, double[] v)
    {
        Array.Clear(_projDiv, 0, _size);
        Array.Clear(_projP, 0, _size);
        var div = _projDiv;
        var p = _projP;

        double h = 1.0 / Math.Max(_nx - 2, _ny - 2);

        // Compute divergence
        for (int j = 1; j < _ny - 1; j++)
        {
            for (int i = 1; i < _nx - 1; i++)
            {
                int idx = i + j * _nx;
                div[idx] = -0.5 * h * (
                    u[idx + 1] - u[idx - 1] +
                    v[idx + _nx] - v[idx - _nx]);
            }
        }
        SetBoundary(0, div);
        SetBoundary(0, p);

        // Solve Poisson equation: ∇²p = div
        for (int k = 0; k < _config.PoissonIterations; k++)
        {
            for (int j = 1; j < _ny - 1; j++)
            {
                for (int i = 1; i < _nx - 1; i++)
                {
                    int idx = i + j * _nx;
                    p[idx] = (div[idx] +
                        p[idx - 1] + p[idx + 1] +
                        p[idx - _nx] + p[idx + _nx]) / 4.0;
                }
            }
            SetBoundary(0, p);
        }

        // Subtract pressure gradient
        for (int j = 1; j < _ny - 1; j++)
        {
            for (int i = 1; i < _nx - 1; i++)
            {
                int idx = i + j * _nx;
                u[idx] -= 0.5 * (p[idx + 1] - p[idx - 1]) / h;
                v[idx] -= 0.5 * (p[idx + _nx] - p[idx - _nx]) / h;
            }
        }
        SetBoundary(1, u);
        SetBoundary(2, v);
    }

    /// <summary>
    /// Sets boundary conditions on ghost cells.
    /// b=0: scalar (copy), b=1: u-velocity (negate at x-walls), b=2: v-velocity (negate at y-walls).
    /// </summary>
    private void SetBoundary(int b, double[] x)
    {
        for (int i = 1; i < _nx - 1; i++)
        {
            x[i] = b == 2 ? -x[i + _nx] : x[i + _nx];                               // bottom
            x[i + (_ny - 1) * _nx] = b == 2 ? -x[i + (_ny - 2) * _nx] : x[i + (_ny - 2) * _nx]; // top
        }
        for (int j = 1; j < _ny - 1; j++)
        {
            x[j * _nx] = b == 1 ? -x[1 + j * _nx] : x[1 + j * _nx];                 // left
            x[(_nx - 1) + j * _nx] = b == 1 ? -x[(_nx - 2) + j * _nx] : x[(_nx - 2) + j * _nx]; // right
        }

        // Corners: average of neighbours
        x[0] = 0.5 * (x[1] + x[_nx]);
        x[_nx - 1] = 0.5 * (x[_nx - 2] + x[2 * _nx - 1]);
        x[(_ny - 1) * _nx] = 0.5 * (x[(_ny - 2) * _nx] + x[(_ny - 1) * _nx + 1]);
        x[_ny * _nx - 1] = 0.5 * (x[_ny * _nx - 2] + x[(_ny - 1) * _nx - 1]);
    }

    private double Interpolate(double[] field, double x, double y)
    {
        int i0 = (int)x, j0 = (int)y;
        int i1 = i0 + 1, j1 = j0 + 1;
        double sx = x - i0, sy = y - j0;

        i0 = Math.Clamp(i0, 0, _nx - 1);
        i1 = Math.Clamp(i1, 0, _nx - 1);
        j0 = Math.Clamp(j0, 0, _ny - 1);
        j1 = Math.Clamp(j1, 0, _ny - 1);

        return (1 - sx) * ((1 - sy) * field[i0 + j0 * _nx] + sy * field[i0 + j1 * _nx])
             + sx * ((1 - sy) * field[i1 + j0 * _nx] + sy * field[i1 + j1 * _nx]);
    }

    private void ApplyEmitterTemperature(FluidEmitter e, double dt)
    {
        int cx = (int)e.Position.x;
        int cy = (int)e.Position.y;
        int r = (int)Math.Ceiling(e.Radius);

        for (int j = Math.Max(1, cy - r); j <= Math.Min(_ny - 2, cy + r); j++)
        {
            for (int i = Math.Max(1, cx - r); i <= Math.Min(_nx - 2, cx + r); i++)
            {
                double dx = i - e.Position.x;
                double dy = j - e.Position.y;
                double dist = Math.Sqrt(dx * dx + dy * dy);
                if (dist > e.Radius) continue;

                double falloff = 1.0 - dist / e.Radius;
                int idx = i + j * _nx;
                _temperature[idx] += (e.Temperature - _config.AmbientTemperature) * falloff * dt;
            }
        }
    }

    private static void Swap(ref double[] a, ref double[] b)
    {
        var tmp = a; a = b; b = tmp;
    }
}
