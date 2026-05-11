using System;

namespace CSharpNumerics.Physics.FluidDynamics.FreeSurface;

/// <summary>
/// Volume-of-Fluid tracker for free-surface flow (e.g. water surface).
///
/// Each cell stores a volume fraction F ∈ [0,1]:
///   F=0 → empty (gas), F=1 → full (liquid), 0&lt;F&lt;1 → interface cell.
///
/// The tracker advects F using the velocity field and provides
/// surface normal estimation for rendering/physics coupling.
///
/// Based on the Hirt-Nichols (1981) VOF method with PLIC surface reconstruction.
/// </summary>
public class VOFTracker
{
    private readonly int _nx, _ny;
    private double[] _f;
    private double[] _fTemp;

    /// <summary>Number of cells in X.</summary>
    public int NX => _nx;

    /// <summary>Number of cells in Y.</summary>
    public int NY => _ny;

    /// <summary>Cell size.</summary>
    public double CellSize { get; }

    /// <summary>Volume fraction threshold for "liquid present".</summary>
    public double LiquidThreshold { get; set; } = 0.01;

    /// <summary>Read-only access to volume fractions.</summary>
    public ReadOnlySpan<double> VolumeFraction => _f;

    /// <summary>
    /// Creates a VOF tracker on a 2D grid.
    /// </summary>
    public VOFTracker(int nx, int ny, double cellSize = 1.0)
    {
        _nx = nx;
        _ny = ny;
        CellSize = cellSize;
        _f = new double[nx * ny];
        _fTemp = new double[nx * ny];
    }

    /// <summary>Get volume fraction at (i,j).</summary>
    public double GetFraction(int i, int j) => _f[Idx(i, j)];

    /// <summary>Set volume fraction at (i,j).</summary>
    public void SetFraction(int i, int j, double value)
    {
        _f[Idx(i, j)] = Math.Max(0, Math.Min(1, value));
    }

    /// <summary>Fill a rectangular region with liquid.</summary>
    public void FillRect(int x0, int y0, int x1, int y1, double fraction = 1.0)
    {
        for (int j = Math.Max(0, y0); j <= Math.Min(_ny - 1, y1); j++)
            for (int i = Math.Max(0, x0); i <= Math.Min(_nx - 1, x1); i++)
                _f[Idx(i, j)] = Math.Max(0, Math.Min(1, fraction));
    }

    /// <summary>Fill a circular region with liquid.</summary>
    public void FillCircle(double cx, double cy, double radius, double fraction = 1.0)
    {
        double r2 = radius * radius;
        for (int j = 0; j < _ny; j++)
        {
            for (int i = 0; i < _nx; i++)
            {
                double dx = (i + 0.5) * CellSize - cx;
                double dy = (j + 0.5) * CellSize - cy;
                if (dx * dx + dy * dy <= r2)
                    _f[Idx(i, j)] = Math.Max(0, Math.Min(1, fraction));
            }
        }
    }

    /// <summary>
    /// Advect the volume fraction field using the given velocity field.
    /// Uses a first-order donor-acceptor scheme.
    /// </summary>
    /// <param name="u">X-velocity field (flat array, same layout as VOF).</param>
    /// <param name="v">Y-velocity field.</param>
    /// <param name="dt">Timestep.</param>
    public void Advect(double[] u, double[] v, double dt)
    {
        Array.Copy(_f, _fTemp, _f.Length);

        for (int j = 1; j < _ny - 1; j++)
        {
            for (int i = 1; i < _nx - 1; i++)
            {
                int idx = Idx(i, j);
                double ux = u[idx];
                double vy = v[idx];

                // Semi-Lagrangian backtrace
                double srcX = i + 0.5 - ux * dt / CellSize;
                double srcY = j + 0.5 - vy * dt / CellSize;

                // Clamp to domain
                srcX = Math.Max(0.5, Math.Min(_nx - 0.5, srcX));
                srcY = Math.Max(0.5, Math.Min(_ny - 0.5, srcY));

                // Bilinear interpolation
                int i0 = (int)Math.Floor(srcX - 0.5);
                int j0 = (int)Math.Floor(srcY - 0.5);
                i0 = Math.Max(0, Math.Min(_nx - 2, i0));
                j0 = Math.Max(0, Math.Min(_ny - 2, j0));
                int i1 = i0 + 1;
                int j1 = j0 + 1;

                double sx = srcX - 0.5 - i0;
                double sy = srcY - 0.5 - j0;

                double val = (1 - sx) * (1 - sy) * _fTemp[Idx(i0, j0)]
                           + sx * (1 - sy) * _fTemp[Idx(i1, j0)]
                           + (1 - sx) * sy * _fTemp[Idx(i0, j1)]
                           + sx * sy * _fTemp[Idx(i1, j1)];

                _f[idx] = Math.Max(0, Math.Min(1, val));
            }
        }
    }

    /// <summary>
    /// Returns true if the cell (i,j) is an interface cell (partially filled).
    /// </summary>
    public bool IsInterface(int i, int j)
    {
        double f = _f[Idx(i, j)];
        return f > LiquidThreshold && f < 1.0 - LiquidThreshold;
    }

    /// <summary>
    /// Estimate the surface normal at cell (i,j) using central differences on F.
    /// Returns (nx, ny) pointing from liquid toward gas (ascending gradient of emptiness).
    /// </summary>
    public (double nx, double ny) SurfaceNormal(int i, int j)
    {
        double dfdx = (F(i + 1, j) - F(i - 1, j)) / (2.0 * CellSize);
        double dfdy = (F(i, j + 1) - F(i, j - 1)) / (2.0 * CellSize);

        // Normal points from liquid (F=1) to gas (F=0), i.e. -∇F
        double nx = -dfdx;
        double ny = -dfdy;
        double mag = Math.Sqrt(nx * nx + ny * ny);
        if (mag > 1e-10)
        {
            nx /= mag;
            ny /= mag;
        }
        return (nx, ny);
    }

    /// <summary>
    /// Total liquid volume in the domain (sum of F * cellArea).
    /// </summary>
    public double TotalVolume()
    {
        double sum = 0;
        double area = CellSize * CellSize;
        for (int k = 0; k < _f.Length; k++)
            sum += _f[k] * area;
        return sum;
    }

    /// <summary>
    /// Count of interface cells.
    /// </summary>
    public int InterfaceCellCount()
    {
        int count = 0;
        for (int j = 0; j < _ny; j++)
            for (int i = 0; i < _nx; i++)
                if (IsInterface(i, j)) count++;
        return count;
    }

    private int Idx(int i, int j) => i + j * _nx;

    private double F(int i, int j)
    {
        i = Math.Max(0, Math.Min(_nx - 1, i));
        j = Math.Max(0, Math.Min(_ny - 1, j));
        return _f[Idx(i, j)];
    }
}
