using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Numerics.FiniteDifference;

/// <summary>
/// A uniform 3D rectangular grid that maps between a flat <see cref="VectorN"/>
/// (row-major order: ix fastest, then iy, then iz) and 3D index space.
/// Extension of <see cref="Grid2D"/> to three dimensions for volumetric PDE solving.
/// </summary>
public class Grid3D
{
    /// <summary>Number of cells in the x-direction.</summary>
    public int Nx { get; }

    /// <summary>Number of cells in the y-direction.</summary>
    public int Ny { get; }

    /// <summary>Number of cells in the z-direction.</summary>
    public int Nz { get; }

    /// <summary>Cell spacing in the x-direction.</summary>
    public double Dx { get; }

    /// <summary>Cell spacing in the y-direction.</summary>
    public double Dy { get; }

    /// <summary>Cell spacing in the z-direction.</summary>
    public double Dz { get; }

    /// <summary>Total number of cells (Nx × Ny × Nz).</summary>
    public int Length => Nx * Ny * Nz;

    /// <summary>
    /// Creates a uniform 3D grid.
    /// </summary>
    public Grid3D(int nx, int ny, int nz, double dx, double dy, double dz)
    {
        if (nx <= 0) throw new ArgumentException("nx must be positive", nameof(nx));
        if (ny <= 0) throw new ArgumentException("ny must be positive", nameof(ny));
        if (nz <= 0) throw new ArgumentException("nz must be positive", nameof(nz));
        if (dx <= 0) throw new ArgumentException("dx must be positive", nameof(dx));
        if (dy <= 0) throw new ArgumentException("dy must be positive", nameof(dy));
        if (dz <= 0) throw new ArgumentException("dz must be positive", nameof(dz));

        Nx = nx;
        Ny = ny;
        Nz = nz;
        Dx = dx;
        Dy = dy;
        Dz = dz;
    }

    /// <summary>
    /// Creates a uniform 3D grid with equal spacing in all directions.
    /// </summary>
    public Grid3D(int nx, int ny, int nz, double d) : this(nx, ny, nz, d, d, d) { }

    /// <summary>
    /// Converts a 3D index (ix, iy, iz) to a flat row-major index.
    /// Layout: ix varies fastest, then iy, then iz.
    /// </summary>
    public int Index(int ix, int iy, int iz) => iz * (Nx * Ny) + iy * Nx + ix;

    /// <summary>
    /// Converts a flat row-major index to 3D indices (ix, iy, iz).
    /// </summary>
    public (int ix, int iy, int iz) Index3D(int flatIndex)
    {
        int nxy = Nx * Ny;
        int iz = flatIndex / nxy;
        int rem = flatIndex % nxy;
        int iy = rem / Nx;
        int ix = rem % Nx;
        return (ix, iy, iz);
    }

    /// <summary>
    /// Packs a 3D array [ix, iy, iz] into a flat <see cref="VectorN"/> in row-major order.
    /// </summary>
    public VectorN ToVector(double[,,] field)
    {
        if (field.GetLength(0) != Nx || field.GetLength(1) != Ny || field.GetLength(2) != Nz)
            throw new ArgumentException(
                $"Expected [{Nx},{Ny},{Nz}] array, got [{field.GetLength(0)},{field.GetLength(1)},{field.GetLength(2)}]");

        var values = new double[Length];
        for (int iz = 0; iz < Nz; iz++)
            for (int iy = 0; iy < Ny; iy++)
                for (int ix = 0; ix < Nx; ix++)
                    values[Index(ix, iy, iz)] = field[ix, iy, iz];

        return new VectorN(values);
    }

    /// <summary>
    /// Unpacks a flat <see cref="VectorN"/> into a 3D array [ix, iy, iz].
    /// </summary>
    public double[,,] ToArray(VectorN v)
    {
        if (v.Length != Length)
            throw new ArgumentException($"Expected length {Length}, got {v.Length}");

        var field = new double[Nx, Ny, Nz];
        for (int iz = 0; iz < Nz; iz++)
            for (int iy = 0; iy < Ny; iy++)
                for (int ix = 0; ix < Nx; ix++)
                    field[ix, iy, iz] = v[Index(ix, iy, iz)];

        return field;
    }

    /// <summary>
    /// Creates a <see cref="VectorN"/> of zeros with length Nx × Ny × Nz.
    /// </summary>
    public VectorN Zeros() => new VectorN(new double[Length]);

    /// <summary>
    /// Creates a <see cref="VectorN"/> initialised by evaluating a function at each cell centre.
    /// Cell centres are at (ix * Dx, iy * Dy, iz * Dz).
    /// </summary>
    public VectorN Initialize(Func<double, double, double, double> f)
    {
        var values = new double[Length];
        for (int iz = 0; iz < Nz; iz++)
            for (int iy = 0; iy < Ny; iy++)
                for (int ix = 0; ix < Nx; ix++)
                    values[Index(ix, iy, iz)] = f(ix * Dx, iy * Dy, iz * Dz);

        return new VectorN(values);
    }
}
