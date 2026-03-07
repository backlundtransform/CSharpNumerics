using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Numerics.FiniteDifference
{
    /// <summary>
    /// A uniform 2D rectangular grid that maps between a flat <see cref="VectorN"/>
    /// (row-major order) and 2D index space. Designed for Method-of-Lines PDE solving:
    /// pack spatial state into a <see cref="VectorN"/>, apply discrete operators, and
    /// feed the result to the existing ODE solvers (RK4, Euler, Verlet, …).
    /// </summary>
    public class Grid2D
    {
        /// <summary>Number of cells in the x-direction.</summary>
        public int Nx { get; }

        /// <summary>Number of cells in the y-direction.</summary>
        public int Ny { get; }

        /// <summary>Cell spacing in the x-direction.</summary>
        public double Dx { get; }

        /// <summary>Cell spacing in the y-direction.</summary>
        public double Dy { get; }

        /// <summary>Total number of cells (Nx × Ny).</summary>
        public int Length => Nx * Ny;

        /// <summary>
        /// Creates a uniform 2D grid.
        /// </summary>
        /// <param name="nx">Number of cells in x.</param>
        /// <param name="ny">Number of cells in y.</param>
        /// <param name="dx">Cell spacing in x.</param>
        /// <param name="dy">Cell spacing in y.</param>
        public Grid2D(int nx, int ny, double dx, double dy)
        {
            if (nx <= 0) throw new ArgumentException("nx must be positive", nameof(nx));
            if (ny <= 0) throw new ArgumentException("ny must be positive", nameof(ny));
            if (dx <= 0) throw new ArgumentException("dx must be positive", nameof(dx));
            if (dy <= 0) throw new ArgumentException("dy must be positive", nameof(dy));

            Nx = nx;
            Ny = ny;
            Dx = dx;
            Dy = dy;
        }

        /// <summary>
        /// Creates a uniform 2D grid with equal spacing in both directions.
        /// </summary>
        public Grid2D(int nx, int ny, double d) : this(nx, ny, d, d) { }

        /// <summary>
        /// Converts a 2D index (ix, iy) to a flat row-major index.
        /// </summary>
        public int Index(int ix, int iy) => iy * Nx + ix;

        /// <summary>
        /// Converts a flat row-major index to 2D indices (ix, iy).
        /// </summary>
        public (int ix, int iy) Index2D(int flatIndex) => (flatIndex % Nx, flatIndex / Nx);

        /// <summary>
        /// Packs a 2D array [ix, iy] into a flat <see cref="VectorN"/> in row-major order.
        /// </summary>
        public VectorN ToVector(double[,] field)
        {
            if (field.GetLength(0) != Nx || field.GetLength(1) != Ny)
                throw new ArgumentException($"Expected [{Nx},{Ny}] array, got [{field.GetLength(0)},{field.GetLength(1)}]");

            var values = new double[Length];
            for (int iy = 0; iy < Ny; iy++)
                for (int ix = 0; ix < Nx; ix++)
                    values[Index(ix, iy)] = field[ix, iy];

            return new VectorN(values);
        }

        /// <summary>
        /// Unpacks a flat <see cref="VectorN"/> into a 2D array [ix, iy].
        /// </summary>
        public double[,] ToArray(VectorN v)
        {
            if (v.Length != Length)
                throw new ArgumentException($"Expected length {Length}, got {v.Length}");

            var field = new double[Nx, Ny];
            for (int iy = 0; iy < Ny; iy++)
                for (int ix = 0; ix < Nx; ix++)
                    field[ix, iy] = v[Index(ix, iy)];

            return field;
        }

        /// <summary>
        /// Creates a <see cref="VectorN"/> of zeros with length Nx × Ny.
        /// </summary>
        public VectorN Zeros() => new VectorN(new double[Length]);

        /// <summary>
        /// Creates a <see cref="VectorN"/> initialised by evaluating a function at each cell centre.
        /// Cell centres are at (ix * Dx, iy * Dy).
        /// </summary>
        /// <param name="f">Function (x, y) → value.</param>
        public VectorN Initialize(Func<double, double, double> f)
        {
            var values = new double[Length];
            for (int iy = 0; iy < Ny; iy++)
                for (int ix = 0; ix < Nx; ix++)
                    values[Index(ix, iy)] = f(ix * Dx, iy * Dy);

            return new VectorN(values);
        }
    }
}
