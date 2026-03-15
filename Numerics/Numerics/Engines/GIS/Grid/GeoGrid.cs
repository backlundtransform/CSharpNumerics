using CSharpNumerics.Engines.GIS.Coordinates;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Grid
{
    /// <summary>
    /// A uniform 3-D spatial grid defined by axis-aligned bounds and a cell step size.
    /// Enumerates cell-centre <see cref="Vector"/> positions and maps between
    /// flat indices and (ix, iy, iz) triples (row-major: ix varies fastest).
    /// <para>
    /// Constructed with the syntax:
    /// <c>new GeoGrid(xMin, xMax, yMin, yMax, zMin, zMax, step)</c>
    /// </para>
    /// </summary>
    public class GeoGrid : IEnumerable<Vector>
    {
        /// <summary>Minimum x coordinate (metres).</summary>
        public double XMin { get; }

        /// <summary>Maximum x coordinate (metres).</summary>
        public double XMax { get; }

        /// <summary>Minimum y coordinate (metres).</summary>
        public double YMin { get; }

        /// <summary>Maximum y coordinate (metres).</summary>
        public double YMax { get; }

        /// <summary>Minimum z coordinate (metres).</summary>
        public double ZMin { get; }

        /// <summary>Maximum z coordinate (metres).</summary>
        public double ZMax { get; }

        /// <summary>Uniform cell spacing in all three dimensions (metres).</summary>
        public double Step { get; }

        /// <summary>Number of cells along the x-axis.</summary>
        public int Nx { get; }

        /// <summary>Number of cells along the y-axis.</summary>
        public int Ny { get; }

        /// <summary>Number of cells along the z-axis.</summary>
        public int Nz { get; }

        /// <summary>Total number of cells: Nx × Ny × Nz.</summary>
        public int CellCount => Nx * Ny * Nz;

        /// <summary>
        /// Optional projection that maps local (x, y, z) back to WGS-84 coordinates.
        /// Null when the grid was created with plain metre coordinates.
        /// </summary>
        public Projection Projection { get; private set; }

        /// <summary>
        /// Creates a uniform 3-D grid.
        /// </summary>
        /// <param name="xMin">Minimum x (metres).</param>
        /// <param name="xMax">Maximum x (metres).</param>
        /// <param name="yMin">Minimum y (metres).</param>
        /// <param name="yMax">Maximum y (metres).</param>
        /// <param name="zMin">Minimum z (metres).</param>
        /// <param name="zMax">Maximum z (metres).</param>
        /// <param name="step">Cell spacing in all three dimensions (metres).</param>
        public GeoGrid(double xMin, double xMax,
                        double yMin, double yMax,
                        double zMin, double zMax,
                        double step)
        {
            if (step <= 0) throw new ArgumentException("Step must be positive.", nameof(step));
            if (xMax <= xMin) throw new ArgumentException("xMax must be greater than xMin.");
            if (yMax <= yMin) throw new ArgumentException("yMax must be greater than yMin.");
            if (zMax < zMin) throw new ArgumentException("zMax must be >= zMin.");

            XMin = xMin;
            XMax = xMax;
            YMin = yMin;
            YMax = yMax;
            ZMin = zMin;
            ZMax = zMax;
            Step = step;

            Nx = Math.Max(1, (int)Math.Floor((xMax - xMin) / step) + 1);
            Ny = Math.Max(1, (int)Math.Floor((yMax - yMin) / step) + 1);
            Nz = Math.Max(1, (int)Math.Floor((zMax - zMin) / step) + 1);
        }

        /// <summary>
        /// Creates a geo-referenced grid from a geographic (lat/lon) bounding box.
        /// The bounding box is projected to local metres using the specified projection type.
        /// The south-west corner becomes the projection origin.
        /// </summary>
        /// <param name="southWest">South-west corner (min lat, min lon).</param>
        /// <param name="northEast">North-east corner (max lat, max lon).</param>
        /// <param name="altMin">Minimum altitude in metres.</param>
        /// <param name="altMax">Maximum altitude in metres.</param>
        /// <param name="step">Cell spacing in metres.</param>
        /// <param name="projectionType">Projection method (default: LocalTangentPlane).</param>
        /// <returns>A GeoGrid with an attached <see cref="Coordinates.Projection"/>.</returns>
        public static GeoGrid FromLatLon(
            GeoCoordinate southWest,
            GeoCoordinate northEast,
            double altMin,
            double altMax,
            double step,
            ProjectionType projectionType = ProjectionType.LocalTangentPlane)
        {
            if (northEast.Latitude <= southWest.Latitude)
                throw new ArgumentException("northEast latitude must be greater than southWest latitude.");
            if (northEast.Longitude <= southWest.Longitude)
                throw new ArgumentException("northEast longitude must be greater than southWest longitude.");

            var projection = new Projection(southWest, projectionType);
            var localNE = projection.ToLocal(northEast);

            var grid = new GeoGrid(0, localNE.x, 0, localNE.y, altMin, altMax, step);
            grid.Projection = projection;
            return grid;
        }

        /// <summary>
        /// Returns the geographic coordinate of the cell at the given flat index.
        /// Requires <see cref="Projection"/> to be set (throws otherwise).
        /// </summary>
        public GeoCoordinate CellCentreGeo(int flatIndex)
        {
            if (Projection == null)
                throw new InvalidOperationException("No projection set. Use FromLatLon to create a geo-referenced grid.");
            return Projection.ToGeo(CellCentre(flatIndex));
        }

        /// <summary>
        /// Returns the geographic coordinate of the cell at (ix, iy, iz).
        /// Requires <see cref="Projection"/> to be set (throws otherwise).
        /// </summary>
        public GeoCoordinate CellCentreGeo(int ix, int iy, int iz)
        {
            if (Projection == null)
                throw new InvalidOperationException("No projection set. Use FromLatLon to create a geo-referenced grid.");
            return Projection.ToGeo(CellCentre(ix, iy, iz));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Index ↔ (ix, iy, iz) mapping  (row-major: ix fastest)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Converts 3-D indices to a flat row-major index.
        /// </summary>
        public int Index(int ix, int iy, int iz) => iz * (Nx * Ny) + iy * Nx + ix;

        /// <summary>
        /// Converts a flat index back to 3-D indices.
        /// </summary>
        public (int ix, int iy, int iz) Index3D(int flatIndex)
        {
            int iz = flatIndex / (Nx * Ny);
            int remainder = flatIndex % (Nx * Ny);
            int iy = remainder / Nx;
            int ix = remainder % Nx;
            return (ix, iy, iz);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Position queries
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Returns the world-space centre of cell (ix, iy, iz).
        /// </summary>
        public Vector CellCentre(int ix, int iy, int iz) =>
            new Vector(XMin + ix * Step, YMin + iy * Step, ZMin + iz * Step);

        /// <summary>
        /// Returns the world-space centre of the cell at the given flat index.
        /// </summary>
        public Vector CellCentre(int flatIndex)
        {
            var (ix, iy, iz) = Index3D(flatIndex);
            return CellCentre(ix, iy, iz);
        }

        /// <summary>
        /// Returns the 3-D index of the cell nearest to <paramref name="position"/>.
        /// Clamps to grid bounds.
        /// </summary>
        public (int ix, int iy, int iz) NearestIndex(Vector position)
        {
            int ix = Clamp((int)Math.Round((position.x - XMin) / Step), 0, Nx - 1);
            int iy = Clamp((int)Math.Round((position.y - YMin) / Step), 0, Ny - 1);
            int iz = Clamp((int)Math.Round((position.z - ZMin) / Step), 0, Nz - 1);
            return (ix, iy, iz);
        }

        /// <summary>
        /// Returns the flat index of the cell nearest to <paramref name="position"/>.
        /// </summary>
        public int NearestFlatIndex(Vector position)
        {
            var (ix, iy, iz) = NearestIndex(position);
            return Index(ix, iy, iz);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Enumeration (all cell centres)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Enumerates the centre positions of all cells in flat-index order.
        /// </summary>
        public IEnumerator<Vector> GetEnumerator()
        {
            for (int iz = 0; iz < Nz; iz++)
                for (int iy = 0; iy < Ny; iy++)
                    for (int ix = 0; ix < Nx; ix++)
                        yield return CellCentre(ix, iy, iz);
        }

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

        // ═══════════════════════════════════════════════════════════════
        //  Helpers
        // ═══════════════════════════════════════════════════════════════

        private static int Clamp(int value, int min, int max) =>
            value < min ? min : value > max ? max : value;
    }
}
