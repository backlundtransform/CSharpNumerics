using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Analysis
{
    /// <summary>
    /// Generates <see cref="ExposurePolygon"/> boundaries from grid snapshots.
    /// Computes peak or time-integrated values on the ground-level (z = zMin) slice
    /// and extracts the contour boundary using marching squares.
    /// </summary>
    public static class ExposurePolygonGenerator
    {
        /// <summary>
        /// Generates a polygon bounding the area where the peak (maximum across
        /// all time steps) of the specified layer exceeds <paramref name="threshold"/>.
        /// </summary>
        /// <param name="snapshots">Time-series snapshots with the requested layer.</param>
        /// <param name="threshold">Exceedance threshold (in layer units).</param>
        /// <param name="layerName">
        /// Layer to use: <c>null</c> for concentration (default values),
        /// or a named layer such as "activity" or "dose".
        /// </param>
        /// <param name="excludeValues">
        /// Optional set of cell values to exclude from the computation.
        /// Cells matching any of these values are treated as zero.
        /// Useful for filtering out non-data states (e.g. firebreaks).
        /// </param>
        public static ExposurePolygon PeakExposure(
            IList<GridSnapshot> snapshots,
            double threshold,
            string layerName = null,
            ISet<double> excludeValues = null)
        {
            if (snapshots == null || snapshots.Count == 0)
                throw new ArgumentException("Snapshots list must not be empty.", nameof(snapshots));

            var grid = snapshots[0].Grid;
            double[] peakGround = ComputePeakGround(snapshots, grid, layerName, excludeValues);

            return BuildPolygon(grid, peakGround, threshold,
                layerName ?? "concentration", ExposureType.Peak);
        }

        /// <summary>
        /// Generates a polygon bounding the area where the time-integrated
        /// (sum over all time steps × step duration) exposure of the specified
        /// layer exceeds <paramref name="threshold"/>.
        /// </summary>
        /// <param name="snapshots">Time-series snapshots.</param>
        /// <param name="stepSeconds">Duration of each time step in seconds.</param>
        /// <param name="threshold">Exceedance threshold (in layer-units × seconds).</param>
        /// <param name="layerName">
        /// Layer to use: <c>null</c> for concentration, or "activity", "dose", etc.
        /// </param>
        /// <param name="excludeValues">
        /// Optional set of cell values to exclude from the computation.
        /// Cells matching any of these values are treated as zero.
        /// </param>
        public static ExposurePolygon IntegratedExposure(
            IList<GridSnapshot> snapshots,
            double stepSeconds,
            double threshold,
            string layerName = null,
            ISet<double> excludeValues = null)
        {
            if (snapshots == null || snapshots.Count == 0)
                throw new ArgumentException("Snapshots list must not be empty.", nameof(snapshots));

            var grid = snapshots[0].Grid;
            double[] integratedGround = ComputeIntegratedGround(snapshots, grid, stepSeconds, layerName, excludeValues);

            return BuildPolygon(grid, integratedGround, threshold,
                layerName ?? "concentration", ExposureType.Integrated);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Ground-level aggregation
        // ═══════════════════════════════════════════════════════════════

        private static double[] ComputePeakGround(
            IList<GridSnapshot> snapshots, GeoGrid grid, string layerName,
            ISet<double> excludeValues = null)
        {
            int groundCells = grid.Nx * grid.Ny;
            double[] peak = new double[groundCells];

            foreach (var snap in snapshots)
            {
                for (int iy = 0; iy < grid.Ny; iy++)
                {
                    for (int ix = 0; ix < grid.Nx; ix++)
                    {
                        int groundIdx = iy * grid.Nx + ix;
                        int flatIdx = grid.Index(ix, iy, 0); // z = 0 slice
                        double val = GetValue(snap, flatIdx, layerName);
                        if (excludeValues != null && excludeValues.Contains(val))
                            continue;
                        if (val > peak[groundIdx])
                            peak[groundIdx] = val;
                    }
                }
            }
            return peak;
        }

        private static double[] ComputeIntegratedGround(
            IList<GridSnapshot> snapshots, GeoGrid grid,
            double stepSeconds, string layerName,
            ISet<double> excludeValues = null)
        {
            int groundCells = grid.Nx * grid.Ny;
            double[] integrated = new double[groundCells];

            foreach (var snap in snapshots)
            {
                for (int iy = 0; iy < grid.Ny; iy++)
                {
                    for (int ix = 0; ix < grid.Nx; ix++)
                    {
                        int groundIdx = iy * grid.Nx + ix;
                        int flatIdx = grid.Index(ix, iy, 0);
                        double val = GetValue(snap, flatIdx, layerName);
                        if (excludeValues != null && excludeValues.Contains(val))
                            continue;
                        integrated[groundIdx] += val * stepSeconds;
                    }
                }
            }
            return integrated;
        }

        private static double GetValue(GridSnapshot snap, int flatIdx, string layerName)
        {
            if (layerName == null)
                return snap[flatIdx];
            return snap.GetLayer(layerName)[flatIdx];
        }

        // ═══════════════════════════════════════════════════════════════
        //  Polygon construction (marching squares on 2-D grid)
        // ═══════════════════════════════════════════════════════════════

        private static ExposurePolygon BuildPolygon(
            GeoGrid grid, double[] groundValues, double threshold,
            string layerName, ExposureType type)
        {
            int nx = grid.Nx;
            int ny = grid.Ny;
            double step = grid.Step;

            // Count exceedance cells
            int exceedCount = 0;
            for (int i = 0; i < groundValues.Length; i++)
                if (groundValues[i] > threshold) exceedCount++;

            if (exceedCount == 0)
            {
                return new ExposurePolygon(
                    new List<Vector>(), threshold, layerName, type, 0, 0);
            }

            // Build boundary using marching-squares edge segments
            var edgeSegments = MarchingSquares(groundValues, nx, ny, threshold, grid);

            // Chain segments into ordered boundary
            var boundary = ChainSegments(edgeSegments);

            double area = exceedCount * step * step;

            return new ExposurePolygon(boundary, threshold, layerName, type, exceedCount, area);
        }

        /// <summary>
        /// Marching squares: for each 2×2 cell quad, classify corners as
        /// above/below threshold and emit interpolated edge segments.
        /// </summary>
        private static List<(Vector a, Vector b)> MarchingSquares(
            double[] values, int nx, int ny, double threshold, GeoGrid grid)
        {
            var segments = new List<(Vector, Vector)>();

            for (int iy = 0; iy < ny - 1; iy++)
            {
                for (int ix = 0; ix < nx - 1; ix++)
                {
                    // Four corners of the quad (bottom-left origin)
                    //  3---2
                    //  |   |
                    //  0---1
                    double v0 = values[iy * nx + ix];
                    double v1 = values[iy * nx + ix + 1];
                    double v2 = values[(iy + 1) * nx + ix + 1];
                    double v3 = values[(iy + 1) * nx + ix];

                    Vector p0 = GroundPos(grid, ix, iy);
                    Vector p1 = GroundPos(grid, ix + 1, iy);
                    Vector p2 = GroundPos(grid, ix + 1, iy + 1);
                    Vector p3 = GroundPos(grid, ix, iy + 1);

                    int code = 0;
                    if (v0 > threshold) code |= 1;
                    if (v1 > threshold) code |= 2;
                    if (v2 > threshold) code |= 4;
                    if (v3 > threshold) code |= 8;

                    // Skip all-inside or all-outside
                    if (code == 0 || code == 15) continue;

                    // Interpolated edge midpoints
                    Vector bot = Lerp(p0, p1, v0, v1, threshold);   // edge 0-1
                    Vector right = Lerp(p1, p2, v1, v2, threshold); // edge 1-2
                    Vector top = Lerp(p3, p2, v3, v2, threshold);   // edge 3-2
                    Vector left = Lerp(p0, p3, v0, v3, threshold);  // edge 0-3

                    switch (code)
                    {
                        case 1: segments.Add((bot, left)); break;
                        case 2: segments.Add((right, bot)); break;
                        case 3: segments.Add((right, left)); break;
                        case 4: segments.Add((top, right)); break;
                        case 5:
                            // Saddle: use average to disambiguate
                            double avg = (v0 + v1 + v2 + v3) / 4.0;
                            if (avg > threshold)
                            { segments.Add((bot, right)); segments.Add((top, left)); }
                            else
                            { segments.Add((bot, left)); segments.Add((top, right)); }
                            break;
                        case 6: segments.Add((top, bot)); break;
                        case 7: segments.Add((top, left)); break;
                        case 8: segments.Add((left, top)); break;
                        case 9: segments.Add((bot, top)); break;
                        case 10:
                            double avg2 = (v0 + v1 + v2 + v3) / 4.0;
                            if (avg2 > threshold)
                            { segments.Add((left, bot)); segments.Add((right, top)); }
                            else
                            { segments.Add((left, top)); segments.Add((right, bot)); }
                            break;
                        case 11: segments.Add((right, top)); break;
                        case 12: segments.Add((left, right)); break;
                        case 13: segments.Add((bot, right)); break;
                        case 14: segments.Add((left, bot)); break;
                    }
                }
            }

            return segments;
        }

        private static Vector GroundPos(GeoGrid grid, int ix, int iy)
        {
            return new Vector(grid.XMin + ix * grid.Step, grid.YMin + iy * grid.Step, grid.ZMin);
        }

        private static Vector Lerp(Vector a, Vector b, double va, double vb, double threshold)
        {
            if (Math.Abs(va - vb) < 1e-15)
                return new Vector((a.x + b.x) / 2, (a.y + b.y) / 2, a.z);

            double t = (threshold - va) / (vb - va);
            t = t < 0 ? 0 : (t > 1 ? 1 : t);
            return new Vector(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                a.z);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Segment chaining
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Chains edge segments into a closed boundary ring by matching
        /// endpoints within a small tolerance.
        /// Falls back to convex-hull of all segment endpoints when
        /// chaining does not produce a clean ring.
        /// </summary>
        private static List<Vector> ChainSegments(List<(Vector a, Vector b)> segments)
        {
            if (segments.Count == 0)
                return new List<Vector>();

            double tol = 1e-6;
            var remaining = new List<(Vector a, Vector b)>(segments);
            var ring = new List<Vector>();

            // Start with the first segment
            ring.Add(remaining[0].a);
            ring.Add(remaining[0].b);
            remaining.RemoveAt(0);

            int maxIter = remaining.Count + 1;
            int iter = 0;
            while (remaining.Count > 0 && iter < maxIter)
            {
                iter++;
                var tail = ring[ring.Count - 1];
                bool found = false;

                for (int i = 0; i < remaining.Count; i++)
                {
                    if (DistSq(tail, remaining[i].a) < tol)
                    {
                        ring.Add(remaining[i].b);
                        remaining.RemoveAt(i);
                        found = true;
                        break;
                    }
                    if (DistSq(tail, remaining[i].b) < tol)
                    {
                        ring.Add(remaining[i].a);
                        remaining.RemoveAt(i);
                        found = true;
                        break;
                    }
                }

                if (!found) break;
            }

            // Close the ring if not already closed
            if (ring.Count > 1 && DistSq(ring[0], ring[ring.Count - 1]) > tol)
                ring.Add(ring[0]);

            // If chaining failed (too many segments left), fall back to convex hull
            if (remaining.Count > segments.Count / 3 && segments.Count > 4)
            {
                return ConvexHull(segments);
            }

            return ring;
        }

        /// <summary>
        /// Convex hull of all segment endpoints (Graham scan).
        /// Used as fallback when line chaining doesn't produce a clean ring.
        /// </summary>
        private static List<Vector> ConvexHull(List<(Vector a, Vector b)> segments)
        {
            var points = new List<Vector>(segments.Count * 2);
            foreach (var seg in segments)
            {
                points.Add(seg.a);
                points.Add(seg.b);
            }

            if (points.Count < 3) return points;

            // Find lowest-leftmost point
            int pivotIdx = 0;
            for (int i = 1; i < points.Count; i++)
            {
                if (points[i].y < points[pivotIdx].y ||
                    (points[i].y == points[pivotIdx].y && points[i].x < points[pivotIdx].x))
                    pivotIdx = i;
            }
            var pivot = points[pivotIdx];
            points.RemoveAt(pivotIdx);

            // Sort by polar angle relative to pivot
            points.Sort((a, b) =>
            {
                double cross = Cross2D(
                    a.x - pivot.x, a.y - pivot.y,
                    b.x - pivot.x, b.y - pivot.y);
                if (Math.Abs(cross) < 1e-12)
                {
                    double da = (a.x - pivot.x) * (a.x - pivot.x) + (a.y - pivot.y) * (a.y - pivot.y);
                    double db = (b.x - pivot.x) * (b.x - pivot.x) + (b.y - pivot.y) * (b.y - pivot.y);
                    return da.CompareTo(db);
                }
                return cross > 0 ? -1 : 1;
            });

            var hull = new List<Vector> { pivot };
            foreach (var pt in points)
            {
                while (hull.Count >= 2)
                {
                    var top = hull[hull.Count - 1];
                    var next = hull[hull.Count - 2];
                    double cross = Cross2D(
                        top.x - next.x, top.y - next.y,
                        pt.x - top.x, pt.y - top.y);
                    if (cross <= 0)
                        hull.RemoveAt(hull.Count - 1);
                    else break;
                }
                hull.Add(pt);
            }

            // Close
            if (hull.Count > 1)
                hull.Add(hull[0]);
            return hull;
        }

        private static double Cross2D(double ax, double ay, double bx, double by) =>
            ax * by - ay * bx;

        private static double DistSq(Vector a, Vector b)
        {
            double dx = a.x - b.x, dy = a.y - b.y;
            return dx * dx + dy * dy;
        }
    }
}
