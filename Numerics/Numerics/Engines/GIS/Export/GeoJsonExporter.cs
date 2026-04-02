using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Coordinates;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.Wildfire;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace CSharpNumerics.Engines.GIS.Export
{
    /// <summary>
    /// Exports grid data to GeoJSON FeatureCollections.
    /// Each cell becomes a Point feature with concentration, probability,
    /// and optional cluster properties.
    /// <para>
    /// No external JSON library required — builds GeoJSON strings manually.
    /// </para>
    /// </summary>
    public static class GeoJsonExporter
    {
        // ═══════════════════════════════════════════════════════════════
        //  Single snapshot
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports a single <see cref="GridSnapshot"/> (concentration only).
        /// </summary>
        public static string ToGeoJson(GridSnapshot snapshot, ExportMetadata metadata = null)
        {
            if (snapshot == null) throw new ArgumentNullException(nameof(snapshot));
            var proj = snapshot.Grid.Projection;

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, snapshot.Time, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int i = 0; i < snapshot.Count; i++)
            {
                var pos = snapshot.Grid.CellCentre(i);
                double val = snapshot[i];
                if (!first) sb.Append(',');
                first = false;

                sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
                AppendCoord(sb, pos, proj);
                sb.Append("]},\"properties\":{");
                AppendProp(sb, "concentration", val, true);
                AppendProp(sb, "timeStep", snapshot.Time, false);
                AppendSnapshotLayers(sb, snapshot, i);
                sb.Append("}}");
            }

            sb.Append("]}");
            return sb.ToString();
        }

        // ═══════════════════════════════════════════════════════════════
        //  ProbabilityMap
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports a single <see cref="ProbabilityMap"/>.
        /// </summary>
        public static string ToGeoJson(ProbabilityMap map, ExportMetadata metadata = null)
        {
            if (map == null) throw new ArgumentNullException(nameof(map));
            var proj = map.Grid.Projection;

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, map.Time, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int i = 0; i < map.CellCount; i++)
            {
                var pos = map.Grid.CellCentre(i);
                if (!first) sb.Append(',');
                first = false;
                AppendFeature(sb, pos, proj,
                    ("probability", map[i]),
                    ("threshold", map.Threshold),
                    ("timeStep", map.Time));
            }

            sb.Append("]}");
            return sb.ToString();
        }

        // ═══════════════════════════════════════════════════════════════
        //  ProbabilityMap + concentration snapshot + cluster label
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports combined concentration, probability, and cluster data
        /// for one time step.
        /// </summary>
        public static string ToGeoJson(
            GridSnapshot snapshot,
            ProbabilityMap probabilityMap,
            int? clusterLabel = null,
            ExportMetadata metadata = null)
        {
            if (snapshot == null) throw new ArgumentNullException(nameof(snapshot));
            if (probabilityMap == null) throw new ArgumentNullException(nameof(probabilityMap));
            var proj = snapshot.Grid.Projection;

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, snapshot.Time, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int i = 0; i < snapshot.Count; i++)
            {
                var pos = snapshot.Grid.CellCentre(i);
                if (!first) sb.Append(',');
                first = false;

                sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
                AppendCoord(sb, pos, proj);
                sb.Append("]},\"properties\":{");
                AppendProp(sb, "concentration", snapshot[i], true);
                AppendProp(sb, "probability", probabilityMap[i], false);
                AppendProp(sb, "timeStep", snapshot.Time, false);
                if (clusterLabel.HasValue)
                    AppendPropInt(sb, "cluster", clusterLabel.Value, false);
                sb.Append("}}");
            }

            sb.Append("]}");
            return sb.ToString();
        }

        // ═══════════════════════════════════════════════════════════════
        //  TimeAnimator — all time steps in one file
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports all time steps from a <see cref="TimeAnimator"/> as a single
        /// FeatureCollection. Each cell × time-step combination is a separate feature.
        /// </summary>
        public static string ToGeoJson(TimeAnimator animation, ExportMetadata metadata = null)
        {
            if (animation == null) throw new ArgumentNullException(nameof(animation));
            var proj = animation.Grid.Projection;

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, null, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int t = 0; t < animation.Count; t++)
            {
                var map = animation[t];
                for (int i = 0; i < map.CellCount; i++)
                {
                    var pos = map.Grid.CellCentre(i);
                    if (!first) sb.Append(',');
                    first = false;
                    AppendFeature(sb, pos, proj,
                        ("probability", map[i]),
                        ("threshold", map.Threshold),
                        ("timeStep", map.Time),
                        ("timeIndex", t));
                }
            }

            sb.Append("]}");
            return sb.ToString();
        }

        // ═══════════════════════════════════════════════════════════════
        //  ExposurePolygon
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports an <see cref="ExposurePolygon"/> as a GeoJSON Feature
        /// with Polygon geometry.
        /// </summary>
        /// <param name="polygon">The exposure polygon to export.</param>
        /// <param name="projection">
        /// Optional projection for converting local coordinates to WGS-84.
        /// </param>
        public static string ToGeoJson(ExposurePolygon polygon, Projection projection = null)
        {
            if (polygon == null) throw new ArgumentNullException(nameof(polygon));
            var sb = new StringBuilder();
            AppendPolygonFeature(sb, polygon, projection);
            return sb.ToString();
        }

        /// <summary>
        /// Exports multiple <see cref="ExposurePolygon"/> objects as a
        /// GeoJSON FeatureCollection.
        /// </summary>
        public static string ToGeoJson(IList<ExposurePolygon> polygons, Projection projection = null)
        {
            if (polygons == null) throw new ArgumentNullException(nameof(polygons));
            var sb = new StringBuilder();
            sb.Append("{\"type\":\"FeatureCollection\",");
            if (projection != null)
                sb.Append("\"metadata\":{\"crs\":\"WGS84\"},");
            else
                sb.Append("\"metadata\":{\"crs\":\"local\"},");
            sb.Append("\"features\":[");
            for (int i = 0; i < polygons.Count; i++)
            {
                if (i > 0) sb.Append(',');
                AppendPolygonFeature(sb, polygons[i], projection);
            }
            sb.Append("]}");
            return sb.ToString();
        }

        /// <summary>Save an exposure polygon to a GeoJSON file.</summary>
        public static void Save(ExposurePolygon polygon, string path, Projection projection = null)
            => File.WriteAllText(path, ToGeoJson(polygon, projection), Encoding.UTF8);

        private static void AppendPolygonFeature(StringBuilder sb, ExposurePolygon polygon, Projection projection)
        {
            sb.Append("{\"type\":\"Feature\",\"geometry\":");

            if (polygon.Boundary == null || polygon.Boundary.Count < 3)
            {
                sb.Append("null");
            }
            else
            {
                sb.Append("{\"type\":\"Polygon\",\"coordinates\":[[");
                for (int i = 0; i < polygon.Boundary.Count; i++)
                {
                    if (i > 0) sb.Append(',');
                    sb.Append('[');
                    AppendCoord(sb, polygon.Boundary[i], projection);
                    sb.Append(']');
                }
                sb.Append("]]}");
            }

            sb.Append(",\"properties\":{");
            AppendProp(sb, "threshold", polygon.Threshold, true);
            sb.Append(",\"layerName\":\"");
            sb.Append(EscapeJson(polygon.LayerName));
            sb.Append('"');
            sb.Append(",\"exposureType\":\"");
            sb.Append(polygon.Type == ExposureType.Peak ? "peak" : "integrated");
            sb.Append('"');
            AppendProp(sb, "exceedanceCellCount", polygon.ExceedanceCellCount, false);
            AppendProp(sb, "areaSquareMetres", polygon.AreaSquareMetres, false);
            sb.Append("}}");
        }

        // ═══════════════════════════════════════════════════════════════
        //  SpreadSnapshot — fire point features
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports a single <see cref="SpreadSnapshot"/> as point features with
        /// fire-specific properties: burnState, flameLength, rateOfSpread.
        /// </summary>
        public static string ToGeoJson(SpreadSnapshot snapshot, ExportMetadata metadata = null)
        {
            if (snapshot == null) throw new ArgumentNullException(nameof(snapshot));
            var proj = snapshot.Grid.Projection;
            var grid = snapshot.Grid;

            var bs = snapshot.Snapshot.GetLayer("burnState");
            var fl = snapshot.Snapshot.GetLayer("flameLength");
            var ros = snapshot.Snapshot.GetLayer("rateOfSpread");

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, snapshot.Time, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int i = 0; i < grid.CellCount; i++)
            {
                var pos = grid.CellCentre(i);
                if (!first) sb.Append(',');
                first = false;

                sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
                AppendCoord(sb, pos, proj);
                sb.Append("]},\"properties\":{");
                AppendPropInt(sb, "burnState", (int)bs[i], true);
                AppendProp(sb, "flameLength", fl[i], false);
                AppendProp(sb, "rateOfSpread", ros[i], false);
                AppendProp(sb, "timeStep", snapshot.Time, false);
                sb.Append("}}");
            }

            sb.Append("]}");
            return sb.ToString();
        }

        /// <summary>
        /// Exports a time series of <see cref="SpreadSnapshot"/> as a single
        /// FeatureCollection. Each cell × time-step is a separate feature.
        /// </summary>
        public static string ToGeoJson(IReadOnlyList<SpreadSnapshot> snapshots, ExportMetadata metadata = null)
        {
            if (snapshots == null) throw new ArgumentNullException(nameof(snapshots));
            if (snapshots.Count == 0) throw new ArgumentException("snapshots must not be empty.");
            var proj = snapshots[0].Grid.Projection;

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, null, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int t = 0; t < snapshots.Count; t++)
            {
                var snap = snapshots[t];
                var grid = snap.Grid;
                var bs = snap.Snapshot.GetLayer("burnState");
                var fl = snap.Snapshot.GetLayer("flameLength");
                var ros = snap.Snapshot.GetLayer("rateOfSpread");

                for (int i = 0; i < grid.CellCount; i++)
                {
                    var pos = grid.CellCentre(i);
                    if (!first) sb.Append(',');
                    first = false;

                    sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
                    AppendCoord(sb, pos, proj);
                    sb.Append("]},\"properties\":{");
                    AppendPropInt(sb, "burnState", (int)bs[i], true);
                    AppendProp(sb, "flameLength", fl[i], false);
                    AppendProp(sb, "rateOfSpread", ros[i], false);
                    AppendProp(sb, "timeStep", snap.Time, false);
                    AppendPropInt(sb, "timeIndex", t, false);
                    sb.Append("}}");
                }
            }

            sb.Append("]}");
            return sb.ToString();
        }

        /// <summary>
        /// Exports fire perimeters (one polygon per time step) as a
        /// GeoJSON FeatureCollection. Each polygon has a <c>timeIndex</c>
        /// and <c>timeStep</c> property.
        /// </summary>
        public static string FirePerimetersToGeoJson(
            IReadOnlyList<ExposurePolygon> perimeters,
            IReadOnlyList<SpreadSnapshot> snapshots,
            Projection projection = null)
        {
            if (perimeters == null) throw new ArgumentNullException(nameof(perimeters));
            if (snapshots == null) throw new ArgumentNullException(nameof(snapshots));

            var sb = new StringBuilder();
            sb.Append("{\"type\":\"FeatureCollection\",");
            sb.Append("\"metadata\":{");
            sb.Append(projection != null ? "\"crs\":\"WGS84\"" : "\"crs\":\"local\"");
            sb.Append(",\"type\":\"firePerimeters\"");
            sb.Append("},\"features\":[");

            bool first = true;
            for (int i = 0; i < perimeters.Count; i++)
            {
                var poly = perimeters[i];
                if (poly == null || poly.Boundary == null || poly.Boundary.Count < 3)
                    continue;

                if (!first) sb.Append(',');
                first = false;

                sb.Append("{\"type\":\"Feature\",\"geometry\":");
                sb.Append("{\"type\":\"Polygon\",\"coordinates\":[[");
                for (int v = 0; v < poly.Boundary.Count; v++)
                {
                    if (v > 0) sb.Append(',');
                    sb.Append('[');
                    AppendCoord(sb, poly.Boundary[v], projection);
                    sb.Append(']');
                }
                sb.Append("]]},\"properties\":{");
                AppendPropInt(sb, "timeIndex", i, true);
                AppendProp(sb, "timeStep", i < snapshots.Count ? snapshots[i].Time : 0, false);
                AppendProp(sb, "areaSquareMetres", poly.AreaSquareMetres, false);
                sb.Append("}}");
            }

            sb.Append("]}");
            return sb.ToString();
        }

        /// <summary>
        /// Exports per-cell burn probability from a Monte Carlo wildfire result
        /// as a GeoJSON FeatureCollection heatmap.
        /// </summary>
        public static string BurnProbabilityToGeoJson(
            WildfireMonteCarloResult mcResult,
            ExportMetadata metadata = null)
        {
            if (mcResult == null) throw new ArgumentNullException(nameof(mcResult));
            var grid = mcResult.Grid;
            var proj = grid.Projection;

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");
            AppendMetadata(sb, metadata, null, proj);
            sb.Append("\"features\":[");

            bool first = true;
            for (int i = 0; i < grid.CellCount; i++)
            {
                var pos = grid.CellCentre(i);
                if (!first) sb.Append(',');
                first = false;

                sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
                AppendCoord(sb, pos, proj);
                sb.Append("]},\"properties\":{");
                AppendProp(sb, "burnProbability", mcResult.BurnProbability[i], true);
                AppendPropInt(sb, "iterations", mcResult.Iterations, false);
                sb.Append("}}");
            }

            sb.Append("]}");
            return sb.ToString();
        }

        /// <summary>Save a wildfire spread snapshot to a GeoJSON file.</summary>
        public static void Save(SpreadSnapshot snapshot, string path, ExportMetadata metadata = null)
            => File.WriteAllText(path, ToGeoJson(snapshot, metadata), Encoding.UTF8);

        /// <summary>Save a wildfire spread time series to a GeoJSON file.</summary>
        public static void Save(IReadOnlyList<SpreadSnapshot> snapshots, string path, ExportMetadata metadata = null)
            => File.WriteAllText(path, ToGeoJson(snapshots, metadata), Encoding.UTF8);

        /// <summary>Save fire perimeters to a GeoJSON file.</summary>
        public static void SaveFirePerimeters(
            IReadOnlyList<ExposurePolygon> perimeters,
            IReadOnlyList<SpreadSnapshot> snapshots,
            string path,
            Projection projection = null)
            => File.WriteAllText(path, FirePerimetersToGeoJson(perimeters, snapshots, projection), Encoding.UTF8);

        /// <summary>Save burn probability heatmap to a GeoJSON file.</summary>
        public static void SaveBurnProbability(
            WildfireMonteCarloResult mcResult,
            string path,
            ExportMetadata metadata = null)
            => File.WriteAllText(path, BurnProbabilityToGeoJson(mcResult, metadata), Encoding.UTF8);

        // ═══════════════════════════════════════════════════════════════
        //  File writers
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Write a single snapshot to a GeoJSON file.</summary>
        public static void Save(GridSnapshot snapshot, string path, ExportMetadata metadata = null)
            => File.WriteAllText(path, ToGeoJson(snapshot, metadata), Encoding.UTF8);

        /// <summary>Write a probability map to a GeoJSON file.</summary>
        public static void Save(ProbabilityMap map, string path, ExportMetadata metadata = null)
            => File.WriteAllText(path, ToGeoJson(map, metadata), Encoding.UTF8);

        /// <summary>Write a full time animation to a GeoJSON file.</summary>
        public static void Save(TimeAnimator animation, string path, ExportMetadata metadata = null)
            => File.WriteAllText(path, ToGeoJson(animation, metadata), Encoding.UTF8);

        /// <summary>
        /// Write one GeoJSON file per time step.
        /// Files are named <c>{basePath}_t{index}.geojson</c>.
        /// Returns the list of written file paths.
        /// </summary>
        public static List<string> SavePerTimeStep(TimeAnimator animation, string basePath, ExportMetadata metadata = null)
        {
            if (animation == null) throw new ArgumentNullException(nameof(animation));
            var paths = new List<string>();
            for (int t = 0; t < animation.Count; t++)
            {
                string filePath = $"{basePath}_t{t}.geojson";
                Save(animation[t], filePath, metadata);
                paths.Add(filePath);
            }
            return paths;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Internal JSON helpers (no external deps)
        // ═══════════════════════════════════════════════════════════════

        private static void AppendMetadata(StringBuilder sb, ExportMetadata meta, double? time, Projection projection = null)
        {
            sb.Append("\"metadata\":{");
            if (projection != null)
                sb.Append("\"crs\":\"WGS84\"");
            else
                sb.Append("\"crs\":\"local\"");
            if (meta != null)
            {
                if (meta.Simulation != null)
                {
                    sb.Append(",\"simulation\":\"");
                    sb.Append(EscapeJson(meta.Simulation));
                    sb.Append('"');
                }
                if (!double.IsNaN(meta.EmissionRate))
                {
                    sb.Append(",\"emissionRate\":");
                    sb.Append(Fmt(meta.EmissionRate));
                }
                if (meta.Unit != null)
                {
                    sb.Append(",\"unit\":\"");
                    sb.Append(EscapeJson(meta.Unit));
                    sb.Append('"');
                }
            }
            if (time.HasValue)
            {
                sb.Append(",\"timeStep\":");
                sb.Append(Fmt(time.Value));
            }
            sb.Append("},");
        }

        private static void AppendFeature(StringBuilder sb, Vector pos, Projection projection, params (string name, double value)[] props)
        {
            sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
            AppendCoord(sb, pos, projection);
            sb.Append("]},\"properties\":{");
            for (int i = 0; i < props.Length; i++)
            {
                if (i > 0) sb.Append(',');
                sb.Append('"');
                sb.Append(props[i].name);
                sb.Append("\":");
                sb.Append(Fmt(props[i].value));
            }
            sb.Append("}}");
        }

        private static void AppendCoord(StringBuilder sb, Vector pos, Projection projection = null)
        {
            if (projection != null)
            {
                var geo = projection.ToGeo(pos);
                // GeoJSON order: [longitude, latitude, altitude]
                sb.Append(Fmt(geo.Longitude));
                sb.Append(',');
                sb.Append(Fmt(geo.Latitude));
                sb.Append(',');
                sb.Append(Fmt(geo.Altitude));
            }
            else
            {
                sb.Append(Fmt(pos.x));
                sb.Append(',');
                sb.Append(Fmt(pos.y));
                sb.Append(',');
                sb.Append(Fmt(pos.z));
            }
        }

        private static void AppendProp(StringBuilder sb, string name, double value, bool first)
        {
            if (!first) sb.Append(',');
            sb.Append('"');
            sb.Append(name);
            sb.Append("\":");
            sb.Append(Fmt(value));
        }

        private static void AppendPropInt(StringBuilder sb, string name, int value, bool first)
        {
            if (!first) sb.Append(',');
            sb.Append('"');
            sb.Append(name);
            sb.Append("\":");
            sb.Append(value);
        }

        private static void AppendSnapshotLayers(StringBuilder sb, GridSnapshot snapshot, int cellIndex)
        {
            foreach (var layerName in snapshot.LayerNames)
            {
                var layer = snapshot.GetLayer(layerName);
                AppendProp(sb, layerName, layer[cellIndex], false);
            }
        }

        private static string Fmt(double v) =>
            v.ToString("R", CultureInfo.InvariantCulture);

        private static string EscapeJson(string s) =>
            s.Replace("\\", "\\\\").Replace("\"", "\\\"");
    }

    /// <summary>
    /// Optional metadata to include in GeoJSON output.
    /// </summary>
    public class ExportMetadata
    {
        /// <summary>Simulation type (e.g. "GaussianPlume").</summary>
        public string Simulation { get; set; }

        /// <summary>Emission rate Q (kg/s). NaN if not set.</summary>
        public double EmissionRate { get; set; } = double.NaN;

        /// <summary>Concentration unit (e.g. "kg/m³").</summary>
        public string Unit { get; set; }
    }
}
