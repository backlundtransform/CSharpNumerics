using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Coordinates;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace CSharpNumerics.Engines.GIS.Export
{
    /// <summary>
    /// Exports time-animated probability data to CZML (Cesium Markup Language),
    /// a JSON format that Cesium can consume for time-dynamic entity visualisation.
    /// <para>
    /// Each cell is represented as a Point entity with time-varying colour
    /// derived from exceedance probability.
    /// </para>
    /// </summary>
    public static class CesiumExporter
    {
        // ═══════════════════════════════════════════════════════════════
        //  CZML from TimeAnimator
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Generates a CZML document from a <see cref="TimeAnimator"/>.
        /// The document contains one packet per cell with a time-sampled colour
        /// property that maps probability [0..1] to a red→green colour ramp.
        /// </summary>
        /// <param name="animation">Time-animated probability maps.</param>
        /// <param name="name">Document name (shown in Cesium viewer).</param>
        /// <param name="minProbability">
        /// Cells below this probability at all time steps are omitted to
        /// keep file size manageable. Default 0 (include all).
        /// </param>
        public static string ToCzml(
            TimeAnimator animation,
            string name = "Plume Risk",
            double minProbability = 0)
        {
            if (animation == null) throw new ArgumentNullException(nameof(animation));

            var sb = new StringBuilder();
            sb.Append('[');

            // Document packet
            sb.Append("{\"id\":\"document\",\"name\":\"");
            sb.Append(EscapeJson(name));
            sb.Append("\",\"version\":\"1.0\",\"clock\":{\"interval\":\"");
            AppendIso(sb, animation.TimeFrame.Start);
            sb.Append('/');
            AppendIso(sb, animation.TimeFrame.End);
            sb.Append("\",\"currentTime\":\"");
            AppendIso(sb, animation.TimeFrame.Start);
            sb.Append("\",\"multiplier\":1}}");

            // Entity packets — one per cell
            int cellCount = animation.Grid.CellCount;
            var proj = animation.Grid.Projection;
            for (int i = 0; i < cellCount; i++)
            {
                if (!CellExceedsMin(animation, i, minProbability))
                    continue;

                var pos = animation.Grid.CellCentre(i);
                sb.Append(',');
                AppendCellPacket(sb, animation, i, pos, proj);
            }

            sb.Append(']');
            return sb.ToString();
        }

        /// <summary>Write CZML to a file.</summary>
        public static void Save(TimeAnimator animation, string path,
            string name = "Plume Risk", double minProbability = 0)
            => File.WriteAllText(path, ToCzml(animation, name, minProbability), Encoding.UTF8);

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON with Cesium metadata
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Exports a single <see cref="ProbabilityMap"/> as GeoJSON with
        /// Cesium-compatible metadata (3D Tiles hints).
        /// </summary>
        public static string ToGeoJsonCesium(ProbabilityMap map, ExportMetadata metadata = null)
        {
            if (map == null) throw new ArgumentNullException(nameof(map));

            var sb = new StringBuilder();
            sb.Append('{');
            sb.Append("\"type\":\"FeatureCollection\",");

            var proj = map.Grid.Projection;

            // Cesium-specific metadata
            sb.Append("\"metadata\":{");
            if (proj != null)
                sb.Append("\"crs\":\"WGS84\"");
            else
                sb.Append("\"crs\":\"local\"");
            sb.Append(",\"cesium\":{\"viewer\":\"3DTiles\",\"heightReference\":\"CLAMP_TO_GROUND\"}");
            if (metadata?.Simulation != null)
            {
                sb.Append(",\"simulation\":\"");
                sb.Append(EscapeJson(metadata.Simulation));
                sb.Append('"');
            }
            if (metadata != null && !double.IsNaN(metadata.EmissionRate))
            {
                sb.Append(",\"emissionRate\":");
                sb.Append(Fmt(metadata.EmissionRate));
            }
            sb.Append(",\"timeStep\":");
            sb.Append(Fmt(map.Time));
            sb.Append("},");

            sb.Append("\"features\":[");
            bool first = true;
            for (int i = 0; i < map.CellCount; i++)
            {
                var pos = map.Grid.CellCentre(i);
                if (!first) sb.Append(',');
                first = false;

                sb.Append("{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\",\"coordinates\":[");
                if (proj != null)
                {
                    var geo = proj.ToGeo(pos);
                    sb.Append(Fmt(geo.Longitude)); sb.Append(',');
                    sb.Append(Fmt(geo.Latitude)); sb.Append(',');
                    sb.Append(Fmt(geo.Altitude));
                }
                else
                {
                    sb.Append(Fmt(pos.x)); sb.Append(',');
                    sb.Append(Fmt(pos.y)); sb.Append(',');
                    sb.Append(Fmt(pos.z));
                }
                sb.Append("]},\"properties\":{");
                sb.Append("\"probability\":"); sb.Append(Fmt(map[i]));
                sb.Append(",\"threshold\":"); sb.Append(Fmt(map.Threshold));
                sb.Append(",\"timeStep\":"); sb.Append(Fmt(map.Time));
                // Cesium styling hints
                var (r, g, b, a) = ProbabilityToColor(map[i]);
                sb.Append(",\"cesium:color\":[");
                sb.Append(r); sb.Append(',');
                sb.Append(g); sb.Append(',');
                sb.Append(b); sb.Append(',');
                sb.Append(a);
                sb.Append(']');
                sb.Append("}}");
            }

            sb.Append("]}");
            return sb.ToString();
        }

        /// <summary>Write Cesium-flavored GeoJSON to file.</summary>
        public static void SaveGeoJsonCesium(ProbabilityMap map, string path, ExportMetadata metadata = null)
            => File.WriteAllText(path, ToGeoJsonCesium(map, metadata), Encoding.UTF8);

        // ═══════════════════════════════════════════════════════════════
        //  CZML from fire spread snapshots
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Generates a CZML document from a time series of <see cref="SpreadSnapshot"/>.
        /// Each cell is a Point entity with time-varying colour:
        /// red = Burning, grey = Burned, transparent = Unburned.
        /// </summary>
        /// <param name="snapshots">Time-ordered fire spread snapshots.</param>
        /// <param name="timeFrame">The simulation time frame.</param>
        /// <param name="name">Document name shown in Cesium viewer.</param>
        public static string ToFireCzml(
            IReadOnlyList<SpreadSnapshot> snapshots,
            TimeFrame timeFrame,
            string name = "Wildfire Spread")
        {
            if (snapshots == null) throw new ArgumentNullException(nameof(snapshots));
            if (timeFrame == null) throw new ArgumentNullException(nameof(timeFrame));
            if (snapshots.Count == 0) throw new ArgumentException("snapshots must not be empty.");

            var grid = snapshots[0].Grid;
            var proj = grid.Projection;
            int cellCount = grid.CellCount;

            var sb = new StringBuilder();
            sb.Append('[');

            // Document packet
            sb.Append("{\"id\":\"document\",\"name\":\"");
            sb.Append(EscapeJson(name));
            sb.Append("\",\"version\":\"1.0\",\"clock\":{\"interval\":\"");
            AppendIso(sb, timeFrame.Start);
            sb.Append('/');
            AppendIso(sb, timeFrame.End);
            sb.Append("\",\"currentTime\":\"");
            AppendIso(sb, timeFrame.Start);
            sb.Append("\",\"multiplier\":1}}");

            // Entity packets — one per cell that ever burns
            for (int i = 0; i < cellCount; i++)
            {
                if (!CellEverBurns(snapshots, i))
                    continue;

                var pos = grid.CellCentre(i);
                sb.Append(',');
                AppendFireCellPacket(sb, snapshots, timeFrame, i, pos, proj);
            }

            sb.Append(']');
            return sb.ToString();
        }

        /// <summary>Write fire CZML to a file.</summary>
        public static void SaveFireCzml(
            IReadOnlyList<SpreadSnapshot> snapshots,
            TimeFrame timeFrame,
            string path,
            string name = "Wildfire Spread")
            => File.WriteAllText(path, ToFireCzml(snapshots, timeFrame, name), Encoding.UTF8);

        /// <summary>
        /// Maps <see cref="CellBurnState"/> to RGBA (0-255).
        /// Burning → red, Burned → grey, Unburned/Firebreak → transparent.
        /// </summary>
        public static (int r, int g, int b, int a) BurnStateToColor(int burnState)
        {
            switch (burnState)
            {
                case (int)CellBurnState.Burning: return (255, 50, 0, 220);   // red-orange
                case (int)CellBurnState.Burned: return (128, 128, 128, 160); // grey
                case (int)CellBurnState.Firebreak: return (60, 60, 60, 100); // dark grey
                default: return (0, 0, 0, 0);                                // transparent
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Internal helpers
        // ═══════════════════════════════════════════════════════════════

        private static bool CellExceedsMin(TimeAnimator animation, int cellIndex, double minP)
        {
            if (minP <= 0) return true;
            for (int t = 0; t < animation.Count; t++)
                if (animation[t][cellIndex] >= minP) return true;
            return false;
        }

        private static void AppendCellPacket(StringBuilder sb, TimeAnimator anim, int cellIndex, Vector pos, Projection proj)
        {
            sb.Append("{\"id\":\"cell_");
            sb.Append(cellIndex);
            sb.Append("\",\"position\":{\"cartographicDegrees\":[");
            if (proj != null)
            {
                var geo = proj.ToGeo(pos);
                sb.Append(Fmt(geo.Longitude)); sb.Append(',');
                sb.Append(Fmt(geo.Latitude)); sb.Append(',');
                sb.Append(Fmt(geo.Altitude));
            }
            else
            {
                sb.Append(Fmt(pos.x)); sb.Append(',');
                sb.Append(Fmt(pos.y)); sb.Append(',');
                sb.Append(Fmt(pos.z));
            }
            sb.Append("]},\"point\":{\"pixelSize\":6,\"color\":{\"epoch\":\"");
            AppendIso(sb, anim.TimeFrame.Start);
            sb.Append("\",\"rgba\":[");

            for (int t = 0; t < anim.Count; t++)
            {
                if (t > 0) sb.Append(',');
                double seconds = anim.TimeFrame.TimeAt(t) - anim.TimeFrame.Start;
                var (r, g, b, a) = ProbabilityToColor(anim[t][cellIndex]);
                sb.Append(Fmt(seconds)); sb.Append(',');
                sb.Append(r); sb.Append(',');
                sb.Append(g); sb.Append(',');
                sb.Append(b); sb.Append(',');
                sb.Append(a);
            }

            sb.Append("]}}}");
        }

        /// <summary>
        /// Maps probability [0..1] to RGBA (0-255).
        /// 0 → transparent, low → green, high → red.
        /// </summary>
        public static (int r, int g, int b, int a) ProbabilityToColor(double p)
        {
            if (p <= 0) return (0, 0, 0, 0);
            p = Math.Min(p, 1.0);
            int r = (int)(255 * Math.Min(1.0, 2.0 * p));
            int g = (int)(255 * Math.Min(1.0, 2.0 * (1.0 - p)));
            int a = (int)(50 + 205 * p); // 50..255
            return (r, g, 0, a);
        }

        private static void AppendIso(StringBuilder sb, double seconds)
        {
            // Represent as "epoch + seconds" for simplicity
            // In real usage the user would convert to ISO 8601 with a real start time
            sb.Append("2000-01-01T00:00:00Z");
        }

        private static string Fmt(double v) =>
            v.ToString("R", CultureInfo.InvariantCulture);

        private static string EscapeJson(string s) =>
            s.Replace("\\", "\\\\").Replace("\"", "\\\"");

        // ═══════════════════════════════════════════════════════════════
        //  Internal — fire helpers
        // ═══════════════════════════════════════════════════════════════

        private static bool CellEverBurns(IReadOnlyList<SpreadSnapshot> snapshots, int cellIndex)
        {
            for (int t = 0; t < snapshots.Count; t++)
            {
                int state = (int)snapshots[t].Snapshot.GetLayer("burnState")[cellIndex];
                if (state == (int)CellBurnState.Burning || state == (int)CellBurnState.Burned)
                    return true;
            }
            return false;
        }

        private static void AppendFireCellPacket(
            StringBuilder sb,
            IReadOnlyList<SpreadSnapshot> snapshots,
            TimeFrame timeFrame,
            int cellIndex,
            Vector pos,
            Projection proj)
        {
            sb.Append("{\"id\":\"fire_");
            sb.Append(cellIndex);
            sb.Append("\",\"position\":{\"cartographicDegrees\":[");
            if (proj != null)
            {
                var geo = proj.ToGeo(pos);
                sb.Append(Fmt(geo.Longitude)); sb.Append(',');
                sb.Append(Fmt(geo.Latitude)); sb.Append(',');
                sb.Append(Fmt(geo.Altitude));
            }
            else
            {
                sb.Append(Fmt(pos.x)); sb.Append(',');
                sb.Append(Fmt(pos.y)); sb.Append(',');
                sb.Append(Fmt(pos.z));
            }
            sb.Append("]},\"point\":{\"pixelSize\":8,\"color\":{\"epoch\":\"");
            AppendIso(sb, timeFrame.Start);
            sb.Append("\",\"rgba\":[");

            for (int t = 0; t < snapshots.Count; t++)
            {
                if (t > 0) sb.Append(',');
                double seconds = snapshots[t].Time - timeFrame.Start;
                int state = (int)snapshots[t].Snapshot.GetLayer("burnState")[cellIndex];
                var (r, g, b, a) = BurnStateToColor(state);
                sb.Append(Fmt(seconds)); sb.Append(',');
                sb.Append(r); sb.Append(',');
                sb.Append(g); sb.Append(',');
                sb.Append(b); sb.Append(',');
                sb.Append(a);
            }

            sb.Append("]}}}");
        }
    }
}
