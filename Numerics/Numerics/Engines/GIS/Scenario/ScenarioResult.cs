using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.GIS.Scenario
{
    /// <summary>
    /// Terminal result of the fluent risk-scenario pipeline.
    /// Holds all layers of output: deterministic snapshots, Monte Carlo data,
    /// clustering analysis, probability maps, and time-animated probabilities.
    /// </summary>
    public class ScenarioResult
    {
        // ── Deterministic (RunSingle) ────────────────────────────────

        /// <summary>
        /// Deterministic snapshots (one per time step).
        /// Only set when <see cref="RiskScenarioBuilder.RunSingle"/> was used.
        /// </summary>
        public List<GridSnapshot> Snapshots { get; }

        // ── Monte Carlo ──────────────────────────────────────────────

        /// <summary>
        /// Full Monte Carlo output. Null for deterministic runs.
        /// </summary>
        public MonteCarloScenarioResult MonteCarloResult { get; }

        // ── Clustering ───────────────────────────────────────────────

        /// <summary>
        /// Clustering analysis result. Null when clustering was not applied.
        /// </summary>
        public ClusterAnalysisResult ClusterAnalysis { get; }

        // ── Probability / time animation ─────────────────────────────

        /// <summary>
        /// Time-animated probability maps. Null for deterministic runs.
        /// When clustering was applied, probabilities are filtered to the
        /// dominant cluster.
        /// </summary>
        public TimeAnimator Animation { get; }

        /// <summary>Concentration threshold used for probability maps.</summary>
        public double Threshold { get; }

        // ── Grid & time ──────────────────────────────────────────────

        /// <summary>The spatial grid.</summary>
        public GeoGrid Grid { get; }

        /// <summary>The time frame.</summary>
        public TimeFrame TimeFrame { get; }

        // ═══════════════════════════════════════════════════════════════
        //  Constructors
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Monte Carlo constructor (with or without clustering).</summary>
        internal ScenarioResult(
            MonteCarloScenarioResult mcResult,
            ClusterAnalysisResult clusterAnalysis,
            TimeAnimator animation,
            double threshold)
        {
            MonteCarloResult = mcResult;
            ClusterAnalysis = clusterAnalysis;
            Animation = animation;
            Threshold = threshold;
            Grid = mcResult.Grid;
            TimeFrame = mcResult.TimeFrame;
        }

        /// <summary>Deterministic (single-scenario) constructor.</summary>
        internal ScenarioResult(
            List<GridSnapshot> snapshots,
            GeoGrid grid,
            TimeFrame timeFrame)
        {
            Snapshots = snapshots;
            Grid = grid;
            TimeFrame = timeFrame;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Convenience queries
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Instantaneous exceedance probability at a position and time.
        /// Requires Monte Carlo data (returns 0 for deterministic runs).
        /// </summary>
        public double ProbabilityAt(Vector position, double timeSeconds)
        {
            return Animation?.ProbabilityAt(position, timeSeconds) ?? 0;
        }

        /// <summary>
        /// Cumulative exceedance probability: fraction of (filtered) scenarios
        /// where concentration exceeded the threshold at the given cell in
        /// any time step up to <paramref name="timeSeconds"/>.
        /// </summary>
        public double CumulativeProbabilityAt(Vector position, double timeSeconds)
        {
            return Animation?.CumulativeProbabilityAt(position, timeSeconds) ?? 0;
        }

        /// <summary>
        /// Returns the <see cref="ProbabilityMap"/> at the given time-step index.
        /// Null for deterministic runs.
        /// </summary>
        public ProbabilityMap ProbabilityMapAt(int timeIndex)
        {
            return Animation?[timeIndex];
        }

        // ═══════════════════════════════════════════════════════════════
        //  Export — GeoJSON
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Export the full time animation (or deterministic snapshots) to GeoJSON.
        /// </summary>
        public string ExportGeoJson(ExportMetadata metadata = null)
        {
            if (Animation != null)
                return GeoJsonExporter.ToGeoJson(Animation, metadata);

            if (Snapshots != null && Snapshots.Count > 0)
                return GeoJsonExporter.ToGeoJson(Snapshots[0], metadata);

            return "{}";
        }

        /// <summary>Export to a GeoJSON file.</summary>
        public void ExportGeoJson(string path, ExportMetadata metadata = null)
        {
            if (Animation != null)
            {
                GeoJsonExporter.Save(Animation, path, metadata);
                return;
            }
            if (Snapshots != null && Snapshots.Count > 0)
                GeoJsonExporter.Save(Snapshots[0], path, metadata);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Export — Unity binary
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Export to the compact Unity binary format (.bin).
        /// </summary>
        public void ExportUnity(string path)
        {
            if (Animation != null && ClusterAnalysis != null)
            {
                var meanSnaps = ClusterAnalysis
                    .GetClusterMeanSnapshots(ClusterAnalysis.DominantCluster)
                    .ToArray();
                UnityBinaryExporter.Save(meanSnaps, Animation, path);
            }
            else if (Animation != null)
            {
                UnityBinaryExporter.Save(Animation, path);
            }
            else if (Snapshots != null)
            {
                UnityBinaryExporter.Save(Snapshots.ToArray(), Grid, TimeFrame, path);
            }
        }

        // ═══════════════════════════════════════════════════════════════
        //  Export — Cesium (CZML)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Export the time animation to CZML for Cesium visualisation.
        /// </summary>
        public void ExportCesium(string path, string name = "Plume Risk", double minProbability = 0)
        {
            if (Animation != null)
                CesiumExporter.Save(Animation, path, name, minProbability);
        }
    }
}
