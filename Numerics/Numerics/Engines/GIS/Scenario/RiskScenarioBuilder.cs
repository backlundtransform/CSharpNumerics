using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using System;

namespace CSharpNumerics.Engines.GIS.Scenario
{
    /// <summary>
    /// Fluent builder that collects all configuration for a risk scenario and
    /// drives execution through the pipeline stages:
    /// configure → Monte Carlo → clustering analysis → result.
    /// </summary>
    public class RiskScenarioBuilder
    {
        // ── Baseline physics ─────────────────────────────────────────
        private readonly double _emissionRate;
        private Vector _sourcePosition;
        private double _windSpeed = 5;
        private Vector _windDirection = new Vector(1, 0, 0);
        private double _stackHeight;
        private StabilityClass _stability = StabilityClass.D;
        private PlumeMode _mode = PlumeMode.SteadyState;
        private double _releaseSeconds;

        // ── Variation ────────────────────────────────────────────────
        private ScenarioVariation _variation;

        // ── Grid & time ──────────────────────────────────────────────
        private GeoGrid _grid;
        private TimeFrame _timeFrame;

        internal RiskScenarioBuilder(double emissionRate)
        {
            _emissionRate = emissionRate;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Configuration steps (return self)
        // ═══════════════════════════════════════════════════════════════

        /// <summary>Set the source position. Stack height is derived from z.</summary>
        public RiskScenarioBuilder FromSource(Vector position)
        {
            _sourcePosition = position;
            _stackHeight = position.z;
            return this;
        }

        /// <summary>Set baseline wind speed and direction.</summary>
        public RiskScenarioBuilder WithWind(double speed, Vector direction)
        {
            _windSpeed = speed;
            _windDirection = direction;
            return this;
        }

        /// <summary>Set the Pasquill–Gifford atmospheric stability class.</summary>
        public RiskScenarioBuilder WithStability(StabilityClass stability)
        {
            _stability = stability;
            return this;
        }

        /// <summary>Set physics mode (SteadyState or Transient).</summary>
        public RiskScenarioBuilder WithMode(PlumeMode mode, double releaseSeconds = 0)
        {
            _mode = mode;
            _releaseSeconds = releaseSeconds;
            return this;
        }

        /// <summary>Configure Monte Carlo parameter variation ranges.</summary>
        public RiskScenarioBuilder WithVariation(Action<ScenarioVariation> configure)
        {
            _variation = new ScenarioVariation();
            configure(_variation);
            return this;
        }

        /// <summary>Define the spatial grid.</summary>
        public RiskScenarioBuilder OverGrid(GeoGrid grid)
        {
            _grid = grid ?? throw new ArgumentNullException(nameof(grid));
            return this;
        }

        /// <summary>Define the time range.</summary>
        public RiskScenarioBuilder OverTime(double startSeconds, double endSeconds, double stepSeconds)
        {
            _timeFrame = new TimeFrame(startSeconds, endSeconds, stepSeconds);
            return this;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Terminal: Run Monte Carlo
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Execute N Monte Carlo iterations and transition to the
        /// Monte Carlo stage where clustering can be applied.
        /// </summary>
        /// <param name="iterations">Number of stochastic trials.</param>
        /// <param name="seed">Optional random seed for reproducibility.</param>
        public MonteCarloStageResult RunMonteCarlo(int iterations, int? seed = null)
        {
            Validate();

            var model = new PlumeMonteCarloModel(
                _emissionRate, _windSpeed, _windDirection, _stackHeight,
                _sourcePosition, _grid, _timeFrame, _variation,
                _stability, _mode, _releaseSeconds);

            var mcResult = model.RunBatch(iterations, seed);
            return new MonteCarloStageResult(mcResult);
        }

        /// <summary>
        /// Run a single deterministic scenario (no Monte Carlo) and return
        /// a <see cref="ScenarioResult"/> directly.
        /// </summary>
        public ScenarioResult RunSingle()
        {
            Validate();

            var sim = new PlumeSimulator(
                _emissionRate, _windSpeed, _windDirection,
                _stackHeight, _sourcePosition, _stability, _mode);
            sim.ReleaseSeconds = _releaseSeconds;

            var snapshots = sim.Run(_grid, _timeFrame);
            return new ScenarioResult(snapshots, _grid, _timeFrame);
        }

        private void Validate()
        {
            if (_grid == null) throw new InvalidOperationException("Grid not set. Call OverGrid().");
            if (_timeFrame == null) throw new InvalidOperationException("Time frame not set. Call OverTime().");
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Stage 2: Monte Carlo completed — can apply clustering
    // ═══════════════════════════════════════════════════════════════════

    /// <summary>
    /// Intermediate result after Monte Carlo execution.
    /// From here you can either analyse with clustering or build a result directly.
    /// </summary>
    public class MonteCarloStageResult
    {
        internal MonteCarloScenarioResult McResult { get; }

        internal MonteCarloStageResult(MonteCarloScenarioResult mcResult)
        {
            McResult = mcResult;
        }

        /// <summary>
        /// Run clustering analysis on the scenario matrix.
        /// </summary>
        /// <param name="algorithm">Clustering algorithm (e.g. <c>new KMeans()</c>).</param>
        /// <param name="evaluator">Quality evaluator (e.g. <c>new SilhouetteEvaluator()</c>).</param>
        /// <param name="minK">Minimum cluster count to try.</param>
        /// <param name="maxK">Maximum cluster count to try.</param>
        public AnalysisStageResult AnalyzeWith(
            IClusteringModel algorithm,
            IClusteringEvaluator evaluator,
            int minK = 2,
            int maxK = 6)
        {
            var analysis = ScenarioClusterAnalyzer
                .For(McResult)
                .WithAlgorithm(algorithm)
                .TryClusterCounts(minK, maxK)
                .WithEvaluator(evaluator)
                .Run();

            return new AnalysisStageResult(McResult, analysis);
        }

        /// <summary>
        /// Run clustering analysis using a full <see cref="ClusteringGrid"/>
        /// for hyper-parameter search.
        /// </summary>
        public AnalysisStageResult AnalyzeWith(
            ClusteringGrid grid,
            IClusteringEvaluator evaluator)
        {
            var analysis = ScenarioClusterAnalyzer
                .For(McResult)
                .WithGrid(grid)
                .WithEvaluator(evaluator)
                .Run();

            return new AnalysisStageResult(McResult, analysis);
        }

        /// <summary>
        /// Skip clustering and build a <see cref="ScenarioResult"/> directly
        /// from the Monte Carlo output.
        /// </summary>
        /// <param name="threshold">Concentration threshold for probability maps (kg/m³).</param>
        public ScenarioResult Build(double threshold = 1e-6)
        {
            var animator = TimeAnimator.Build(McResult, threshold);
            return new ScenarioResult(McResult, null, animator, threshold);
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    //  Stage 3: Clustering completed — build final result
    // ═══════════════════════════════════════════════════════════════════

    /// <summary>
    /// Intermediate result after clustering analysis.
    /// Call <see cref="Build"/> to produce the final <see cref="ScenarioResult"/>.
    /// </summary>
    public class AnalysisStageResult
    {
        internal MonteCarloScenarioResult McResult { get; }
        internal ClusterAnalysisResult Analysis { get; }

        internal AnalysisStageResult(
            MonteCarloScenarioResult mcResult,
            ClusterAnalysisResult analysis)
        {
            McResult = mcResult;
            Analysis = analysis;
        }

        /// <summary>
        /// Build the final <see cref="ScenarioResult"/> with probability maps
        /// computed over the dominant cluster.
        /// </summary>
        /// <param name="threshold">Concentration threshold for probability maps (kg/m³).</param>
        public ScenarioResult Build(double threshold = 1e-6)
        {
            int[] dominantIters = Analysis.GetClusterIterations(Analysis.DominantCluster);
            var animator = TimeAnimator.Build(McResult, threshold, dominantIters);
            return new ScenarioResult(McResult, Analysis, animator, threshold);
        }
    }
}
