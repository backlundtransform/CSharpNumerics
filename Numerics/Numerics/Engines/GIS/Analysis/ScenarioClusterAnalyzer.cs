using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Clustering.Results;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.GIS.Analysis
{
    /// <summary>
    /// Fluent entry point for clustering Monte Carlo scenario results.
    /// Wraps <see cref="ClusteringExperiment"/> to operate directly on
    /// <see cref="MonteCarloScenarioResult.ScenarioMatrix"/>.
    /// <para>
    /// Usage:
    /// <code>
    /// var analysis = ScenarioClusterAnalyzer
    ///     .For(mcResult)
    ///     .WithAlgorithm(new KMeans())
    ///     .TryClusterCounts(2, 6)
    ///     .WithEvaluator(new SilhouetteEvaluator())
    ///     .Run();
    /// </code>
    /// </para>
    /// </summary>
    public static class ScenarioClusterAnalyzer
    {
        /// <summary>Start building a cluster analysis for Monte Carlo results.</summary>
        public static ScenarioClusterBuilder For(MonteCarloScenarioResult scenarios)
            => new ScenarioClusterBuilder(scenarios);
    }

    /// <summary>
    /// Builder that configures and runs a clustering experiment on the
    /// scenario matrix from a Monte Carlo plume simulation.
    /// </summary>
    public class ScenarioClusterBuilder
    {
        private readonly MonteCarloScenarioResult _scenarios;
        private readonly List<IClusteringModel> _algorithms = new();
        private readonly List<IClusteringEvaluator> _evaluators = new();
        private (int min, int max)? _clusterRange;
        private ClusteringGrid _grid;

        internal ScenarioClusterBuilder(MonteCarloScenarioResult scenarios)
        {
            _scenarios = scenarios ?? throw new ArgumentNullException(nameof(scenarios));
        }

        /// <summary>Add a clustering algorithm.</summary>
        public ScenarioClusterBuilder WithAlgorithm(IClusteringModel model)
        {
            _algorithms.Add(model);
            return this;
        }

        /// <summary>Try cluster counts from min to max (inclusive).</summary>
        public ScenarioClusterBuilder TryClusterCounts(int min, int max)
        {
            _clusterRange = (min, max);
            return this;
        }

        /// <summary>Add a cluster-quality evaluator.</summary>
        public ScenarioClusterBuilder WithEvaluator(IClusteringEvaluator evaluator)
        {
            _evaluators.Add(evaluator);
            return this;
        }

        /// <summary>Use a full <see cref="ClusteringGrid"/> for hyper-parameter search.</summary>
        public ScenarioClusterBuilder WithGrid(ClusteringGrid grid)
        {
            _grid = grid;
            return this;
        }

        /// <summary>Execute the clustering experiment and return a domain-specific result.</summary>
        public ClusterAnalysisResult Run()
        {
            var builder = ClusteringExperiment.For(_scenarios.ScenarioMatrix);

            if (_grid != null)
            {
                builder = builder.WithGrid(_grid);
            }
            else
            {
                foreach (var alg in _algorithms)
                    builder = builder.WithAlgorithm(alg);
                if (_clusterRange.HasValue)
                    builder = builder.TryClusterCounts(_clusterRange.Value.min, _clusterRange.Value.max);
            }

            foreach (var eval in _evaluators)
                builder = builder.WithEvaluator(eval);

            var experimentResult = builder.Run();

            return new ClusterAnalysisResult(_scenarios, experimentResult);
        }
    }

    /// <summary>
    /// Result of clustering Monte Carlo plume scenarios.
    /// Provides cluster labels, the dominant cluster, and iteration look-ups.
    /// </summary>
    public class ClusterAnalysisResult
    {
        /// <summary>The Monte Carlo scenario data that was clustered.</summary>
        public MonteCarloScenarioResult Scenarios { get; }

        /// <summary>The full clustering experiment result (rankings, scores, elbow curve).</summary>
        public ClusteringExperimentResult ExperimentResult { get; }

        /// <summary>Cluster label assigned to each scenario iteration.</summary>
        public VectorN Labels => ExperimentResult.BestLabels;

        /// <summary>Number of clusters in the best solution.</summary>
        public int BestClusterCount => ExperimentResult.BestClusterCount;

        /// <summary>Primary evaluator score of the best solution.</summary>
        public double BestScore => ExperimentResult.BestScore;

        /// <summary>
        /// The cluster with the most assigned iterations.
        /// Represents the most probable outcome group.
        /// </summary>
        public int DominantCluster { get; }

        internal ClusterAnalysisResult(
            MonteCarloScenarioResult scenarios,
            ClusteringExperimentResult experimentResult)
        {
            Scenarios = scenarios;
            ExperimentResult = experimentResult;
            DominantCluster = FindDominantCluster();
        }

        /// <summary>
        /// Returns the iteration indices assigned to the given cluster.
        /// </summary>
        public int[] GetClusterIterations(int cluster)
        {
            var list = new List<int>();
            for (int i = 0; i < Labels.Length; i++)
                if ((int)Labels[i] == cluster)
                    list.Add(i);
            return list.ToArray();
        }

        /// <summary>
        /// Returns the mean <see cref="GridSnapshot"/> for each time step
        /// across all iterations in the given cluster.
        /// </summary>
        public List<GridSnapshot> GetClusterMeanSnapshots(int cluster)
        {
            var iterations = GetClusterIterations(cluster);
            if (iterations.Length == 0)
                return new List<GridSnapshot>();

            var grid = Scenarios.Grid;
            int cellCount = grid.CellCount;
            int timeCount = Scenarios.TimeFrame.Count;
            var snapshots = new List<GridSnapshot>(timeCount);

            for (int t = 0; t < timeCount; t++)
            {
                var mean = new double[cellCount];
                int colOffset = t * cellCount;

                foreach (int iter in iterations)
                {
                    for (int c = 0; c < cellCount; c++)
                        mean[c] += Scenarios.ScenarioMatrix.values[iter, colOffset + c];
                }

                for (int c = 0; c < cellCount; c++)
                    mean[c] /= iterations.Length;

                snapshots.Add(new GridSnapshot(grid, mean, Scenarios.TimeFrame.TimeAt(t), t));
            }

            return snapshots;
        }

        private int FindDominantCluster()
        {
            var counts = new Dictionary<int, int>();
            for (int i = 0; i < Labels.Length; i++)
            {
                int c = (int)Labels[i];
                if (!counts.ContainsKey(c)) counts[c] = 0;
                counts[c]++;
            }

            int dominant = 0;
            int maxCount = 0;
            foreach (var kv in counts)
            {
                if (kv.Value > maxCount)
                {
                    maxCount = kv.Value;
                    dominant = kv.Key;
                }
            }
            return dominant;
        }
    }
}
