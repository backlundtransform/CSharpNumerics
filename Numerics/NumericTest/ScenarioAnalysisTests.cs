using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.Numerics.Objects;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class ScenarioAnalysisTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Helper: build a synthetic MonteCarloScenarioResult
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Creates a simple 3×3×1 grid (9 cells), 2 time steps,
        /// with N iterations of hand-crafted concentration data.
        /// Group A (first half): high concentrations at cell 0.
        /// Group B (second half): high concentrations at cell 8.
        /// </summary>
        private static MonteCarloScenarioResult CreateSynthetic(int iterations = 20)
        {
            var grid = new GeoGrid(0, 20, 0, 20, 0, 0, 10); // 3×3×1 = 9 cells
            var tf = new TimeFrame(0, 60, 60);               // 2 time steps: t=0, t=60

            int cells = grid.CellCount;  // 9
            int times = tf.Count;        // 2
            int cols = cells * times;    // 18

            var data = new double[iterations, cols];
            var snapshots = new List<List<GridSnapshot>>(iterations);

            int halfA = iterations / 2;

            for (int i = 0; i < iterations; i++)
            {
                // Group A: high at cell 0, low elsewhere
                // Group B: high at cell 8, low elsewhere
                bool groupA = i < halfA;

                for (int t = 0; t < times; t++)
                {
                    int offset = t * cells;
                    for (int c = 0; c < cells; c++)
                    {
                        double val;
                        if (groupA)
                            val = c == 0 ? 10.0 + t * 5.0 : 0.1;
                        else
                            val = c == 8 ? 8.0 + t * 4.0 : 0.2;

                        data[i, offset + c] = val;
                    }
                }

                // Build matching GridSnapshot list
                var iterSnaps = new List<GridSnapshot>();
                for (int t = 0; t < times; t++)
                {
                    var vals = new double[cells];
                    for (int c = 0; c < cells; c++)
                        vals[c] = data[i, t * cells + c];
                    iterSnaps.Add(new GridSnapshot(grid, vals, tf.TimeAt(t), t));
                }
                snapshots.Add(iterSnaps);
            }

            return new MonteCarloScenarioResult(
                new Matrix(data), snapshots, grid, tf);
        }

        // ═══════════════════════════════════════════════════════════════
        //  ScenarioClusterAnalyzer tests
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ClusterAnalyzer_Finds_Two_Groups()
        {
            var mc = CreateSynthetic(20);

            var analysis = ScenarioClusterAnalyzer
                .For(mc)
                .WithAlgorithm(new KMeans { K = 2, Seed = 42 })
                .WithEvaluator(new SilhouetteEvaluator())
                .Run();

            Assert.AreEqual(2, analysis.BestClusterCount);
            Assert.AreEqual(20, analysis.Labels.Length);
        }

        [TestMethod]
        public void ClusterAnalyzer_Labels_Match_Iterations()
        {
            var mc = CreateSynthetic(30);

            var analysis = ScenarioClusterAnalyzer
                .For(mc)
                .WithAlgorithm(new KMeans { K = 2, Seed = 42 })
                .WithEvaluator(new SilhouetteEvaluator())
                .Run();

            // Every iteration should be in exactly one cluster
            var allIters = new HashSet<int>();
            for (int c = 0; c < analysis.BestClusterCount; c++)
                foreach (int idx in analysis.GetClusterIterations(c))
                    allIters.Add(idx);

            Assert.AreEqual(30, allIters.Count, "All iterations accounted for");
        }

        [TestMethod]
        public void ClusterAnalyzer_DominantCluster_Is_Largest()
        {
            // 30 iterations: 15 in each group → both equal, dominant is one of them
            var mc = CreateSynthetic(30);

            var analysis = ScenarioClusterAnalyzer
                .For(mc)
                .WithAlgorithm(new KMeans { K = 2, Seed = 42 })
                .WithEvaluator(new SilhouetteEvaluator())
                .Run();

            int dom = analysis.DominantCluster;
            int domCount = analysis.GetClusterIterations(dom).Length;

            // Dominant should have the most iterations (or tied)
            for (int c = 0; c < analysis.BestClusterCount; c++)
            {
                Assert.IsTrue(domCount >= analysis.GetClusterIterations(c).Length,
                    $"Cluster {dom} should have >= iterations than cluster {c}");
            }
        }

        [TestMethod]
        public void ClusterAnalyzer_TryClusterCounts_Selects_Best_K()
        {
            var mc = CreateSynthetic(20);

            var analysis = ScenarioClusterAnalyzer
                .For(mc)
                .WithAlgorithm(new KMeans { Seed = 42 })
                .TryClusterCounts(2, 4)
                .WithEvaluator(new SilhouetteEvaluator())
                .Run();

            // With 2 clear groups, K=2 should be optimal
            Assert.AreEqual(2, analysis.BestClusterCount);
        }

        [TestMethod]
        public void ClusterAnalyzer_ExperimentResult_Is_Accessible()
        {
            var mc = CreateSynthetic(20);

            var analysis = ScenarioClusterAnalyzer
                .For(mc)
                .WithAlgorithm(new KMeans { K = 2, Seed = 42 })
                .WithEvaluator(new SilhouetteEvaluator())
                .Run();

            Assert.IsNotNull(analysis.ExperimentResult);
            Assert.IsTrue(analysis.ExperimentResult.Rankings.Count > 0);
            Assert.IsTrue(analysis.BestScore > 0, "Silhouette score should be positive for 2 clear clusters");
        }

        [TestMethod]
        public void ClusterAnalyzer_GetClusterMeanSnapshots()
        {
            var mc = CreateSynthetic(20);

            var analysis = ScenarioClusterAnalyzer
                .For(mc)
                .WithAlgorithm(new KMeans { K = 2, Seed = 42 })
                .WithEvaluator(new SilhouetteEvaluator())
                .Run();

            var meanSnaps = analysis.GetClusterMeanSnapshots(analysis.DominantCluster);

            Assert.AreEqual(2, meanSnaps.Count, "One snapshot per time step");
            Assert.AreEqual(9, meanSnaps[0].Count, "9 cells");
            Assert.IsTrue(meanSnaps[0].Max() > 1.0, "Mean of high-concentration group should be significant");
        }

        // ═══════════════════════════════════════════════════════════════
        //  ProbabilityMap tests
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ProbabilityMap_AllScenarios_ValuesInZeroOne()
        {
            var mc = CreateSynthetic(20);

            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 1.0);

            for (int i = 0; i < pmap.CellCount; i++)
            {
                Assert.IsTrue(pmap[i] >= 0.0 && pmap[i] <= 1.0,
                    $"Cell {i} probability {pmap[i]} out of range");
            }
        }

        [TestMethod]
        public void ProbabilityMap_HighThreshold_Fewer_Exceedances()
        {
            var mc = CreateSynthetic(20);

            var lowT = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 0.05);
            var highT = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 100.0);

            // With threshold 0.05, nearly every cell exceeds in every scenario
            Assert.IsTrue(lowT.Max() > 0.5, "Low threshold → high probability");
            // With threshold 100, nothing exceeds
            Assert.AreEqual(0.0, highT.Max(), "Very high threshold → zero probability");
        }

        [TestMethod]
        public void ProbabilityMap_Cell0_High_For_GroupA()
        {
            var mc = CreateSynthetic(20);

            // No filter: cell 0 has value 10 for first 10 iterations → P = 10/20 = 0.5
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            Assert.AreEqual(0.5, pmap[0], 1e-10, "Half of iterations exceed threshold at cell 0");
        }

        [TestMethod]
        public void ProbabilityMap_WithIterationFilter()
        {
            var mc = CreateSynthetic(20);

            // Filter to only the first 10 iterations (group A)
            int[] groupA = Enumerable.Range(0, 10).ToArray();
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0, groupA);

            // Cell 0: all 10 iterations have value 10 → P = 1.0
            Assert.AreEqual(1.0, pmap[0], 1e-10, "All group-A iterations exceed at cell 0");
            Assert.AreEqual(10, pmap.ScenarioCount);
        }

        [TestMethod]
        public void ProbabilityMap_At_Position_Nearest()
        {
            var mc = CreateSynthetic(20);
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            // Position (0,0,0) is cell 0 centre. Value should match pmap[0].
            double p = pmap.At(new Vector(0, 0, 0));
            Assert.AreEqual(pmap[0], p, 1e-10);
        }

        [TestMethod]
        public void ProbabilityMap_CellsAbove_Returns_High_Prob_Cells()
        {
            var mc = CreateSynthetic(20);
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            // Only cell 0 (groupA) and cell 8 (groupB) exceed 5.0
            var danger = pmap.CellsAbove(0.3);
            Assert.IsTrue(danger.Count >= 1 && danger.Count <= 2,
                $"Expected 1-2 cells above P=0.3, got {danger.Count}");
        }

        [TestMethod]
        public void ProbabilityMap_ThreeDimensional_Indexer()
        {
            var mc = CreateSynthetic(20);
            var pmap = ProbabilityMap.Build(mc, timeIndex: 0, threshold: 5.0);

            // Cell (0,0,0) should equal flat index 0
            Assert.AreEqual(pmap[0], pmap[0, 0, 0], 1e-10);
        }

        // ═══════════════════════════════════════════════════════════════
        //  TimeAnimator tests
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void TimeAnimator_Build_ProducesCorrectCount()
        {
            var mc = CreateSynthetic(20);
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            Assert.AreEqual(2, anim.Count, "2 time steps → 2 maps");
            Assert.IsNotNull(anim[0]);
            Assert.IsNotNull(anim[1]);
        }

        [TestMethod]
        public void TimeAnimator_ProbabilityAt_Matches_Map()
        {
            var mc = CreateSynthetic(20);
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            var pos = new Vector(0, 0, 0);
            double fromAnim = anim.ProbabilityAt(pos, timeSeconds: 0);
            double fromMap = anim[0].At(pos);

            Assert.AreEqual(fromMap, fromAnim, 1e-10);
        }

        [TestMethod]
        public void TimeAnimator_CumulativeProbability_GrowsOrStays()
        {
            var mc = CreateSynthetic(20);
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            // At cell 0: group A always exceeds → cumulative at t=0 is 0.5, at t=60 is 0.5
            var pos = new Vector(0, 0, 0);
            double cumT0 = anim.CumulativeProbabilityAt(pos, timeSeconds: 0);
            double cumT1 = anim.CumulativeProbabilityAt(pos, timeSeconds: 60);

            Assert.IsTrue(cumT1 >= cumT0,
                "Cumulative probability must be non-decreasing over time");
        }

        [TestMethod]
        public void TimeAnimator_WithClusterFilter()
        {
            var mc = CreateSynthetic(20);

            // Filter to group A only (iterations 0..9)
            int[] groupA = Enumerable.Range(0, 10).ToArray();
            var anim = TimeAnimator.Build(mc, threshold: 5.0, groupA);

            // Cell 0 should have P=1.0 at t=0 (all group-A iterations have value 10 at cell 0)
            double p = anim.ProbabilityAt(new Vector(0, 0, 0), timeSeconds: 0);
            Assert.AreEqual(1.0, p, 1e-10);
        }

        [TestMethod]
        public void TimeAnimator_CumulativeProbability_Cell8_GroupB()
        {
            var mc = CreateSynthetic(20);

            // Cell 8 has high values only in group B (iterations 10..19)
            // Across all 20: 10 exceed → cumulative should be 0.5
            var pos = new Vector(20, 20, 0);
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            double cum = anim.CumulativeProbabilityAt(pos, timeSeconds: 60);
            Assert.AreEqual(0.5, cum, 1e-10,
                "Group B (10/20) exceeds threshold at cell 8");
        }

        [TestMethod]
        public void TimeAnimator_ToArray_ReturnsClone()
        {
            var mc = CreateSynthetic(20);
            var anim = TimeAnimator.Build(mc, threshold: 5.0);

            var arr = anim.ToArray();
            Assert.AreEqual(anim.Count, arr.Length);
            Assert.AreSame(anim[0], arr[0], "Same ProbabilityMap object");
        }
    }
}
