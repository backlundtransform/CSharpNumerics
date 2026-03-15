using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using CSharpNumerics.Statistics.MonteCarlo;
using CSharpNumerics.Statistics.Random;
using System;
using System.Linq;

namespace NumericsTests
{
    [TestClass]
    public class PlumeMonteCarloTests
    {
        private GeoGrid SmallGrid => new GeoGrid(-100, 100, -100, 100, 0, 0, 50);
        private TimeFrame ShortTime => new TimeFrame(0, 0, 1);

        // ════════════════════════════════════════════════════════════════
        //  RunBatch — basic
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void RunBatch_ProducesCorrectMatrixShape()
        {
            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime);

            var result = model.RunBatch(iterations: 10, seed: 42);

            Assert.AreEqual(10, result.Iterations);
            Assert.AreEqual(SmallGrid.CellCount * ShortTime.Count, result.FeatureCount);
            Assert.AreEqual(10, result.Snapshots.Count);
        }

        [TestMethod]
        public void RunBatch_Deterministic_SameResults()
        {
            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime);

            // Without variation, all rows should be identical (deterministic)
            var result = model.RunBatch(5, seed: 42);

            var row0 = result.GetScenarioVector(0);
            for (int i = 1; i < result.Iterations; i++)
            {
                var row = result.GetScenarioVector(i);
                for (int j = 0; j < row.Length; j++)
                    Assert.AreEqual(row0[j], row[j], 1e-20,
                        $"Row {i}, col {j} should match row 0 (no variation).");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  RunBatch — with variation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void RunBatch_WithVariation_ProducesDistinctScenarios()
        {
            var variation = new ScenarioVariation()
                .WindSpeed(5, 15)
                .WindDirectionJitter(20)
                .EmissionRate(3, 7);

            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime,
                variation: variation);

            var result = model.RunBatch(20, seed: 42);

            // At least some rows should differ
            var row0 = result.GetScenarioVector(0);
            bool anyDiff = false;
            for (int i = 1; i < result.Iterations && !anyDiff; i++)
            {
                var row = result.GetScenarioVector(i);
                for (int j = 0; j < row.Length; j++)
                    if (Math.Abs(row0[j] - row[j]) > 1e-20) { anyDiff = true; break; }
            }
            Assert.IsTrue(anyDiff, "With variation, scenarios should differ.");
        }

        [TestMethod]
        public void RunBatch_WithStabilityWeights_ProducesDistinctScenarios()
        {
            var variation = new ScenarioVariation()
                .SetStabilityWeights(a: 0.3, d: 0.4, f: 0.3);

            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime,
                variation: variation);

            var result = model.RunBatch(30, seed: 42);

            // Different stability classes produce different dispersion → different values
            var row0 = result.GetScenarioVector(0);
            bool anyDiff = false;
            for (int i = 1; i < result.Iterations && !anyDiff; i++)
            {
                var row = result.GetScenarioVector(i);
                for (int j = 0; j < row.Length; j++)
                    if (Math.Abs(row0[j] - row[j]) > 1e-20) { anyDiff = true; break; }
            }
            Assert.IsTrue(anyDiff, "With stability weights, scenarios should differ.");
        }

        // ════════════════════════════════════════════════════════════════
        //  RunBatch — with multiple time steps (transient)
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void RunBatch_Transient_MatrixIncludesAllTimeSteps()
        {
            var tf = new TimeFrame(10, 30, 10); // 3 time steps
            var grid = new GeoGrid(-50, 50, -50, 50, 0, 0, 50); // Nx=3, Ny=3 → 9 cells

            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 50,
                sourcePosition: new Vector(0, 0, 50),
                grid: grid,
                timeFrame: tf,
                mode: PlumeMode.Transient,
                releaseSeconds: 5);

            var result = model.RunBatch(5, seed: 42);

            // Columns = cells × timeSteps = 9 × 3 = 27
            Assert.AreEqual(9 * 3, result.FeatureCount);
            Assert.AreEqual(3, result.Snapshots[0].Count);
        }

        // ════════════════════════════════════════════════════════════════
        //  GetCellDistribution
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void GetCellDistribution_ReturnsCorrectValues()
        {
            var variation = new ScenarioVariation().WindSpeed(5, 15);

            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime,
                variation: variation);

            var result = model.RunBatch(10, seed: 42);

            // Pick an arbitrary cell
            int cellIdx = SmallGrid.NearestFlatIndex(new Vector(50, 0, 0));
            var dist = result.GetCellDistribution(cellIdx, 0);

            Assert.AreEqual(10, dist.Length);

            // Values should match the scenario matrix
            for (int i = 0; i < 10; i++)
                Assert.AreEqual(result.ScenarioMatrix.values[i, cellIdx], dist[i], 1e-20);
        }

        // ════════════════════════════════════════════════════════════════
        //  IMonteCarloModel (scalar) — compatible with MonteCarloSimulator
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void Evaluate_ReturnsPositiveScalar()
        {
            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime);

            var rng = new RandomGenerator(42);
            double value = model.Evaluate(rng);

            Assert.IsTrue(value >= 0, "Peak concentration should be non-negative.");
        }

        [TestMethod]
        public void Evaluate_WorksWithMonteCarloSimulator()
        {
            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime,
                variation: new ScenarioVariation().WindSpeed(5, 15));

            var sim = new MonteCarloSimulator(seed: 42);
            var result = sim.Run(model, iterations: 20);

            Assert.AreEqual(20, result.Count);
            Assert.IsTrue(result.Mean > 0, "Mean peak concentration should be positive.");
            Assert.IsTrue(result.StandardDeviation >= 0);
        }

        // ════════════════════════════════════════════════════════════════
        //  Seeded reproducibility
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        public void RunBatch_SameSeed_Reproducible()
        {
            var variation = new ScenarioVariation().WindSpeed(5, 15).WindDirectionJitter(10);

            PlumeMonteCarloModel MakeModel() => new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime,
                variation: variation);

            var r1 = MakeModel().RunBatch(10, seed: 123);
            var r2 = MakeModel().RunBatch(10, seed: 123);

            for (int i = 0; i < r1.Iterations; i++)
            {
                var v1 = r1.GetScenarioVector(i);
                var v2 = r2.GetScenarioVector(i);
                for (int j = 0; j < v1.Length; j++)
                    Assert.AreEqual(v1[j], v2[j], 1e-20,
                        $"Iteration {i}, feature {j}: same seed should produce identical results.");
            }
        }

        // ════════════════════════════════════════════════════════════════
        //  Validation
        // ════════════════════════════════════════════════════════════════

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void RunBatch_ZeroIterations_Throws()
        {
            var model = new PlumeMonteCarloModel(
                emissionRate: 5.0,
                windSpeed: 10,
                windDirection: new Vector(1, 0, 0),
                stackHeight: 0,
                sourcePosition: new Vector(0, 0, 0),
                grid: SmallGrid,
                timeFrame: ShortTime);

            model.RunBatch(0);
        }
    }
}
