using CSharpNumerics.Engines.GIS.Analysis;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Simulation;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class RiskScenarioTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Full fluent chain: MC + clustering
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void FullChain_MC_And_Clustering_Produces_ScenarioResult()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .WithVariation(v => v
                    .WindSpeed(8, 12)
                    .WindDirectionJitter(10))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(20, seed: 42)
                .AnalyzeWith(new KMeans { Seed = 42 }, new SilhouetteEvaluator(), 2, 4)
                .Build(threshold: 1e-8);

            Assert.IsNotNull(result);
            Assert.IsNotNull(result.MonteCarloResult);
            Assert.IsNotNull(result.ClusterAnalysis);
            Assert.IsNotNull(result.Animation);
            Assert.IsNotNull(result.Grid);
            Assert.IsNotNull(result.TimeFrame);
            Assert.AreEqual(1e-8, result.Threshold);
        }

        // ═══════════════════════════════════════════════════════════════
        //  MC without clustering
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void MC_Without_Clustering_SkipsAnalysis()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 1)
                .Build(threshold: 1e-6);

            Assert.IsNotNull(result.MonteCarloResult);
            Assert.IsNull(result.ClusterAnalysis);
            Assert.IsNotNull(result.Animation);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Deterministic single scenario
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void RunSingle_Produces_Deterministic_Snapshots()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithStability(StabilityClass.D)
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunSingle();

            Assert.IsNotNull(result.Snapshots);
            Assert.AreEqual(2, result.Snapshots.Count, "2 time steps");
            Assert.IsNull(result.MonteCarloResult);
            Assert.IsNull(result.ClusterAnalysis);
            Assert.IsNull(result.Animation);
        }

        // ═══════════════════════════════════════════════════════════════
        //  ProbabilityAt / CumulativeProbabilityAt
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ProbabilityAt_Returns_Values_In_Range()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(5, 15))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(20, seed: 42)
                .Build(threshold: 1e-8);

            double p = result.ProbabilityAt(new Vector(50, 0, 0), timeSeconds: 0);
            Assert.IsTrue(p >= 0.0 && p <= 1.0, $"Probability {p} out of range");
        }

        [TestMethod]
        public void CumulativeProbability_NonDecreasing()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(5, 15))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 120, 60)
                .RunMonteCarlo(20, seed: 42)
                .Build(threshold: 1e-8);

            var pos = new Vector(50, 0, 0);
            double cum0 = result.CumulativeProbabilityAt(pos, 0);
            double cum60 = result.CumulativeProbabilityAt(pos, 60);
            double cum120 = result.CumulativeProbabilityAt(pos, 120);

            Assert.IsTrue(cum60 >= cum0, "Cumulative probability must not decrease");
            Assert.IsTrue(cum120 >= cum60, "Cumulative probability must not decrease");
        }

        // ═══════════════════════════════════════════════════════════════
        //  Deterministic ProbabilityAt returns 0
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Deterministic_ProbabilityAt_Returns_Zero()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunSingle();

            Assert.AreEqual(0.0, result.ProbabilityAt(new Vector(50, 0, 0), 30));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Transient mode
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Transient_Mode_RunSingle()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithMode(PlumeMode.Transient, releaseSeconds: 10)
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunSingle();

            Assert.AreEqual(2, result.Snapshots.Count);
        }

        [TestMethod]
        public void Transient_Mode_MC()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithMode(PlumeMode.Transient, releaseSeconds: 10)
                .WithVariation(v => v.WindSpeed(8, 12))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build(threshold: 1e-8);

            Assert.IsNotNull(result.MonteCarloResult);
            Assert.AreEqual(10, result.MonteCarloResult.Iterations);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Validation
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void MissingGrid_Throws()
        {
            RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10);
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void MissingTimeFrame_Throws()
        {
            RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .RunMonteCarlo(10);
        }

        // ═══════════════════════════════════════════════════════════════
        //  ProbabilityMapAt
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ProbabilityMapAt_Returns_Correct_TimeStep()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build(threshold: 1e-8);

            var map0 = result.ProbabilityMapAt(0);
            var map1 = result.ProbabilityMapAt(1);

            Assert.IsNotNull(map0);
            Assert.IsNotNull(map1);
            Assert.AreEqual(0, map0.TimeIndex);
            Assert.AreEqual(1, map1.TimeIndex);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Reproducibility
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Seeded_MC_Is_Reproducible()
        {
            ScenarioResult Run() => RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(5, 15))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(10, seed: 123)
                .Build(threshold: 1e-8);

            var r1 = Run();
            var r2 = Run();

            var pos = new Vector(50, 0, 0);
            Assert.AreEqual(r1.ProbabilityAt(pos, 0), r2.ProbabilityAt(pos, 0), 1e-10);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Grid & TimeFrame propagation
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Result_Carries_Grid_And_TimeFrame()
        {
            var grid = new GeoGrid(-100, 100, -100, 100, 0, 0, 50);
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .OverGrid(grid)
                .OverTime(0, 120, 60)
                .RunMonteCarlo(10, seed: 42)
                .Build();

            Assert.AreSame(grid, result.Grid);
            Assert.AreEqual(3, result.TimeFrame.Count);
        }

        // ═══════════════════════════════════════════════════════════════
        //  ClusterAnalysis accessible from result
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void ClusterAnalysis_Accessible_After_AnalyzeWith()
        {
            var result = RiskScenario
                .ForGaussianPlume(5.0)
                .FromSource(new Vector(0, 0, 50))
                .WithWind(10, new Vector(1, 0, 0))
                .WithVariation(v => v.WindSpeed(5, 15).WindDirectionJitter(15))
                .OverGrid(new GeoGrid(-100, 100, -100, 100, 0, 0, 50))
                .OverTime(0, 60, 60)
                .RunMonteCarlo(20, seed: 42)
                .AnalyzeWith(new KMeans { Seed = 42 }, new SilhouetteEvaluator(), 2, 4)
                .Build(threshold: 1e-8);

            Assert.IsNotNull(result.ClusterAnalysis);
            Assert.IsTrue(result.ClusterAnalysis.BestClusterCount >= 2);
            Assert.IsTrue(result.ClusterAnalysis.Labels.Length == 20);
        }
    }
}
