using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;
using CSharpNumerics.Physics.Materials.Nuclear.DecayChains;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using DecayCalc = CSharpNumerics.Physics.Materials.Nuclear.Decay.Decay;

namespace NumericTest
{
    [TestClass]
    public class DecayTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  Decay static class
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Activity_FromMass()
        {
            // 1 kg of Cs-137 should give ~3.2e12 Bq
            double a = DecayCalc.Activity(1.0, Isotope.Cs137);
            Assert.AreEqual(Isotope.Cs137.SpecificActivity, a, 1e6);
        }

        [TestMethod]
        public void Activity_StableIsZero()
        {
            Assert.AreEqual(0.0, DecayCalc.Activity(1.0, Isotope.Xe131));
        }

        [TestMethod]
        public void RemainingMass_AfterOneHalfLife()
        {
            double m = DecayCalc.RemainingMass(1.0, Isotope.Cs137.HalfLife, Isotope.Cs137);
            Assert.AreEqual(0.5, m, 1e-6);
        }

        [TestMethod]
        public void RemainingMass_Stable_NoDecay()
        {
            double m = DecayCalc.RemainingMass(5.0, 1e10, Isotope.Xe131);
            Assert.AreEqual(5.0, m);
        }

        [TestMethod]
        public void ActivityAtTime_DecaysExponentially()
        {
            double a0 = 1e6;
            double t = Isotope.I131.HalfLife; // ~8 days
            double a = DecayCalc.ActivityAtTime(a0, t, Isotope.I131);
            Assert.AreEqual(a0 / 2, a, 1.0);
        }

        [TestMethod]
        public void ConcentrationToActivity_And_Back()
        {
            double c = 1e-6; // kg/m³
            double a = DecayCalc.ConcentrationToActivity(c, Isotope.Cs137);
            Assert.IsTrue(a > 0);

            double cBack = DecayCalc.ActivityToConcentration(a, Isotope.Cs137);
            Assert.AreEqual(c, cBack, 1e-18);
        }

        [TestMethod]
        public void HalfLivesElapsed()
        {
            double n = DecayCalc.HalfLivesElapsed(Isotope.Cs137.HalfLife * 3, Isotope.Cs137);
            Assert.AreEqual(3.0, n, 1e-10);
        }

        [TestMethod]
        public void DecayConstant()
        {
            double lambda = DecayCalc.DecayConstant(Isotope.Cs137);
            Assert.AreEqual(Isotope.Cs137.Lambda, lambda);
        }

        // ═══════════════════════════════════════════════════════════════
        //  DecayChain
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void DecayChain_Cs137_HasTwoSteps()
        {
            var chain = DecayChain.Cs137Chain();
            Assert.AreEqual(2, chain.StepCount);
            var isotopes = chain.AllIsotopes();
            Assert.AreEqual("Cs137", isotopes[0].Name);
            Assert.AreEqual("Ba137m", isotopes[1].Name);
            Assert.AreEqual("Ba137", isotopes[2].Name);
        }

        [TestMethod]
        public void DecayChain_I131_HasOneStep()
        {
            var chain = DecayChain.I131Chain();
            Assert.AreEqual(1, chain.StepCount);
        }

        [TestMethod]
        public void DecayChain_Evolve_MassConservation()
        {
            // Total mass should be roughly conserved.
            // With 94.7% branching to Ba-137m, ~5.3% goes to other branches
            // not modelled here, so expect ~97% mass recovery.
            var chain = DecayChain.Cs137Chain();
            var initial = new double[] { 1.0, 0, 0 }; // 1 kg Cs-137

            var result = chain.Evolve(initial, Isotope.Cs137.HalfLife);

            double totalMass = result[0] + result[1] + result[2];
            Assert.IsTrue(totalMass > 0.90 && totalMass <= 1.01,
                $"Total mass={totalMass:F4}, Cs137={result[0]:E3}, Ba137m={result[1]:E3}, Ba137={result[2]:E3}");
        }

        [TestMethod]
        public void DecayChain_Evolve_ParentDecaysToHalf()
        {
            var chain = DecayChain.Cs137Chain();
            var initial = new double[] { 1.0, 0, 0 };

            var result = chain.Evolve(initial, Isotope.Cs137.HalfLife);
            Assert.AreEqual(0.5, result[0], 0.01, $"Cs137 after 1 half-life: {result[0]}");
        }

        [TestMethod]
        public void DecayChain_Evolve_ZeroTime_ReturnsInitial()
        {
            var chain = DecayChain.Cs137Chain();
            var initial = new double[] { 1.0, 0, 0 };

            var result = chain.Evolve(initial, 0);
            Assert.AreEqual(1.0, result[0], 1e-10);
            Assert.AreEqual(0.0, result[1], 1e-10);
            Assert.AreEqual(0.0, result[2], 1e-10);
        }

        [TestMethod]
        public void DecayChain_I131_AfterManyHalfLives_AllXe131()
        {
            var chain = DecayChain.I131Chain();
            var initial = new double[] { 1.0, 0 };

            // After 100 half-lives, essentially all should be Xe-131
            var result = chain.Evolve(initial, Isotope.I131.HalfLife * 100);
            Assert.IsTrue(result[0] < 1e-20, $"I131 remaining: {result[0]:E3}");
            Assert.AreEqual(1.0, result[1], 0.01);
        }

        [TestMethod]
        public void DecayChain_EvolveActivity()
        {
            var chain = DecayChain.Cs137Chain();
            var initial = new double[] { 1e-6, 0, 0 }; // 1 µg Cs-137
            var activities = chain.EvolveActivity(initial, 0);
            // At t=0, only Cs-137 has activity
            Assert.IsTrue(activities[0] > 0);
            Assert.AreEqual(0.0, activities[1], 1e-10); // no Ba-137m yet (mass=0)
            Assert.AreEqual(0.0, activities[2]); // stable
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void DecayChain_Evolve_WrongArrayLength_Throws()
        {
            var chain = DecayChain.Cs137Chain();
            chain.Evolve(new double[] { 1.0 }, 100); // needs 3
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void DecayChain_AddStep_InvalidBranching_Throws()
        {
            new DecayChain(Isotope.Cs137).AddStep(Isotope.Ba137m, 1.5);
        }
    }
}
