using CSharpNumerics.Physics.Thermodynamics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest
{
    [TestClass]
    public class HeatSlabTests
    {
        // Mild steel.
        private const double K = 45;        // W/(m·K)
        private const double Rho = 7850;    // kg/m³
        private const double Cp = 490;      // J/(kg·K)
        private const double Ambient = 293.15;

        private static HeatSlab SteelPlate(double thickness = 0.010, int nodes = 5) =>
            new HeatSlab(thickness, K, Rho, Cp, Ambient, nodes);

        // ════════════════════════════════════════════
        //  Analytical checks
        // ════════════════════════════════════════════

        [TestMethod]
        public void HeatSlab_SteadyState_MatchesSeriesResistance()
        {
            // At steady state the flux through convection, conduction and
            // convection in series is q = ΔT / (1/h₁ + L/k + 1/h₂), and each
            // surface temperature follows from its own film drop.
            const double hotGas = 800, hHot = 25, hCold = 8;
            var slab = SteelPlate();

            double resistance = 1 / hHot + 0.010 / K + 1 / hCold;
            double q = (hotGas - Ambient) / resistance;
            double expectedHotSurface = hotGas - q / hHot;
            double expectedColdSurface = Ambient + q / hCold;

            for (int i = 0; i < 6000; i++)
                slab.Step(1.0, hotGas, hHot, Ambient, hCold);

            Assert.AreEqual(expectedHotSurface, slab.HotSideTemperature, 0.5);
            Assert.AreEqual(expectedColdSurface, slab.ColdSideTemperature, 0.5);
        }

        [TestMethod]
        public void HeatSlab_ThinPlate_FollowsTheLumpedCapacitanceSolution()
        {
            // Biot = hL/k ≈ 0.001, so the slab is effectively isothermal and
            // T(t) = T_gas − (T_gas − T₀)·exp(−t/τ) with τ = ρ·c·L / (2h)
            // when both faces see the same gas.
            const double gas = 400, h = 10, L = 0.005;
            var slab = new HeatSlab(L, K, Rho, Cp, Ambient);

            double tau = Rho * Cp * L / (2 * h);
            int steps = (int)Math.Round(tau);
            for (int i = 0; i < steps; i++)
                slab.Step(1.0, gas, h, gas, h);

            double expected = gas - (gas - Ambient) / Math.E;
            Assert.AreEqual(expected, slab.AverageTemperature, (gas - Ambient) * 0.01);
        }

        [TestMethod]
        public void HeatSlab_AbsorbedEnergy_MatchesTheTemperatureRise()
        {
            // With the cold side insulated, every joule reported as absorbed
            // must show up as stored heat: ∫q·dt = ρ·c·L·ΔT̄. Many nodes so the
            // unweighted node mean is close to the mass-weighted one.
            var slab = new HeatSlab(0.010, K, Rho, Cp, Ambient, nodeCount: 21);

            double absorbed = 0;
            for (int i = 0; i < 600; i++)
                absorbed += slab.Step(1.0, 700, 30, Ambient, 0);

            double stored = Rho * Cp * 0.010 * (slab.AverageTemperature - Ambient);
            Assert.AreEqual(absorbed, stored, absorbed * 0.02);
        }

        // ════════════════════════════════════════════
        //  Behaviour
        // ════════════════════════════════════════════

        [TestMethod]
        public void HeatSlab_ColdSideLagsTheHotSide()
        {
            var slab = SteelPlate();
            slab.Step(5.0, 900, 50, Ambient, 8);

            Assert.IsTrue(slab.HotSideTemperature > slab.ColdSideTemperature,
                "the exposed face must lead while heat is flowing inward");
            Assert.IsTrue(slab.ColdSideTemperature >= Ambient);
        }

        [TestMethod]
        public void HeatSlab_InsulatedColdSide_ApproachesTheGasTemperature()
        {
            var slab = SteelPlate();
            for (int i = 0; i < 20000; i++)
                slab.Step(1.0, 600, 25, Ambient, 0);

            Assert.AreEqual(600, slab.HotSideTemperature, 1.0);
            Assert.AreEqual(600, slab.ColdSideTemperature, 1.0);
        }

        [TestMethod]
        public void HeatSlab_NoExchange_HoldsItsTemperature()
        {
            var slab = SteelPlate();
            double absorbed = slab.Step(60, 1000, 0, 100, 0);

            Assert.AreEqual(0, absorbed, 1e-12);
            Assert.AreEqual(Ambient, slab.HotSideTemperature, 1e-9);
            Assert.AreEqual(Ambient, slab.ColdSideTemperature, 1e-9);
        }

        [TestMethod]
        public void HeatSlab_LongTimeStep_IsSubSteppedNotUnstable()
        {
            // A 60 s step is far beyond the explicit stability limit for a
            // 2.5 mm node spacing; the slab must sub-step rather than blow up.
            var slab = SteelPlate();
            slab.Step(60, 800, 25, Ambient, 8);

            Assert.IsFalse(double.IsNaN(slab.HotSideTemperature));
            Assert.IsTrue(slab.HotSideTemperature > Ambient);
            Assert.IsTrue(slab.HotSideTemperature < 800,
                "the plate cannot overshoot the gas heating it");
        }

        // ════════════════════════════════════════════
        //  The number that matters for a ship
        // ════════════════════════════════════════════

        [TestMethod]
        public void HeatSlab_BareSteelDeck_FailsTheInsulationCriterionInMinutes()
        {
            // The A-class fire test criterion limits the unexposed-side rise to
            // 140 K on average. A bare 10 mm deck plate over an 800 K fire
            // compartment blows through that in well under fifteen minutes —
            // which is exactly why A-60 divisions carry insulation, and the
            // mechanism by which a fire climbs a ship without any opening.
            var slab = SteelPlate();

            double t = 0;
            while (slab.ColdSideTemperature < Ambient + 140 && t < 3600)
            {
                slab.Step(1.0, 800, 25, Ambient, 8);
                t += 1.0;
            }

            Assert.IsTrue(t < 900,
                "bare steel should fail the 140 K rise within 15 minutes, took " + t + " s");
        }
    }
}
