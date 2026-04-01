using CSharpNumerics.Physics.Materials.Fire;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using CSharpNumerics.Physics.Environmental.Fire;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest
{
    [TestClass]
    public class WildfireTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  FuelModel & FuelLibrary
        // ═══════════════════════════════════════════════════════════════

        #region FuelModel & FuelLibrary

        [TestMethod]
        public void FuelLibrary_ReturnsAll13AndersonModels()
        {
            var all = FuelLibrary.All();
            Assert.AreEqual(13, all.Count, "Should contain all 13 Anderson fuel models.");
        }

        [TestMethod]
        public void FuelLibrary_Get_ShortGrass_HasCorrectProperties()
        {
            var fuel = FuelLibrary.Get(FuelModelType.ShortGrass);

            Assert.AreEqual(FuelModelType.ShortGrass, fuel.Type);
            Assert.AreEqual("Short Grass", fuel.Name);
            Assert.AreEqual(11483, fuel.SurfaceAreaToVolumeRatio, 1.0,
                "σ = 3500 1/ft ≈ 11483 1/m");
            Assert.AreEqual(0.305, fuel.FuelBedDepth, 0.001,
                "δ = 1.0 ft ≈ 0.305 m");
            Assert.AreEqual(0.12, fuel.MoistureOfExtinction, 0.001);
            Assert.AreEqual(18608, fuel.LowHeatContent, 1.0);
            Assert.AreEqual(513.0, fuel.ParticleDensity, 0.1);
        }

        [TestMethod]
        public void FuelLibrary_Get_AllModelsRoundTrip()
        {
            foreach (FuelModelType type in Enum.GetValues(typeof(FuelModelType)))
            {
                if (type == FuelModelType.NoFuel) continue;

                Assert.IsTrue(FuelLibrary.TryGet(type, out var fuel),
                    $"FuelLibrary should contain model {type}.");
                Assert.AreEqual(type, fuel.Type);
                Assert.IsTrue(fuel.SurfaceAreaToVolumeRatio > 0, $"{type} σ > 0");
                Assert.IsTrue(fuel.FuelBedDepth > 0, $"{type} δ > 0");
                Assert.IsTrue(fuel.OvendryFuelLoad > 0, $"{type} w₀ > 0");
            }
        }

        [TestMethod]
        public void FuelLibrary_Get_NoFuel_ThrowsKeyNotFound()
        {
            Assert.ThrowsException<System.Collections.Generic.KeyNotFoundException>(
                () => FuelLibrary.Get(FuelModelType.NoFuel));
        }

        [TestMethod]
        public void FuelModel_Equality_ByType()
        {
            var a = FuelModel.ShortGrass;
            var b = FuelLibrary.Get(FuelModelType.ShortGrass);
            Assert.AreEqual(a, b);
            Assert.IsTrue(a == b);
        }

        [TestMethod]
        public void FuelModel_ToString_ContainsName()
        {
            var fuel = FuelModel.Chaparral;
            var str = fuel.ToString();
            Assert.IsTrue(str.Contains("Chaparral"), $"ToString should contain name: {str}");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  RothermelModel — Sub-equations
        // ═══════════════════════════════════════════════════════════════

        #region Rothermel Sub-equations

        [TestMethod]
        public void PackingRatio_ShortGrass_ReasonableValue()
        {
            double beta = RothermelModel.PackingRatio(FuelModel.ShortGrass);
            // β = (w₀/δ) / ρp = (0.166/0.305) / 513 ≈ 0.00106
            Assert.IsTrue(beta > 0 && beta < 0.1,
                $"Packing ratio should be small for grass: {beta}");
        }

        [TestMethod]
        public void OptimalPackingRatio_ShortGrass_ReasonableValue()
        {
            double betaOp = RothermelModel.OptimalPackingRatio(FuelModel.ShortGrass);
            // β_op = 3.348 · σ^(-0.8189), for σ=11483 → small value
            Assert.IsTrue(betaOp > 0 && betaOp < 0.1,
                $"Optimal packing ratio should be small: {betaOp}");
        }

        [TestMethod]
        public void EffectiveHeatingNumber_ShortGrass_CloseToOne()
        {
            double eps = RothermelModel.EffectiveHeatingNumber(FuelModel.ShortGrass);
            // ε = exp(-138/11483) ≈ 0.988
            Assert.AreEqual(0.988, eps, 0.005,
                "High σ gives ε close to 1.");
        }

        [TestMethod]
        public void HeatOfPreignition_DryFuel()
        {
            double Qig = RothermelModel.HeatOfPreignition(0.05);
            // Qig = 581 + 2594 · 0.05 = 710.7
            Assert.AreEqual(710.7, Qig, 0.1);
        }

        [TestMethod]
        public void HeatOfPreignition_ZeroMoisture()
        {
            double Qig = RothermelModel.HeatOfPreignition(0.0);
            Assert.AreEqual(581.0, Qig, 0.01);
        }

        [TestMethod]
        public void SlopeFactor_FlatTerrain_IsZero()
        {
            double beta = RothermelModel.PackingRatio(FuelModel.ShortGrass);
            double phiS = RothermelModel.SlopeFactor(beta, 0.0);
            Assert.AreEqual(0.0, phiS, 1e-10, "Flat terrain → φs = 0.");
        }

        [TestMethod]
        public void SlopeFactor_30Degrees_IsPositive()
        {
            double beta = RothermelModel.PackingRatio(FuelModel.ShortGrass);
            double slope30 = Math.PI / 6.0; // 30°
            double phiS = RothermelModel.SlopeFactor(beta, slope30);
            Assert.IsTrue(phiS > 0, $"30° slope should give φs > 0: {phiS}");
        }

        [TestMethod]
        public void SlopeFactor_IncreasesWithSlope()
        {
            double beta = RothermelModel.PackingRatio(FuelModel.ShortGrass);
            double phiS_20 = RothermelModel.SlopeFactor(beta, 20.0 * Math.PI / 180.0);
            double phiS_40 = RothermelModel.SlopeFactor(beta, 40.0 * Math.PI / 180.0);
            Assert.IsTrue(phiS_40 > phiS_20,
                $"Steeper slope should give larger φs: {phiS_20} < {phiS_40}");
        }

        [TestMethod]
        public void WindFactor_ZeroWind_IsZero()
        {
            double phiW = RothermelModel.WindFactor(FuelModel.ShortGrass, 0.0);
            Assert.AreEqual(0.0, phiW, 1e-10, "Zero wind → φw = 0.");
        }

        [TestMethod]
        public void WindFactor_Positive_IncreasesWithWind()
        {
            double phiW_low = RothermelModel.WindFactor(FuelModel.ShortGrass, 2.0);
            double phiW_high = RothermelModel.WindFactor(FuelModel.ShortGrass, 5.0);
            Assert.IsTrue(phiW_low > 0, "Wind factor should be positive.");
            Assert.IsTrue(phiW_high > phiW_low,
                $"More wind → larger φw: {phiW_low} < {phiW_high}");
        }

        [TestMethod]
        public void PropagatingFluxRatio_ShortGrass_Positive()
        {
            double xi = RothermelModel.PropagatingFluxRatio(FuelModel.ShortGrass);
            Assert.IsTrue(xi > 0 && xi < 1,
                $"Propagating flux ratio should be in (0,1): {xi}");
        }

        [TestMethod]
        public void ReactionIntensity_ShortGrass_DryFuel_Positive()
        {
            double IR = RothermelModel.ReactionIntensity(FuelModel.ShortGrass, 0.05);
            Assert.IsTrue(IR > 0,
                $"Reaction intensity should be positive for dry fuel: {IR}");
        }

        [TestMethod]
        public void ReactionIntensity_AtExtinction_IsZero()
        {
            var fuel = FuelModel.ShortGrass;
            double IR = RothermelModel.ReactionIntensity(fuel, fuel.MoistureOfExtinction);
            // At extinction moisture, ηM → 0 so IR should be ≈ 0
            Assert.AreEqual(0, IR, 1.0,
                "Reaction intensity at moisture extinction should be near zero.");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  RothermelModel — Rate of Spread
        // ═══════════════════════════════════════════════════════════════

        #region Rate of Spread

        [TestMethod]
        public void RateOfSpread_ShortGrass_DryFuel_Wind5mph_Flat()
        {
            // BehavePlus reference: Short Grass, M=0.05, wind 5 mph (≈2.24 m/s), flat
            // Expected R ≈ 23–25 m/min (BehavePlus validation range)
            var fuel = FuelModel.ShortGrass;
            double R = RothermelModel.RateOfSpread(fuel, 0.05, 2.24, 0.0);

            Assert.IsTrue(R > 5 && R < 80,
                $"Short Grass ROS with 5 mph wind should be in reasonable range: {R} m/min");
        }

        [TestMethod]
        public void RateOfSpread_ShortGrass_NoWind_Flat()
        {
            // No wind, flat → only base spread (no φw, no φs)
            var fuel = FuelModel.ShortGrass;
            double R = RothermelModel.RateOfSpread(fuel, 0.05, 0.0, 0.0);

            Assert.IsTrue(R > 0,
                $"ROS should be positive even without wind: {R} m/min");
        }

        [TestMethod]
        public void RateOfSpread_Chaparral_Moisture010_Wind10mph_Flat()
        {
            // Chaparral (model 4), M=0.10, wind 10 mph (≈4.47 m/s), flat
            var fuel = FuelModel.Chaparral;
            double R = RothermelModel.RateOfSpread(fuel, 0.10, 4.47, 0.0);

            Assert.IsTrue(R > 0,
                $"Chaparral ROS should be positive: {R} m/min");
        }

        [TestMethod]
        public void RateOfSpread_MoistureAtExtinction_IsZero()
        {
            var fuel = FuelModel.ShortGrass;
            double R = RothermelModel.RateOfSpread(fuel, fuel.MoistureOfExtinction, 2.24, 0.0);
            Assert.AreEqual(0, R, 1e-10, "ROS at moisture extinction should be 0.");
        }

        [TestMethod]
        public void RateOfSpread_MoistureAboveExtinction_IsZero()
        {
            var fuel = FuelModel.ShortGrass;
            double R = RothermelModel.RateOfSpread(fuel, 0.50, 2.24, 0.0);
            Assert.AreEqual(0, R, 1e-10, "ROS above moisture extinction should be 0.");
        }

        [TestMethod]
        public void RateOfSpread_WindIncreasesSpread()
        {
            var fuel = FuelModel.ShortGrass;
            double R_noWind = RothermelModel.RateOfSpread(fuel, 0.05, 0.0, 0.0);
            double R_wind = RothermelModel.RateOfSpread(fuel, 0.05, 3.0, 0.0);

            Assert.IsTrue(R_wind > R_noWind,
                $"Wind should increase ROS: no-wind={R_noWind}, wind={R_wind} m/min");
        }

        [TestMethod]
        public void RateOfSpread_SlopeIncreasesSpread()
        {
            var fuel = FuelModel.ShortGrass;
            double R_flat = RothermelModel.RateOfSpread(fuel, 0.05, 0.0, 0.0);
            double R_slope = RothermelModel.RateOfSpread(fuel, 0.05, 0.0, Math.PI / 6.0);

            Assert.IsTrue(R_slope > R_flat,
                $"Uphill slope should increase ROS: flat={R_flat}, slope={R_slope} m/min");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  RothermelModel — Flame Length
        // ═══════════════════════════════════════════════════════════════

        #region Flame Length

        [TestMethod]
        public void FlameLength_ZeroIntensity_IsZero()
        {
            Assert.AreEqual(0, RothermelModel.FlameLength(0, 10));
        }

        [TestMethod]
        public void FlameLength_ZeroSpread_IsZero()
        {
            Assert.AreEqual(0, RothermelModel.FlameLength(1000, 0));
        }

        [TestMethod]
        public void FlameLength_IncreasesWithIntensity()
        {
            double fl_low = RothermelModel.FlameLength(500, 5);
            double fl_high = RothermelModel.FlameLength(2000, 5);
            Assert.IsTrue(fl_high > fl_low,
                $"Higher IR should give longer flames: {fl_low} < {fl_high}");
        }

        [TestMethod]
        public void FlameLength_IncreasesWithSpreadRate()
        {
            double fl_slow = RothermelModel.FlameLength(1000, 2);
            double fl_fast = RothermelModel.FlameLength(1000, 20);
            Assert.IsTrue(fl_fast > fl_slow,
                $"Faster spread should give longer flames: {fl_slow} < {fl_fast}");
        }

        [TestMethod]
        public void FlameLength_ShortGrass_ReasonableValue()
        {
            var fuel = FuelModel.ShortGrass;
            double IR = RothermelModel.ReactionIntensity(fuel, 0.05);
            double R = RothermelModel.RateOfSpread(fuel, 0.05, 2.24, 0.0);
            double fl = RothermelModel.FlameLength(IR, R);

            // Short grass flame lengths typically 0.5–3 m
            Assert.IsTrue(fl > 0.1 && fl < 10,
                $"Short Grass flame length should be reasonable: {fl} m");
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  RothermelModel — Cross-model consistency
        // ═══════════════════════════════════════════════════════════════

        #region Cross-model Consistency

        [TestMethod]
        public void AllFuelModels_ProducePositiveROS_WhenDry()
        {
            foreach (var fuel in FuelLibrary.All())
            {
                double R = RothermelModel.RateOfSpread(fuel, 0.04, 2.0, 0.0);
                Assert.IsTrue(R > 0,
                    $"{fuel.Name} should produce positive ROS at low moisture: {R}");
            }
        }

        [TestMethod]
        public void AllFuelModels_ZeroROS_WhenSaturated()
        {
            foreach (var fuel in FuelLibrary.All())
            {
                double R = RothermelModel.RateOfSpread(fuel, 0.99, 2.0, 0.0);
                Assert.AreEqual(0, R, 1e-10,
                    $"{fuel.Name} should have zero ROS when saturated.");
            }
        }

        #endregion
    }
}
