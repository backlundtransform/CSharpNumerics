using CSharpNumerics.Physics.Environmental.Enums;
using CSharpNumerics.Physics.Environmental.Fire;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest
{
    [TestClass]
    public class CompartmentFireTests
    {
        // ════════════════════════════════════════════
        //  Heat release rate — t-squared growth
        // ════════════════════════════════════════════

        [TestMethod]
        public void HeatReleaseRate_GrowthClassesMatchTheStandardCoefficients()
        {
            // The published values: 0.00293, 0.01172, 0.0469 and 0.1876 kW/s².
            Assert.AreEqual(0.002930, HeatReleaseRate.GrowthCoefficient(FireGrowthRate.Slow), 1e-6);
            Assert.AreEqual(0.011722, HeatReleaseRate.GrowthCoefficient(FireGrowthRate.Medium), 1e-6);
            Assert.AreEqual(0.046889, HeatReleaseRate.GrowthCoefficient(FireGrowthRate.Fast), 1e-6);
            Assert.AreEqual(0.187556, HeatReleaseRate.GrowthCoefficient(FireGrowthRate.UltraFast), 1e-6);
        }

        [TestMethod]
        public void HeatReleaseRate_EachClassReachesTheReferenceAtItsDefiningTime()
        {
            foreach (FireGrowthRate rate in Enum.GetValues(typeof(FireGrowthRate)))
            {
                double alpha = HeatReleaseRate.GrowthCoefficient(rate);
                double t = HeatReleaseRate.TimeToReference(rate);

                Assert.AreEqual(1055.0, HeatReleaseRate.Growth(alpha, t), 1e-9,
                    rate + " should reach 1055 kW at " + t + " s");
                Assert.AreEqual(t, HeatReleaseRate.TimeToReach(alpha, 1055.0), 1e-9);
            }
        }

        [TestMethod]
        public void HeatReleaseRate_GrowthIsQuadraticAndZeroBeforeIgnition()
        {
            double alpha = HeatReleaseRate.GrowthCoefficient(FireGrowthRate.Medium);

            Assert.AreEqual(0, HeatReleaseRate.Growth(alpha, -10), 1e-12);
            Assert.AreEqual(0, HeatReleaseRate.Growth(alpha, 0), 1e-12);

            // Doubling the time quadruples the heat release rate.
            double at60 = HeatReleaseRate.Growth(alpha, 60);
            double at120 = HeatReleaseRate.Growth(alpha, 120);
            Assert.AreEqual(4.0, at120 / at60, 1e-9);
        }

        [TestMethod]
        public void HeatReleaseRate_CurveRunsGrowthThenSteadyThenDecay()
        {
            double alpha = HeatReleaseRate.GrowthCoefficient(FireGrowthRate.Fast);
            const double peak = 2000;
            const double steady = 300;
            const double decay = 200;

            double growthEnd = HeatReleaseRate.TimeToReach(alpha, peak);

            Assert.AreEqual(0, HeatReleaseRate.Curve(alpha, peak, steady, decay, 0), 1e-12);

            // Still growing halfway to the peak.
            double half = HeatReleaseRate.Curve(alpha, peak, steady, decay, growthEnd / 2);
            Assert.IsTrue(half > 0 && half < peak);

            // At the peak and through the steady period.
            Assert.AreEqual(peak, HeatReleaseRate.Curve(alpha, peak, steady, decay, growthEnd), 1e-6);
            Assert.AreEqual(peak, HeatReleaseRate.Curve(alpha, peak, steady, decay, growthEnd + steady), 1e-6);

            // Halfway through decay, and burnt out after it.
            Assert.AreEqual(peak / 2, HeatReleaseRate.Curve(alpha, peak, steady, decay, growthEnd + steady + decay / 2), 1e-6);
            Assert.AreEqual(0, HeatReleaseRate.Curve(alpha, peak, steady, decay, growthEnd + steady + decay), 1e-12);
            Assert.AreEqual(0, HeatReleaseRate.Curve(alpha, peak, steady, decay, 1e6), 1e-12);
        }

        [TestMethod]
        public void HeatReleaseRate_ConvectivePartIsAFractionOfTheTotal()
        {
            Assert.AreEqual(700, HeatReleaseRate.Convective(1000), 1e-9);
            Assert.AreEqual(600, HeatReleaseRate.Convective(1000, 0.6), 1e-9);
            Assert.AreEqual(0, HeatReleaseRate.Convective(0), 1e-12);
        }

        // ════════════════════════════════════════════
        //  Heskestad plume
        // ════════════════════════════════════════════

        [TestMethod]
        public void FirePlume_FlameHeightMatchesTheCorrelation()
        {
            // L = 0.235 Q^(2/5) − 1.02 D, so a 1 MW fire 1 m across stands 2.70 m.
            Assert.AreEqual(2.7045, FirePlume.FlameHeight(1000, 1.0), 1e-3);

            // A wide, weak fire has no coherent flame.
            Assert.AreEqual(0, FirePlume.FlameHeight(50, 5.0), 1e-12);
            Assert.AreEqual(0, FirePlume.FlameHeight(0, 1.0), 1e-12);
        }

        [TestMethod]
        public void FirePlume_FlameHeightGrowsWithFireSize()
        {
            double small = FirePlume.FlameHeight(500, 1.0);
            double large = FirePlume.FlameHeight(5000, 1.0);
            Assert.IsTrue(large > small);
        }

        [TestMethod]
        public void FirePlume_VirtualOriginMatchesTheCorrelation()
        {
            Assert.AreEqual(0.2955, FirePlume.VirtualOrigin(1000, 1.0), 1e-3);
        }

        [TestMethod]
        public void FirePlume_CentrelineTemperatureReducesToTheStandardAirForm()
        {
            // For standard air the leading factor collapses to about 25, so
            // ΔT ≈ 25 · Qc^(2/3) · z^(−5/3).
            double expected = 25.0 * Math.Pow(700, 2.0 / 3.0) * Math.Pow(3.0, -5.0 / 3.0);
            double actual = FirePlume.CentrelineTemperatureRise(700, 3.0, 0.0);

            Assert.AreEqual(expected, actual, expected * 0.002);
        }

        [TestMethod]
        public void FirePlume_CentrelineTemperatureFallsWithHeight()
        {
            double low = FirePlume.CentrelineTemperatureRise(700, 2.0, 0.0);
            double high = FirePlume.CentrelineTemperatureRise(700, 6.0, 0.0);

            Assert.IsTrue(low > high, "the plume cools as it rises and entrains air");
            Assert.AreEqual(0, FirePlume.CentrelineTemperatureRise(700, 0.5, 1.0), 1e-12);
        }

        [TestMethod]
        public void FirePlume_CentrelineVelocityReducesToTheStandardAirForm()
        {
            // u ≈ 1.03 · (Qc / z)^(1/3)
            double expected = 1.031 * Math.Pow(700 / 3.0, 1.0 / 3.0);
            double actual = FirePlume.CentrelineVelocity(700, 3.0, 0.0);

            Assert.AreEqual(expected, actual, expected * 0.002);
        }

        [TestMethod]
        public void FirePlume_CentrelineVelocityIsPlausibleForACabinFire()
        {
            // A 1 MW fire, 3 m up: metres per second, not centimetres and not
            // tens of metres. This is the number a buoyancy model consumes, so
            // an order-of-magnitude error here would be invisible but ruinous.
            double u = FirePlume.CentrelineVelocity(700, 3.0, 0.0);

            Assert.IsTrue(u > 3 && u < 12, "implausible plume velocity: " + u + " m/s");
        }

        [TestMethod]
        public void FirePlume_EntrainsFarMoreAirThanTheFireProduces()
        {
            // The point of the correlation: a compartment fills with smoke much
            // faster than the fire makes combustion products, because the plume
            // drags room air up with it.
            double atThreeMetres = FirePlume.MassFlowRate(700, 3.0, 0.0);
            double atSixMetres = FirePlume.MassFlowRate(700, 6.0, 0.0);

            Assert.IsTrue(atThreeMetres > 1.0, "a 1 MW plume moves kilograms per second");
            Assert.IsTrue(atSixMetres > atThreeMetres, "entrainment grows with height");
            Assert.AreEqual(0, FirePlume.MassFlowRate(700, 0.0, 1.0), 1e-12);
        }

        // ════════════════════════════════════════════
        //  Fractional effective dose
        // ════════════════════════════════════════════

        [TestMethod]
        public void Fed_CarbonMonoxideMatchesTheStandardSimplification()
        {
            // At the default breathing rate and incapacitating COHb the rate is
            // [CO]^1.036 / 36180 per minute.
            double expected = Math.Pow(1000, 1.036) / 36180.0;
            Assert.AreEqual(expected, FractionalEffectiveDose.CarbonMonoxideRate(1000), expected * 0.005);

            Assert.AreEqual(0, FractionalEffectiveDose.CarbonMonoxideRate(0), 1e-12);
        }

        [TestMethod]
        public void Fed_CarbonMonoxideAtOneThousandPpmIncapacitatesInAboutHalfAnHour()
        {
            double minutes = FractionalEffectiveDose.TimeToIncapacitation(1000);
            Assert.IsTrue(minutes > 20 && minutes < 40,
                "1000 ppm CO should incapacitate in roughly half an hour, got " + minutes + " min");
        }

        [TestMethod]
        public void Fed_CyanideActsFarFasterThanCarbonMonoxide()
        {
            // 200 ppm HCN is rapidly incapacitating, while 200 ppm CO is not.
            double hcn = FractionalEffectiveDose.TimeToIncapacitation(0, 200);
            double co = FractionalEffectiveDose.TimeToIncapacitation(200);

            Assert.IsTrue(hcn < 5, "200 ppm HCN should incapacitate within minutes, got " + hcn);
            Assert.IsTrue(co > hcn * 10, "200 ppm CO should take far longer than the same HCN");
        }

        [TestMethod]
        public void Fed_CarbonDioxideMakesEverythingElseWorse()
        {
            // 5 % CO2 roughly triples uptake by driving breathing rate up.
            double factor = FractionalEffectiveDose.HyperventilationFactor(5.0);
            Assert.IsTrue(factor > 2.4 && factor < 3.0, "unexpected hyperventilation factor: " + factor);

            Assert.AreEqual(1.0, FractionalEffectiveDose.HyperventilationFactor(0), 1e-12);

            double alone = FractionalEffectiveDose.Rate(500);
            double withCo2 = FractionalEffectiveDose.Rate(500, 0, 5.0);
            Assert.IsTrue(withCo2 > alone);
        }

        [TestMethod]
        public void Fed_OxygenDepletionOnlyCountsBelowAmbient()
        {
            Assert.AreEqual(0, FractionalEffectiveDose.OxygenDepletionRate(20.9), 1e-12);
            Assert.AreEqual(0, FractionalEffectiveDose.OxygenDepletionRate(21.5), 1e-12);
            Assert.IsTrue(FractionalEffectiveDose.OxygenDepletionRate(12.0) > 0);
        }

        [TestMethod]
        public void Fed_CleanAirNeverIncapacitates()
        {
            Assert.AreEqual(0, FractionalEffectiveDose.Rate(0), 1e-12);
            Assert.IsTrue(double.IsPositiveInfinity(FractionalEffectiveDose.TimeToIncapacitation(0)));
        }

        [TestMethod]
        public void Fed_DoseAccumulatesOverTime()
        {
            // Held at a constant exposure, the dose after the predicted time to
            // incapacitation must be 1.0.
            double rate = FractionalEffectiveDose.Rate(1000);
            double minutes = FractionalEffectiveDose.TimeToIncapacitation(1000);

            const double stepSeconds = 1.0;
            int steps = (int)Math.Round(minutes * 60.0 / stepSeconds);
            var history = new double[steps];
            for (int i = 0; i < steps; i++) history[i] = rate;

            Assert.AreEqual(1.0, FractionalEffectiveDose.Accumulate(history, stepSeconds), 0.01);
        }

        [TestMethod]
        public void Fed_HeatDoseMatchesPublishedTimes()
        {
            // Purser's convective-heat correlation for a lightly clothed person:
            // about an hour at 40 °C, about twelve minutes at 100 °C.
            double at40 = 1.0 / FractionalEffectiveDose.ConvectiveHeatRate(40);
            double at100 = 1.0 / FractionalEffectiveDose.ConvectiveHeatRate(100);

            Assert.IsTrue(at40 > 50 && at40 < 70, "40 °C should take about an hour, got " + at40 + " min");
            Assert.IsTrue(at100 > 9 && at100 < 14, "100 °C should take about twelve minutes, got " + at100 + " min");
        }

        [TestMethod]
        public void Fed_OrdinaryWarmthIsNotADose()
        {
            Assert.AreEqual(0, FractionalEffectiveDose.ConvectiveHeatRate(20), 1e-12);
            Assert.AreEqual(0, FractionalEffectiveDose.ConvectiveHeatRate(30), 1e-12);
            Assert.IsTrue(FractionalEffectiveDose.ConvectiveHeatRate(31) > 0);
        }

        // ════════════════════════════════════════════
        //  Visibility
        // ════════════════════════════════════════════

        [TestMethod]
        public void Visibility_FallsAsSmokeThickens()
        {
            double thin = FractionalEffectiveDose.Visibility(1e-5);
            double thick = FractionalEffectiveDose.Visibility(1e-3);

            Assert.IsTrue(thin > thick);
            Assert.IsTrue(double.IsPositiveInfinity(FractionalEffectiveDose.Visibility(0)));
        }

        [TestMethod]
        public void Visibility_TenMetreLimitIsReachedAtThinSmoke()
        {
            // S = C / (Km · ρ). With C = 3 and Km = 7600, ten metres visibility
            // corresponds to about 39 mg of soot per cubic metre — smoke thin
            // enough that the toxic dose is still negligible. This is why
            // visibility, not toxicity, is what usually stops escape.
            double concentration = 3.0 / (7600.0 * 10.0);

            Assert.AreEqual(10.0, FractionalEffectiveDose.Visibility(concentration), 1e-9);
            Assert.IsTrue(concentration < 5e-5);
        }
    }
}
