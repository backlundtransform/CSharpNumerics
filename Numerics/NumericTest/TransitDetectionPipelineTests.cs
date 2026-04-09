using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Pipeline;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class TransitDetectionPipelineTests
{
    /// <summary>
    /// Creates a synthetic light curve with injected transit dips.
    /// </summary>
    private static LightCurve CreateSyntheticLightCurve(
        double baseFlux, double noise, double[] periods, double[] depths,
        double[] epochs, double[] durations, int numPoints, double cadenceDays,
        int seed = 42)
    {
        var rng = new Random(seed);
        var time = new double[numPoints];
        var flux = new double[numPoints];
        var err = new double[numPoints];
        var qual = new int[numPoints];

        for (int i = 0; i < numPoints; i++)
        {
            time[i] = i * cadenceDays;
            flux[i] = baseFlux + noise * (rng.NextDouble() * 2 - 1);
            err[i] = noise;
        }

        // Inject transits
        for (int p = 0; p < periods.Length; p++)
        {
            double period = periods[p];
            double depth = depths[p];
            double epoch = epochs[p];
            double duration = durations[p];
            double halfDur = duration / 2.0;

            for (int i = 0; i < numPoints; i++)
            {
                double phase = ((time[i] - epoch) % period + period) % period;
                if (phase > period / 2.0) phase -= period;

                if (Math.Abs(phase) <= halfDur)
                {
                    flux[i] -= depth * baseFlux;
                }
            }
        }

        var meta = new LightCurveMetadata("SYNTHETIC", "SIM", 0, CadenceType.Long);
        return new LightCurve(time, flux, err, qual, meta);
    }

    #region LightCurvePreprocessor Tests

    [TestMethod]
    public void Preprocessor_ShouldNormalizeAndClean()
    {
        var lc = CreateSyntheticLightCurve(1000.0, 1.0, new double[0], new double[0],
            new double[0], new double[0], 500, 0.02);

        var config = new TransitDetectionConfig();
        var processed = LightCurvePreprocessor.Preprocess(lc, config);

        // Processed flux should be close to 1.0 (normalized)
        double mean = processed.Flux.Average();
        Assert.IsTrue(Math.Abs(mean - 1.0) < 0.1,
            $"Preprocessed mean flux {mean} should be close to 1.0");
    }

    [TestMethod]
    public void Preprocessor_Detrend_ShouldRemoveTrend()
    {
        int n = 500;
        var time = new double[n];
        var flux = new double[n];
        var err = new double[n];
        var qual = new int[n];

        // Add a linear trend
        for (int i = 0; i < n; i++)
        {
            time[i] = i * 0.02;
            flux[i] = 1000.0 + 0.5 * i; // strong upward trend
            err[i] = 1.0;
        }

        var lc = new LightCurve(time, flux, err, qual, new LightCurveMetadata());
        var config = new TransitDetectionConfig
        {
            DetrendingMethod = DetrendingMethod.Polynomial,
            DetrendingPolyDegree = 1
        };

        var detrended = LightCurvePreprocessor.Detrend(lc, config);

        // Detrended values should have much less range than original
        double range = detrended.Flux.Max() - detrended.Flux.Min();
        double origRange = flux.Max() - flux.Min();

        Assert.IsTrue(range < origRange * 0.1,
            $"Detrended range ({range}) should be much smaller than original ({origRange})");
    }

    #endregion

    #region PeriodSearcher Tests

    [TestMethod]
    public void PeriodSearcher_BLS_ShouldFindInjectedPeriod()
    {
        double injectedPeriod = 5.0;
        double depth = 0.01;
        double duration = 0.15;

        var lc = CreateSyntheticLightCurve(1.0, 0.001, new[] { injectedPeriod }, new[] { depth },
            new[] { 0.0 }, new[] { duration }, 2000, 0.02);

        var result = PeriodSearcher.Search(lc, 1.0, 20.0, PeriodSearchMethod.BLS, 3000);

        double periodError = Math.Abs(result.BestPeriod - injectedPeriod) / injectedPeriod;
        Assert.IsTrue(periodError < 0.02,
            $"BLS found period {result.BestPeriod}, expected ~{injectedPeriod} (error={periodError:P1})");
    }

    [TestMethod]
    public void PeriodSearcher_LombScargle_ShouldFindSignal()
    {
        // Create a sinusoidal signal (LS is optimal for sinusoids)
        int n = 1000;
        double period = 3.0;
        var time = new double[n];
        var flux = new double[n];
        var err = new double[n];
        var qual = new int[n];
        var rng = new Random(42);

        for (int i = 0; i < n; i++)
        {
            time[i] = i * 0.02;
            flux[i] = 1.0 + 0.01 * Math.Sin(2 * Math.PI * time[i] / period) + 0.001 * (rng.NextDouble() * 2 - 1);
            err[i] = 0.001;
        }

        var lc = new LightCurve(time, flux, err, qual, new LightCurveMetadata());
        var result = PeriodSearcher.Search(lc, 1.0, 10.0, PeriodSearchMethod.LombScargle, 2000);

        double periodError = Math.Abs(result.BestPeriod - period) / period;
        Assert.IsTrue(periodError < 0.05,
            $"LS found period {result.BestPeriod}, expected ~{period} (error={periodError:P1})");
    }

    #endregion

    #region TransitFitter Tests

    [TestMethod]
    public void TransitFitter_ShouldRecoverTransitDepth()
    {
        double period = 5.0;
        double depth = 0.01;
        double duration = 0.15;
        double epoch = 0.0;

        var lc = CreateSyntheticLightCurve(1.0, 0.0005, new[] { period }, new[] { depth },
            new[] { epoch }, new[] { duration }, 3000, 0.02, seed: 123);

        var result = TransitFitter.Fit(lc, period, epoch);

        Assert.IsTrue(Math.Abs(result.Parameters.Depth - depth) < depth * 0.5,
            $"Fitted depth {result.Parameters.Depth} should be close to {depth}");
        Assert.IsTrue(result.Parameters.Period == period);
    }

    #endregion

    #region TransitValidator Tests

    [TestMethod]
    public void Validator_GoodCandidate_ShouldBeValid()
    {
        double period = 5.0;
        double depth = 0.01;
        var lc = CreateSyntheticLightCurve(1.0, 0.0005, new[] { period }, new[] { depth },
            new[] { 0.0 }, new[] { 0.15 }, 3000, 0.02);

        var parameters = new TransitParameters(period, 0.0, depth, 0.15, 0.1, 0.0, 0.02);
        var phaseFolded = LightCurve.FromArrays(new[] { -0.1, 0.0, 0.1 }, new[] { 1.0, 0.99, 1.0 });
        var candidate = new TransitCandidate(parameters, 0.9, TransitDisposition.Unknown, phaseFolded, new TransitFeatureSet());

        var config = new TransitDetectionConfig { SnrThreshold = 3.0, MinTransitDepthPpm = 50 };
        var result = TransitValidator.Validate(candidate, lc, config);

        Assert.IsTrue(result.IsValid, $"Good candidate should validate. Warnings: {string.Join("; ", result.Warnings)}");
        Assert.IsTrue(result.Score > 0, $"Score should be positive: {result.Score}");
        Assert.IsTrue(result.Snr > 0, $"SNR should be positive: {result.Snr}");
    }

    [TestMethod]
    public void Validator_ShallowTransit_ShouldFail()
    {
        var lc = CreateSyntheticLightCurve(1.0, 0.001, new double[0], new double[0],
            new double[0], new double[0], 1000, 0.02);

        var parameters = new TransitParameters(5.0, 0.0, 0.00001, 0.1, 0.003, 0.0, 0.01);
        var phaseFolded = LightCurve.FromArrays(new[] { 0.0 }, new[] { 1.0 });
        var candidate = new TransitCandidate(parameters, 0.5, TransitDisposition.Unknown, phaseFolded, new TransitFeatureSet());

        var config = new TransitDetectionConfig { MinTransitDepthPpm = 100, SnrThreshold = 7.0 };
        var result = TransitValidator.Validate(candidate, lc, config);

        Assert.IsFalse(result.IsValid, "Shallow transit should not validate.");
    }

    #endregion

    #region Full Pipeline Tests

    [TestMethod]
    public void Pipeline_SinglePlanet_ShouldDetectTransit()
    {
        double period = 5.0;
        double depth = 0.015;
        double duration = 0.15;

        var lc = CreateSyntheticLightCurve(1.0, 0.0003, new[] { period }, new[] { depth },
            new[] { 0.0 }, new[] { duration }, 5000, 0.02, seed: 99);

        var config = new TransitDetectionConfig
        {
            MinPeriodDays = 2.0,
            MaxPeriodDays = 15.0,
            MinTransitDepthPpm = 50,
            SnrThreshold = 3.0,
            MaxPlanets = 1,
            PeriodSearchMethod = PeriodSearchMethod.BLS,
            DetrendingMethod = DetrendingMethod.MedianFilter,
            DetrendingWindowSize = 101,
            NumTrialPeriods = 3000
        };

        // Step 1: Preprocess
        var processed = LightCurvePreprocessor.Preprocess(lc, config);
        Assert.IsTrue(processed.Length > 100, $"Preprocessed LC has {processed.Length} points");

        // Step 2: BLS
        var searchResult = PeriodSearcher.Search(processed, config);
        double periodError = Math.Abs(searchResult.BestPeriod - period) / period;
        Assert.IsTrue(periodError < 0.05,
            $"BLS should find period ~{period}, got {searchResult.BestPeriod}");

        // Step 3: Fit
        var fitResult = TransitFitter.Fit(processed, searchResult.BestPeriod, searchResult.BestEpoch);
        Assert.IsTrue(fitResult.Parameters.Depth > 0.001,
            $"Fitted depth should be > 0.001, got {fitResult.Parameters.Depth}");

        // Step 4: Validate
        var phaseFolded = LightCurve.FromArrays(new[] { 0.0 }, new[] { 1.0 });
        var candidate = new TransitCandidate(
            fitResult.Parameters, 0.0, TransitDisposition.Unknown, phaseFolded, new TransitFeatureSet());
        var validation = TransitValidator.Validate(candidate, processed, config);
        Assert.IsTrue(validation.IsValid,
            $"Candidate should validate. SNR={validation.Snr:F2}, Score={validation.Score:F2}, Warnings: {string.Join("; ", validation.Warnings)}");

        // Step 5: Full pipeline
        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.IsTrue(candidates.Length >= 1, "Should detect at least 1 planet.");

        double detectedPeriodError = Math.Abs(candidates[0].Parameters.Period - period) / period;
        Assert.IsTrue(detectedPeriodError < 0.02,
            $"Detected period {candidates[0].Parameters.Period:F4} should be within 2% of {period} (error={detectedPeriodError:P1})");
    }

    [TestMethod]
    public void Pipeline_NoPlanet_ShouldDetectNothing()
    {
        // Pure noise, no transit
        var lc = CreateSyntheticLightCurve(1.0, 0.001, new double[0], new double[0],
            new double[0], new double[0], 2000, 0.02);

        var config = new TransitDetectionConfig
        {
            MinPeriodDays = 1.0,
            MaxPeriodDays = 20.0,
            MinTransitDepthPpm = 500,
            SnrThreshold = 10.0,
            MaxPlanets = 3,
            NumTrialPeriods = 1000
        };

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.AreEqual(0, candidates.Length, "Should detect 0 planets in pure noise with high threshold.");
    }

    [TestMethod]
    public void Pipeline_TwoPlanets_ShouldDetectBoth()
    {
        double p1 = 5.0, p2 = 8.0;
        double d1 = 0.012, d2 = 0.008;
        double dur1 = 0.15, dur2 = 0.18;

        var lc = CreateSyntheticLightCurve(1.0, 0.0003, new[] { p1, p2 }, new[] { d1, d2 },
            new[] { 0.0, 1.5 }, new[] { dur1, dur2 }, 5000, 0.02, seed: 77);

        var config = new TransitDetectionConfig
        {
            MinPeriodDays = 2.0,
            MaxPeriodDays = 20.0,
            MinTransitDepthPpm = 50,
            SnrThreshold = 3.0,
            MaxPlanets = 3,
            PeriodSearchMethod = PeriodSearchMethod.BLS,
            NumTrialPeriods = 3000
        };

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.IsTrue(candidates.Length >= 2,
            $"Should detect at least 2 planets, found {candidates.Length}.");

        // Verify both periods are recovered (order may vary)
        var detectedPeriods = candidates.Select(c => c.Parameters.Period).ToArray();

        bool foundP1 = detectedPeriods.Any(p => Math.Abs(p - p1) / p1 < 0.02);
        bool foundP2 = detectedPeriods.Any(p => Math.Abs(p - p2) / p2 < 0.02);

        Assert.IsTrue(foundP1, $"Should detect planet with period ~{p1}. Detected: {string.Join(", ", detectedPeriods.Select(p => p.ToString("F3")))}");
        Assert.IsTrue(foundP2, $"Should detect planet with period ~{p2}. Detected: {string.Join(", ", detectedPeriods.Select(p => p.ToString("F3")))}");
    }

    #endregion
}
