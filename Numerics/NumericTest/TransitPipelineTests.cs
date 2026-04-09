using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Pipeline;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class TransitPipelineTests
{
    #region Synthetic light curve helpers

    /// <summary>
    /// Generates a synthetic light curve with one or more injected transits.
    /// Uses a simple box transit model (flat bottom) for testability.
    /// </summary>
    private static LightCurve GenerateSyntheticLightCurve(
        double totalDays, double cadenceMinutes,
        TransitSignal[] signals, double noiseLevel = 0.0003, int seed = 42)
    {
        int n = (int)(totalDays * 24 * 60 / cadenceMinutes);
        double[] time = new double[n];
        double[] flux = new double[n];
        double[] err = new double[n];
        int[] qual = new int[n];

        var rng = new Random(seed);
        double cadenceDays = cadenceMinutes / (24.0 * 60.0);

        for (int i = 0; i < n; i++)
        {
            time[i] = i * cadenceDays;
            flux[i] = 1.0;
            err[i] = noiseLevel;
        }

        // Inject transit signals
        foreach (var sig in signals)
        {
            for (int i = 0; i < n; i++)
            {
                double phase = ((time[i] - sig.Epoch) % sig.Period + sig.Period) % sig.Period;
                if (phase > sig.Period / 2.0) phase -= sig.Period;
                double absPhase = Math.Abs(phase);
                double halfDur = sig.Duration / 2.0;

                if (absPhase < halfDur)
                {
                    flux[i] -= sig.Depth;
                }
            }
        }

        // Add Gaussian noise
        for (int i = 0; i < n; i++)
        {
            double u1 = 1.0 - rng.NextDouble();
            double u2 = rng.NextDouble();
            double noise = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            flux[i] += noise * noiseLevel;
        }

        return new LightCurve(time, flux, err, qual,
            new LightCurveMetadata("SYNTH-001", "Synthetic", 0.0, CadenceType.Long));
    }

    /// <summary>
    /// Generates a light curve with a linear trend added on top of transits.
    /// </summary>
    private static LightCurve GenerateWithTrend(
        double totalDays, double cadenceMinutes,
        TransitSignal[] signals, double trendSlope, double noiseLevel = 0.0003, int seed = 42)
    {
        var lc = GenerateSyntheticLightCurve(totalDays, cadenceMinutes, signals, noiseLevel, seed);

        double[] flux = new double[lc.Length];
        for (int i = 0; i < lc.Length; i++)
        {
            double trendFactor = 1.0 + trendSlope * lc.Time[i];
            flux[i] = lc.Flux[i] * trendFactor;
        }

        return new LightCurve((double[])lc.Time.Clone(), flux,
            (double[])lc.FluxError.Clone(), (int[])lc.QualityFlags.Clone(), lc.Metadata);
    }

    private class TransitSignal
    {
        public double Period { get; set; }
        public double Epoch { get; set; }
        public double Depth { get; set; }
        public double Duration { get; set; }
    }

    private static TransitDetectionConfig DefaultConfig(double maxPeriod = 20.0) => new TransitDetectionConfig
    {
        MinPeriodDays = 1.0,
        MaxPeriodDays = maxPeriod,
        MinTransitDepthPpm = 500,
        SnrThreshold = 5.0,
        MaxPlanets = 3,
        DetrendingMethod = DetrendingMethod.MedianFilter,
        PeriodSearchMethod = PeriodSearchMethod.BLS,
        DetrendingWindowSize = 101,
        DetrendingPolyDegree = 3,
        OutlierSigmaThreshold = 5.0,
        NumTrialPeriods = 3000,
        PhaseBins = 200
    };

    #endregion

    #region LightCurvePreprocessor Tests

    [TestMethod]
    public void Preprocessor_ShouldPreserveFlatLightCurve()
    {
        var lc = GenerateSyntheticLightCurve(50, 30, new TransitSignal[0], 0.0001);
        var config = DefaultConfig();

        var result = LightCurvePreprocessor.Preprocess(lc, config);

        Assert.IsTrue(result.Length > 0);
        double mean = result.Flux.Average();
        Assert.AreEqual(1.0, mean, 0.01, "Preprocessed flat curve should be near 1.0");
    }

    [TestMethod]
    public void Preprocessor_ShouldRemoveLinearTrend()
    {
        var signals = new[] { new TransitSignal { Period = 5.0, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateWithTrend(50, 30, signals, 0.001, 0.0003);
        var config = DefaultConfig();
        config.DetrendingMethod = DetrendingMethod.Polynomial;
        config.DetrendingPolyDegree = 1;

        var result = LightCurvePreprocessor.Preprocess(lc, config);

        // After detrending + normalization, out-of-transit flux should be near 1.0
        double mean = result.Flux.Where(f => f > 0.99).Average();
        Assert.AreEqual(1.0, mean, 0.01, $"Detrended mean = {mean}");
    }

    [TestMethod]
    public void Preprocessor_ShouldPreserveTransitDips()
    {
        var signals = new[] { new TransitSignal { Period = 5.0, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(50, 30, signals, 0.0003);
        var config = DefaultConfig();

        var result = LightCurvePreprocessor.Preprocess(lc, config);

        double minFlux = result.Flux.Min();
        Assert.IsTrue(minFlux < 0.995, $"Transit dip should be preserved, min flux = {minFlux}");
    }

    #endregion

    #region PeriodSearcher Tests

    [TestMethod]
    public void PeriodSearcher_BLS_ShouldFindCorrectPeriod()
    {
        double truePeriod = 5.0;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);

        var result = PeriodSearcher.Search(lc, 1.0, 15.0, PeriodSearchMethod.BLS, 3000);

        double periodError = Math.Abs(result.BestPeriod - truePeriod) / truePeriod;
        Assert.IsTrue(periodError < 0.02,
            $"BLS period = {result.BestPeriod:F4}, expected {truePeriod}, error = {periodError:P1}");
        Assert.IsTrue(result.BestPower > 0, "BLS power should be positive");
        Assert.IsTrue(result.BestDepth > 0, "BLS depth should be positive");
    }

    [TestMethod]
    public void PeriodSearcher_LombScargle_ShouldFindPeriod()
    {
        double truePeriod = 5.0;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);

        var result = PeriodSearcher.Search(lc, 1.0, 15.0, PeriodSearchMethod.LombScargle, 3000);

        // Lomb-Scargle may find the period or a harmonic — just verify it runs cleanly
        Assert.IsTrue(result.BestPower > 0, "LS power should be positive");
        Assert.IsTrue(result.BestPeriod > 0, "LS period should be positive");
    }

    [TestMethod]
    public void PeriodSearcher_Both_ShouldReturnResult()
    {
        double truePeriod = 5.0;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);

        var result = PeriodSearcher.Search(lc, 1.0, 15.0, PeriodSearchMethod.Both, 3000);

        Assert.IsTrue(result.BestPower > 0);
        Assert.IsTrue(result.BestPeriod > 0);
    }

    #endregion

    #region TransitFitter Tests

    [TestMethod]
    public void TransitFitter_ShouldRecoverDepth()
    {
        double truePeriod = 5.0;
        double trueDepth = 0.01;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = trueDepth, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);

        var result = TransitFitter.Fit(lc, truePeriod, 1.0);

        double depthError = Math.Abs(result.Parameters.Depth - trueDepth) / trueDepth;
        Assert.IsTrue(depthError < 0.3,
            $"Fitted depth = {result.Parameters.Depth:F5}, expected {trueDepth}, error = {depthError:P0}");
        Assert.IsTrue(result.Parameters.Duration > 0, "Duration should be positive");
        Assert.AreEqual(truePeriod, result.Parameters.Period, 1e-10);
    }

    [TestMethod]
    public void TransitFitter_ShouldRecoverDuration()
    {
        double truePeriod = 5.0;
        double trueDuration = 0.1;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = 0.01, Duration = trueDuration } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);

        var result = TransitFitter.Fit(lc, truePeriod, 1.0);

        double durationError = Math.Abs(result.Parameters.Duration - trueDuration) / trueDuration;
        Assert.IsTrue(durationError < 0.5,
            $"Fitted duration = {result.Parameters.Duration:F4} d, expected {trueDuration}, error = {durationError:P0}");
    }

    [TestMethod]
    public void TransitFitter_ShouldReturnQualityMetrics()
    {
        var signals = new[] { new TransitSignal { Period = 5.0, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);

        var result = TransitFitter.Fit(lc, 5.0, 1.0);

        Assert.IsTrue(result.ChiSquared >= 0, "Chi-squared should be non-negative");
        Assert.IsTrue(!double.IsNaN(result.BIC), "BIC should not be NaN");
        Assert.IsTrue(result.ReducedChiSquared >= 0, "Reduced chi-squared should be non-negative");
        Assert.IsTrue(result.FittedFlux.Length > 0, "Fitted flux should have data");
        Assert.IsTrue(result.Residuals.Length > 0, "Residuals should have data");
    }

    #endregion

    #region TransitValidator Tests

    [TestMethod]
    public void TransitValidator_GoodCandidate_ShouldBeValid()
    {
        double truePeriod = 5.0;
        double trueDepth = 0.01;
        double trueDuration = 0.1;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = trueDepth, Duration = trueDuration } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);
        var config = DefaultConfig();

        var parameters = new TransitParameters(truePeriod, 1.0, trueDepth, trueDuration, 0.1, 0.0, 0.015);
        var candidate = new TransitCandidate(
            parameters, 0.0, TransitDisposition.Unknown,
            LightCurve.FromArrays(new double[] { 0 }, new double[] { 1 }),
            new TransitFeatureSet());

        var result = TransitValidator.Validate(candidate, lc, config);

        Assert.IsTrue(result.IsValid, $"Good candidate should be valid. Warnings: {string.Join("; ", result.Warnings)}");
        Assert.IsTrue(result.Snr > config.SnrThreshold,
            $"SNR ({result.Snr:F1}) should exceed threshold ({config.SnrThreshold})");
        Assert.IsTrue(result.Score > 0.0, $"Score should be > 0, got {result.Score:F3}");
    }

    [TestMethod]
    public void TransitValidator_ShallowTransit_ShouldFailDepthCheck()
    {
        var lc = GenerateSyntheticLightCurve(80, 30, new TransitSignal[0], 0.0003);
        var config = DefaultConfig();
        config.MinTransitDepthPpm = 1000;

        // Tiny depth: 10 ppm
        var parameters = new TransitParameters(5.0, 1.0, 0.00001, 0.1, 0.003, 0.0, 0.015);
        var candidate = new TransitCandidate(
            parameters, 0.0, TransitDisposition.Unknown,
            LightCurve.FromArrays(new double[] { 0 }, new double[] { 1 }),
            new TransitFeatureSet());

        var result = TransitValidator.Validate(candidate, lc, config);

        Assert.IsFalse(result.IsValid, "Shallow transit should be rejected.");
    }

    [TestMethod]
    public void TransitValidator_OddEvenRatio_ShouldBeNearOne()
    {
        double truePeriod = 5.0;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = 0.01, Duration = 0.1 } };
        var lc = GenerateSyntheticLightCurve(80, 30, signals, 0.0003);
        var config = DefaultConfig();

        var parameters = new TransitParameters(truePeriod, 1.0, 0.01, 0.1, 0.1, 0.0, 0.015);
        var candidate = new TransitCandidate(
            parameters, 0.0, TransitDisposition.Unknown,
            LightCurve.FromArrays(new double[] { 0 }, new double[] { 1 }),
            new TransitFeatureSet());

        var result = TransitValidator.Validate(candidate, lc, config);

        Assert.IsTrue(result.OddEvenDepthRatio > 0.5 && result.OddEvenDepthRatio < 2.0,
            $"Odd/even ratio = {result.OddEvenDepthRatio:F2}, expected near 1.0 for planet transit");
    }

    #endregion

    #region TransitDetectionPipeline Tests

    [TestMethod]
    public void Pipeline_SinglePlanet_ShouldDetectTransit()
    {
        double truePeriod = 5.0;
        double trueDepth = 0.015;
        var signals = new[] { new TransitSignal { Period = truePeriod, Epoch = 1.0, Depth = trueDepth, Duration = 0.12 } };
        var lc = GenerateSyntheticLightCurve(100, 30, signals, 0.0003);
        var config = DefaultConfig();
        config.SnrThreshold = 5.0;

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.IsTrue(candidates.Length >= 1,
            $"Should detect at least 1 planet, found {candidates.Length}");

        double periodError = Math.Abs(candidates[0].Parameters.Period - truePeriod) / truePeriod;
        Assert.IsTrue(periodError < 0.02,
            $"Detected period = {candidates[0].Parameters.Period:F4}, expected {truePeriod}, error = {periodError:P1}");
    }

    [TestMethod]
    public void Pipeline_NoPlanet_ShouldReturnEmpty()
    {
        var lc = GenerateSyntheticLightCurve(80, 30, new TransitSignal[0], 0.0003);
        var config = DefaultConfig();

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.AreEqual(0, candidates.Length,
            $"Should detect 0 planets in noise-only data, found {candidates.Length}");
    }

    [TestMethod]
    public void Pipeline_TwoPlanets_ShouldDetectBoth()
    {
        var signals = new[]
        {
            new TransitSignal { Period = 4.0, Epoch = 0.5, Depth = 0.015, Duration = 0.12 },
            new TransitSignal { Period = 7.3, Epoch = 2.0, Depth = 0.008, Duration = 0.10 }
        };
        var lc = GenerateSyntheticLightCurve(120, 30, signals, 0.0003);
        var config = DefaultConfig();
        config.MaxPlanets = 3;
        config.SnrThreshold = 5.0;

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.IsTrue(candidates.Length >= 2,
            $"Should detect at least 2 planets, found {candidates.Length}");

        // Verify both periods are recovered (in any order)
        var detectedPeriods = candidates.Select(c => c.Parameters.Period).ToArray();
        bool foundFirst = detectedPeriods.Any(p => Math.Abs(p - 4.0) / 4.0 < 0.03);
        bool foundSecond = detectedPeriods.Any(p => Math.Abs(p - 7.3) / 7.3 < 0.03);

        Assert.IsTrue(foundFirst,
            $"Should find P≈4.0d. Detected: {string.Join(", ", detectedPeriods.Select(p => $"{p:F3}"))}");
        Assert.IsTrue(foundSecond,
            $"Should find P≈7.3d. Detected: {string.Join(", ", detectedPeriods.Select(p => $"{p:F3}"))}");
    }

    [TestMethod]
    public void Pipeline_CandidateDisposition_ShouldBeSet()
    {
        var signals = new[] { new TransitSignal { Period = 5.0, Epoch = 1.0, Depth = 0.015, Duration = 0.12 } };
        var lc = GenerateSyntheticLightCurve(100, 30, signals, 0.0003);
        var config = DefaultConfig();

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.IsTrue(candidates.Length >= 1);
        Assert.AreNotEqual(TransitDisposition.FalsePositive, candidates[0].Disposition);
        Assert.IsTrue(candidates[0].Score > 0);
    }

    [TestMethod]
    public void Pipeline_PhaseFoldedCurve_ShouldBeProvided()
    {
        var signals = new[] { new TransitSignal { Period = 5.0, Epoch = 1.0, Depth = 0.015, Duration = 0.12 } };
        var lc = GenerateSyntheticLightCurve(100, 30, signals, 0.0003);
        var config = DefaultConfig();

        var candidates = TransitDetectionPipeline.Detect(lc, config);

        Assert.IsTrue(candidates.Length >= 1);
        Assert.IsNotNull(candidates[0].PhaseFoldedCurve);
        Assert.IsTrue(candidates[0].PhaseFoldedCurve.Length > 0,
            "Phase-folded curve should have data points.");
    }

    #endregion

    #region TransitDetectionConfig Tests

    [TestMethod]
    public void TransitDetectionConfig_DefaultValues_ShouldBeReasonable()
    {
        var config = new TransitDetectionConfig();

        Assert.AreEqual(0.5, config.MinPeriodDays);
        Assert.AreEqual(100.0, config.MaxPeriodDays);
        Assert.AreEqual(7.0, config.SnrThreshold);
        Assert.AreEqual(3, config.MaxPlanets);
        Assert.AreEqual(DetrendingMethod.MedianFilter, config.DetrendingMethod);
        Assert.AreEqual(PeriodSearchMethod.BLS, config.PeriodSearchMethod);
    }

    #endregion
}
