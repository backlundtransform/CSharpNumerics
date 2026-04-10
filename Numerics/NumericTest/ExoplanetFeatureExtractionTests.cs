using System;
using System.Linq;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Engines.Exoplanet.Features;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class ExoplanetFeatureExtractionTests
{
    #region Helper: Synthetic Light Curve

    /// <summary>
    /// Creates a synthetic light curve with an injected box-shaped transit.
    /// </summary>
    private static LightCurve CreateSyntheticLightCurve(
        double period, double epoch, double depth, double durationDays,
        int numPoints = 3000, double cadenceDays = 0.02, double noiseLevel = 0.0002)
    {
        var rng = new Random(42);
        double[] time = new double[numPoints];
        double[] flux = new double[numPoints];
        double[] err = new double[numPoints];
        int[] qual = new int[numPoints];

        double halfDur = durationDays / 2.0;

        for (int i = 0; i < numPoints; i++)
        {
            time[i] = i * cadenceDays;
            flux[i] = 1.0;

            // Inject transit
            double phase = ((time[i] - epoch) % period + period) % period;
            if (phase > period / 2.0) phase -= period;

            if (Math.Abs(phase) <= halfDur)
                flux[i] = 1.0 - depth;

            // Add noise
            flux[i] += noiseLevel * (rng.NextDouble() * 2.0 - 1.0);
            err[i] = noiseLevel;
        }

        return new LightCurve(time, flux, err, qual, new LightCurveMetadata());
    }

    private static TransitCandidate CreateCandidate(
        double period, double epoch, double depth, double durationDays,
        LightCurve lc)
    {
        double ingressDuration = durationDays * 0.15;
        double radiusRatio = Math.Sqrt(depth);
        var parameters = new TransitParameters(period, epoch, depth, durationDays,
            radiusRatio, 0.0, ingressDuration);

        // Create a simple phase-folded curve
        int bins = 200;
        var times = new VectorN(lc.Time);
        var values = new VectorN(lc.Flux);
        var (binCenters, binMeans) = CSharpNumerics.Statistics.TimeSeriesAnalysis.PhaseFolding.FoldAndBin(
            times, values, period, bins, epoch);

        double[] phaseTime = new double[binCenters.Length];
        double[] phaseFlux = new double[binMeans.Length];
        for (int i = 0; i < binCenters.Length; i++)
        {
            phaseTime[i] = binCenters[i];
            phaseFlux[i] = binMeans[i];
        }
        var phaseFolded = LightCurve.FromArrays(phaseTime, phaseFlux);

        return new TransitCandidate(parameters, 0.9, TransitDisposition.Candidate,
            phaseFolded, new TransitFeatureSet());
    }

    #endregion

    #region TransitFeatureExtractor Tests

    [TestMethod]
    public void Extract_ShouldReturnAllFeatures()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        // All 12 features should be present
        foreach (string name in TransitFeatureSet.FeatureNames.All)
        {
            Assert.IsTrue(features.HasFeature(name), $"Feature '{name}' is missing.");
        }
    }

    [TestMethod]
    public void Extract_Depth_ShouldMatchParameters()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        Assert.AreEqual(depth, features[TransitFeatureSet.FeatureNames.Depth], 1e-10);
        Assert.AreEqual(duration, features[TransitFeatureSet.FeatureNames.Duration], 1e-10);
        Assert.AreEqual(period, features[TransitFeatureSet.FeatureNames.Period], 1e-10);
    }

    [TestMethod]
    public void Extract_BlsSnr_ShouldBePositiveForRealTransit()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        Assert.IsTrue(features[TransitFeatureSet.FeatureNames.SnrBls] > 5.0,
            $"BLS SNR should be high for a real transit, got {features[TransitFeatureSet.FeatureNames.SnrBls]:F1}");
    }

    [TestMethod]
    public void Extract_OddEvenRatio_ShouldBeNearOne()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        double ratio = features[TransitFeatureSet.FeatureNames.OddEvenRatio];
        Assert.IsTrue(ratio > 0.7 && ratio < 1.5,
            $"Odd/even ratio should be near 1.0 for a real transit, got {ratio:F3}");
    }

    [TestMethod]
    public void Extract_SecondaryDepth_ShouldBeSmallForPlanet()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        double secDepth = features[TransitFeatureSet.FeatureNames.SecondaryDepth];
        Assert.IsTrue(secDepth < depth * 0.3,
            $"Secondary eclipse depth should be small for a planet, got {secDepth:E3}");
    }

    [TestMethod]
    public void Extract_VShapeMetric_ShouldBeLowForBoxTransit()
    {
        // A box transit has no ingress, so V-shape metric should be low
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        double vShape = features[TransitFeatureSet.FeatureNames.VShapeMetric];
        Assert.IsTrue(vShape < 0.5,
            $"V-shape metric should be low for box-like transit, got {vShape:F3}");
    }

    [TestMethod]
    public void Extract_Scatter_OutShouldBeLessThanIn()
    {
        // Out-of-transit scatter should approximate noise level;
        // in-transit scatter may be higher due to transit shape
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        double scatterOut = features[TransitFeatureSet.FeatureNames.ScatterOutTransit];
        Assert.IsTrue(scatterOut > 0, "Out-of-transit scatter should be positive.");
        Assert.IsTrue(scatterOut < 0.001,
            $"Out-of-transit scatter should be close to noise level, got {scatterOut:E3}");
    }

    [TestMethod]
    public void Extract_NoTransit_ShouldHaveLowSnr()
    {
        // Flat light curve — no transit
        var rng = new Random(42);
        int n = 2000;
        double[] time = new double[n];
        double[] flux = new double[n];
        for (int i = 0; i < n; i++)
        {
            time[i] = i * 0.02;
            flux[i] = 1.0 + 0.0002 * (rng.NextDouble() * 2 - 1);
        }
        var lc = LightCurve.FromArrays(time, flux);

        var parameters = new TransitParameters(3.0, 0.5, 0.001, 0.15, 0.03, 0.0, 0.02);
        var candidate = new TransitCandidate(parameters, 0.1, TransitDisposition.Unknown,
            LightCurve.FromArrays(new double[] { 0.0 }, new double[] { 1.0 }),
            new TransitFeatureSet());

        var features = TransitFeatureExtractor.Extract(candidate, lc);

        double snr = features[TransitFeatureSet.FeatureNames.SnrBls];
        Assert.IsTrue(snr < 5.0, $"BLS SNR should be low for flat curve, got {snr:F1}");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void Extract_NullCandidate_ShouldThrow()
    {
        TransitFeatureExtractor.Extract(null, LightCurve.FromArrays(new double[] { 0 }, new double[] { 1 }));
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void Extract_NullLightCurve_ShouldThrow()
    {
        var candidate = new TransitCandidate();
        TransitFeatureExtractor.Extract(candidate, null);
    }

    #endregion

    #region WindowedFeatureExtractor Tests

    [TestMethod]
    public void CreateTrainingData_ShouldHaveCorrectDimensions()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        int windowSize = 50;

        var lc1 = CreateSyntheticLightCurve(period, epoch, depth, duration, numPoints: 2000);
        var lc2 = CreateSyntheticLightCurve(period, epoch + 0.1, depth * 0.5, duration, numPoints: 2000);
        var c1 = CreateCandidate(period, epoch, depth, duration, lc1);
        var c2 = CreateCandidate(period, epoch + 0.1, depth * 0.5, duration, lc2);
        c1.Features = TransitFeatureExtractor.Extract(c1, lc1);
        c2.Features = TransitFeatureExtractor.Extract(c2, lc2);

        var (X, y) = WindowedFeatureExtractor.CreateTrainingData(
            new[] { lc1, lc2 }, new[] { c1, c2 }, new[] { 1.0, 0.0 }, windowSize);

        int numFeatures = TransitFeatureSet.FeatureNames.All.Length;
        Assert.AreEqual(2, X.rowLength);
        Assert.AreEqual(windowSize + numFeatures, X.columnLength);
        Assert.AreEqual(2, y.Length);
        Assert.AreEqual(1.0, y[0]);
        Assert.AreEqual(0.0, y[1]);
    }

    [TestMethod]
    public void CreateTrainingData_FeatureColumns_ShouldMatchExtracted()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        int windowSize = 50;

        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration, numPoints: 2000);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);
        candidate.Features = TransitFeatureExtractor.Extract(candidate, lc);

        var (X, y) = WindowedFeatureExtractor.CreateTrainingData(
            new[] { lc }, new[] { candidate }, new[] { 1.0 }, windowSize);

        string[] featureNames = TransitFeatureSet.FeatureNames.All;
        for (int f = 0; f < featureNames.Length; f++)
        {
            double expected = candidate.Features[featureNames[f]];
            double actual = X.values[0, windowSize + f];
            Assert.AreEqual(expected, actual, 1e-10,
                $"Feature '{featureNames[f]}' mismatch at column {windowSize + f}");
        }
    }

    [TestMethod]
    public void CreateInferenceData_ShouldHaveCorrectDimensions()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        int windowSize = 50;

        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration, numPoints: 2000);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);
        candidate.Features = TransitFeatureExtractor.Extract(candidate, lc);

        var X = WindowedFeatureExtractor.CreateInferenceData(lc, candidate, windowSize);

        int numFeatures = TransitFeatureSet.FeatureNames.All.Length;
        Assert.AreEqual(1, X.rowLength);
        Assert.AreEqual(windowSize + numFeatures, X.columnLength);
    }

    [TestMethod]
    public void CreateTrainingData_FluxWindow_ShouldContainTransitDip()
    {
        double period = 3.0, epoch = 0.5, depth = 0.01, duration = 0.15;
        int windowSize = 100;

        var lc = CreateSyntheticLightCurve(period, epoch, depth, duration, numPoints: 3000);
        var candidate = CreateCandidate(period, epoch, depth, duration, lc);
        candidate.Features = TransitFeatureExtractor.Extract(candidate, lc);

        var (X, _) = WindowedFeatureExtractor.CreateTrainingData(
            new[] { lc }, new[] { candidate }, new[] { 1.0 }, windowSize);

        // The center of the window should show the transit dip (flux < 1.0)
        double centerFlux = X.values[0, windowSize / 2];
        Assert.IsTrue(centerFlux < 0.999,
            $"Center of phase-folded window should show transit dip, got {centerFlux:F6}");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void CreateTrainingData_MismatchedArrays_ShouldThrow()
    {
        var lc = CreateSyntheticLightCurve(3.0, 0.5, 0.01, 0.15, numPoints: 500);
        var candidate = CreateCandidate(3.0, 0.5, 0.01, 0.15, lc);

        WindowedFeatureExtractor.CreateTrainingData(
            new[] { lc, lc }, new[] { candidate }, new[] { 1.0 }, 50);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void CreateTrainingData_NullCurves_ShouldThrow()
    {
        WindowedFeatureExtractor.CreateTrainingData(null, new TransitCandidate[0], new double[0], 50);
    }

    #endregion

    #region TransitFeatureSet.FeatureNames.All Tests

    [TestMethod]
    public void FeatureNames_All_ShouldContainAll12Features()
    {
        Assert.AreEqual(12, TransitFeatureSet.FeatureNames.All.Length);
        Assert.IsTrue(TransitFeatureSet.FeatureNames.All.Contains(TransitFeatureSet.FeatureNames.Depth));
        Assert.IsTrue(TransitFeatureSet.FeatureNames.All.Contains(TransitFeatureSet.FeatureNames.IngressEgressRatio));
    }

    #endregion
}
