using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Statistics.Data;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class ExoplanetDataModelTests
{
    #region LightCurve Tests

    [TestMethod]
    public void LightCurve_Constructor_ShouldStoreAllFields()
    {
        var time = new double[] { 1.0, 2.0, 3.0 };
        var flux = new double[] { 1.0, 0.99, 1.0 };
        var err = new double[] { 0.001, 0.001, 0.001 };
        var qual = new int[] { 0, 0, 0 };
        var meta = new LightCurveMetadata("TIC-12345", "TESS", 2457000.0, CadenceType.Short);

        var lc = new LightCurve(time, flux, err, qual, meta);

        Assert.AreEqual(3, lc.Length);
        Assert.AreEqual(0.99, lc.Flux[1]);
        Assert.AreEqual("TIC-12345", lc.Metadata.TargetId);
        Assert.AreEqual("TESS", lc.Metadata.Mission);
        Assert.AreEqual(CadenceType.Short, lc.Metadata.Cadence);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void LightCurve_MismatchedArrayLengths_ShouldThrow()
    {
        var time = new double[] { 1.0, 2.0 };
        var flux = new double[] { 1.0 };
        var err = new double[] { 0.001, 0.001 };
        var qual = new int[] { 0, 0 };

        new LightCurve(time, flux, err, qual, new LightCurveMetadata());
    }

    [TestMethod]
    public void LightCurve_FromArrays_WithFluxOnly_ShouldCreateDefaultErrors()
    {
        var time = new double[] { 1.0, 2.0, 3.0 };
        var flux = new double[] { 1.0, 0.99, 1.0 };

        var lc = LightCurve.FromArrays(time, flux);

        Assert.AreEqual(3, lc.Length);
        Assert.AreEqual(0.0, lc.FluxError[0]);
        Assert.AreEqual(0, lc.QualityFlags[0]);
    }

    [TestMethod]
    public void LightCurve_FromTimeSeries_ShouldExtractColumns()
    {
        var times = new DateTime[]
        {
            new DateTime(2020, 1, 1),
            new DateTime(2020, 1, 2),
            new DateTime(2020, 1, 3)
        };
        var bjd = new double[] { 2458849.5, 2458850.5, 2458851.5 };
        var flux = new double[] { 1.0, 0.99, 1.0 };
        var err = new double[] { 0.001, 0.002, 0.001 };

        var ts = new TimeSeries(times, new[] { bjd, flux, err }, new[] { "TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR" });

        var lc = LightCurve.FromTimeSeries(ts, "TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR");

        Assert.AreEqual(3, lc.Length);
        Assert.AreEqual(2458849.5, lc.Time[0]);
        Assert.AreEqual(0.99, lc.Flux[1]);
        Assert.AreEqual(0.002, lc.FluxError[1]);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void LightCurve_FromTimeSeries_MissingColumn_ShouldThrow()
    {
        var times = new DateTime[] { new DateTime(2020, 1, 1) };
        var ts = new TimeSeries(times, new[] { new double[] { 1.0 } }, new[] { "Flux" });

        LightCurve.FromTimeSeries(ts, "TIME", "Flux", "Error");
    }

    #endregion

    #region TransitParameters Tests

    [TestMethod]
    public void TransitParameters_ShouldStoreValues()
    {
        var p = new TransitParameters(3.5, 100.0, 0.01, 2.5, 0.1, 0.3, 0.25);

        Assert.AreEqual(3.5, p.Period);
        Assert.AreEqual(100.0, p.Epoch);
        Assert.AreEqual(0.01, p.Depth);
        Assert.AreEqual(2.5, p.Duration);
        Assert.AreEqual(0.1, p.RadiusRatio);
        Assert.AreEqual(0.3, p.ImpactParameter);
        Assert.AreEqual(0.25, p.IngressDuration);
    }

    #endregion

    #region StellarProperties Tests

    [TestMethod]
    public void StellarProperties_ShouldStoreValues()
    {
        var star = new StellarProperties(5778.0, 1.0, 1.0, 4.44, 0.0, SpectralType.G);

        Assert.AreEqual(5778.0, star.EffectiveTemp);
        Assert.AreEqual(1.0, star.Radius);
        Assert.AreEqual(1.0, star.Mass);
        Assert.AreEqual(4.44, star.SurfaceGravity);
        Assert.AreEqual(0.0, star.Metallicity);
        Assert.AreEqual(SpectralType.G, star.Type);
    }

    #endregion

    #region TransitFeatureSet Tests

    [TestMethod]
    public void TransitFeatureSet_ShouldStoreAndRetrieveFeatures()
    {
        var fs = new TransitFeatureSet();
        fs[TransitFeatureSet.FeatureNames.Depth] = 0.01;
        fs[TransitFeatureSet.FeatureNames.Period] = 3.5;
        fs[TransitFeatureSet.FeatureNames.SnrBls] = 12.5;

        Assert.AreEqual(0.01, fs[TransitFeatureSet.FeatureNames.Depth]);
        Assert.AreEqual(3.5, fs[TransitFeatureSet.FeatureNames.Period]);
        Assert.IsTrue(fs.HasFeature(TransitFeatureSet.FeatureNames.Depth));
        Assert.IsFalse(fs.HasFeature(TransitFeatureSet.FeatureNames.VShapeMetric));
    }

    #endregion

    #region TransitCandidate Tests

    [TestMethod]
    public void TransitCandidate_ShouldStoreAllFields()
    {
        var parameters = new TransitParameters(3.5, 100.0, 0.01, 2.5, 0.1, 0.3, 0.25);
        var features = new TransitFeatureSet();
        features[TransitFeatureSet.FeatureNames.Depth] = 0.01;
        var phaseFolded = LightCurve.FromArrays(new double[] { -0.5, 0.0, 0.5 }, new double[] { 1.0, 0.99, 1.0 });

        var candidate = new TransitCandidate(parameters, 0.95, TransitDisposition.Candidate, phaseFolded, features);

        Assert.AreEqual(0.95, candidate.Score);
        Assert.AreEqual(TransitDisposition.Candidate, candidate.Disposition);
        Assert.AreEqual(3.5, candidate.Parameters.Period);
        Assert.AreEqual(3, candidate.PhaseFoldedCurve.Length);
    }

    #endregion

    #region LightCurveMetadata Tests

    [TestMethod]
    public void LightCurveMetadata_DefaultConstructor_ShouldHaveDefaults()
    {
        var meta = new LightCurveMetadata();

        Assert.AreEqual(string.Empty, meta.TargetId);
        Assert.AreEqual(string.Empty, meta.Mission);
        Assert.AreEqual(0.0, meta.TimeOffset);
    }

    #endregion

    #region LightCurveSanitizer Tests

    [TestMethod]
    public void RemoveBadQuality_ShouldFilterFlaggedPoints()
    {
        var lc = new LightCurve(
            new double[] { 1, 2, 3, 4, 5 },
            new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 },
            new double[] { 0.001, 0.001, 0.001, 0.001, 0.001 },
            new int[] { 0, 0, 128, 0, 256 },
            new LightCurveMetadata());

        var cleaned = LightCurveSanitizer.RemoveBadQuality(lc);

        Assert.AreEqual(3, cleaned.Length);
        Assert.AreEqual(1.0, cleaned.Time[0]);
        Assert.AreEqual(2.0, cleaned.Time[1]);
        Assert.AreEqual(4.0, cleaned.Time[2]);
    }

    [TestMethod]
    public void RemoveOutliers_ShouldRemoveSigmaOutliers()
    {
        int n = 100;
        var time = new double[n];
        var flux = new double[n];
        var err = new double[n];
        var qual = new int[n];

        for (int i = 0; i < n; i++)
        {
            time[i] = i;
            flux[i] = 1.0;
            err[i] = 0.001;
        }

        flux[50] = 100.0;
        flux[75] = -50.0;

        var lc = new LightCurve(time, flux, err, qual, new LightCurveMetadata());
        var cleaned = LightCurveSanitizer.RemoveOutliers(lc, 3.0);

        Assert.IsTrue(cleaned.Length < n);
        Assert.IsTrue(cleaned.Flux.All(f => Math.Abs(f) < 10));
    }

    [TestMethod]
    public void NormalizeFlux_ShouldDivideByMedian()
    {
        var lc = new LightCurve(
            new double[] { 1, 2, 3, 4, 5 },
            new double[] { 1000, 1001, 999, 1000, 998 },
            new double[] { 1, 1, 1, 1, 1 },
            new int[] { 0, 0, 0, 0, 0 },
            new LightCurveMetadata());

        var normalized = LightCurveSanitizer.NormalizeFlux(lc);

        Assert.AreEqual(5, normalized.Length);
        Assert.AreEqual(1.0, normalized.Flux[0], 0.01);
        for (int i = 0; i < normalized.Length; i++)
        {
            Assert.IsTrue(Math.Abs(normalized.Flux[i] - 1.0) < 0.01);
        }
    }

    [TestMethod]
    public void FillGaps_ShouldInterpolateSmallGaps()
    {
        var lc = new LightCurve(
            new double[] { 1.0, 2.0, 3.0, 6.0, 7.0 },
            new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 },
            new double[] { 0.001, 0.001, 0.001, 0.001, 0.001 },
            new int[] { 0, 0, 0, 0, 0 },
            new LightCurveMetadata());

        var filled = LightCurveSanitizer.FillGaps(lc, 5.0);

        Assert.IsTrue(filled.Length > 5);
        for (int i = 0; i < filled.Length - 1; i++)
        {
            double dt = filled.Time[i + 1] - filled.Time[i];
            Assert.IsTrue(dt < 2.0, $"Gap at index {i}: dt={dt}");
        }
    }

    [TestMethod]
    public void FillGaps_ShouldNotFillLargeGaps()
    {
        var lc = new LightCurve(
            new double[] { 1.0, 2.0, 100.0, 101.0 },
            new double[] { 1.0, 1.0, 1.0, 1.0 },
            new double[] { 0.001, 0.001, 0.001, 0.001 },
            new int[] { 0, 0, 0, 0 },
            new LightCurveMetadata());

        var filled = LightCurveSanitizer.FillGaps(lc, 5.0);

        Assert.AreEqual(4, filled.Length);
    }

    [TestMethod]
    public void Sanitizer_FullPipeline_ShouldChainOperations()
    {
        int n = 200;
        var time = new double[n];
        var flux = new double[n];
        var err = new double[n];
        var qual = new int[n];

        for (int i = 0; i < n; i++)
        {
            time[i] = i * 0.02;
            flux[i] = 1000.0;
            err[i] = 1.0;
        }

        flux[10] = 5000.0;
        qual[20] = 128;
        qual[21] = 256;

        var lc = new LightCurve(time, flux, err, qual, new LightCurveMetadata("TIC-999", "TESS", 2457000, CadenceType.Short));

        var result = LightCurveSanitizer.RemoveBadQuality(lc);
        result = LightCurveSanitizer.RemoveOutliers(result, 3.0);
        result = LightCurveSanitizer.NormalizeFlux(result);

        Assert.IsTrue(result.Length < n);
        Assert.AreEqual("TIC-999", result.Metadata.TargetId);
        Assert.IsTrue(result.Flux.All(f => Math.Abs(f - 1.0) < 0.01));
    }

    #endregion

    #region Enum Tests

    [TestMethod]
    public void TransitDisposition_ShouldHaveExpectedValues()
    {
        Assert.AreEqual(0, (int)TransitDisposition.Unknown);
        Assert.AreEqual(1, (int)TransitDisposition.Candidate);
        Assert.AreEqual(2, (int)TransitDisposition.Confirmed);
        Assert.AreEqual(3, (int)TransitDisposition.FalsePositive);
        Assert.AreEqual(4, (int)TransitDisposition.AstrophysicalFalsePositive);
    }

    [TestMethod]
    public void CadenceType_ShouldHaveExpectedValues()
    {
        Assert.AreEqual(0, (int)CadenceType.Short);
        Assert.AreEqual(1, (int)CadenceType.Long);
        Assert.AreEqual(2, (int)CadenceType.Fast);
    }

    #endregion
}
