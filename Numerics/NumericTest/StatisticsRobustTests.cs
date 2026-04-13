using System;
using System.Linq;
using CSharpNumerics.Statistics.Robust;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;
using CSharpNumerics.Numerics.Objects;

namespace NumericTest;

[TestClass]
public class StatisticsRobustTests
{
    #region SigmaClipping Tests

    [TestMethod]
    public void SigmaClipping_NoOutliers_ShouldRetainAll()
    {
        var data = new double[] { 1.0, 1.1, 0.9, 1.05, 0.95, 1.02, 0.98 };

        var result = SigmaClipping.Clip(data, 3.0, 3.0);

        Assert.AreEqual(data.Length, result.RetainedCount);
        Assert.AreEqual(0, result.ClippedCount);
        Assert.IsTrue(result.Mask.All(m => m));
    }

    [TestMethod]
    public void SigmaClipping_WithOutliers_ShouldClipThem()
    {
        // 100 normal points around 0, plus 2 extreme outliers
        var rng = new Random(42);
        var data = new double[102];
        for (int i = 0; i < 100; i++)
            data[i] = 0.001 * (rng.NextDouble() * 2 - 1);
        data[100] = 10.0;   // extreme high outlier
        data[101] = -10.0;  // extreme low outlier

        var result = SigmaClipping.Clip(data, 3.0, 3.0, maxIter: 5);

        Assert.IsTrue(result.ClippedCount >= 2, $"Expected at least 2 clipped, got {result.ClippedCount}");
        Assert.IsFalse(result.Mask[100], "High outlier should be clipped.");
        Assert.IsFalse(result.Mask[101], "Low outlier should be clipped.");
    }

    [TestMethod]
    public void SigmaClipping_AsymmetricThresholds_ShouldClipOneSide()
    {
        // Data centered around 0, with positive outlier
        var data = new double[] { 0.0, 0.1, -0.1, 0.05, -0.05, 5.0 };

        // Very strict upper clip (1σ), generous lower clip (100σ)
        var result = SigmaClipping.Clip(data, sigmaLow: 100.0, sigmaHigh: 1.0, maxIter: 10);

        Assert.IsFalse(result.Mask[5], "Positive outlier should be clipped with strict upper threshold.");
    }

    [TestMethod]
    public void SigmaClipping_Convergence_ShouldStopWhenNoMoreClipping()
    {
        var data = new double[] { 1.0, 1.01, 0.99, 1.02, 0.98, 1.005, 0.995 };

        var result = SigmaClipping.Clip(data, 5.0, 5.0, maxIter: 100);

        // No outliers at 5σ — should converge in 1 iteration
        Assert.AreEqual(1, result.Iterations);
        Assert.AreEqual(data.Length, result.RetainedCount);
    }

    [TestMethod]
    public void SigmaClipping_MeanAndStd_ShouldBeCorrectAfterClipping()
    {
        var data = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0, 100.0 };

        var result = SigmaClipping.Clip(data, 2.0, 2.0, maxIter: 10);

        // After clipping the outlier, mean should be close to 1.0
        Assert.AreEqual(1.0, result.Mean, 0.001);
        Assert.IsTrue(result.Std < 0.01, $"Std should be near 0, got {result.Std}");
    }

    [TestMethod]
    public void SigmaClipping_SymmetricOverload_ShouldWork()
    {
        // 20 points near 0, plus 2 extreme outliers that will be clipped
        var data = new double[22];
        var rng = new Random(42);
        for (int i = 0; i < 20; i++)
            data[i] = 0.01 * (rng.NextDouble() * 2 - 1);
        data[20] = 50.0;
        data[21] = -50.0;

        var result = SigmaClipping.Clip(data, sigma: 3.0, maxIter: 10);

        Assert.IsTrue(result.ClippedCount >= 2);
    }

    [TestMethod]
    public void SigmaClipping_Apply_ShouldReturnOnlyRetainedValues()
    {
        // 20 points at 1.0, one extreme outlier
        var data = new double[21];
        for (int i = 0; i < 20; i++)
            data[i] = 1.0;
        data[20] = 100.0;

        var retained = SigmaClipping.Apply(data, 2.0, 2.0);

        Assert.AreEqual(20, retained.Length);
        Assert.IsTrue(retained.All(v => Math.Abs(v - 1.0) < 0.01));
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void SigmaClipping_NullData_ShouldThrow()
    {
        SigmaClipping.Clip(null);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void SigmaClipping_EmptyData_ShouldThrow()
    {
        SigmaClipping.Clip(new double[0]);
    }

    #endregion

    #region SlidingWindowStatistics Tests

    [TestMethod]
    public void RunningMean_ConstantData_ShouldReturnSameValue()
    {
        var data = new double[] { 5.0, 5.0, 5.0, 5.0, 5.0 };

        var result = SlidingWindowStatistics.RunningMean(data, 3);

        for (int i = 0; i < data.Length; i++)
            Assert.AreEqual(5.0, result[i], 1e-10);
    }

    [TestMethod]
    public void RunningMean_KnownValues_ShouldBeCorrect()
    {
        var data = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };

        var result = SlidingWindowStatistics.RunningMean(data, 3);

        // Window size 3 (already odd): centered windows
        // i=0: [1,2] → 1.5 (edge)
        // i=1: [1,2,3] → 2.0
        // i=2: [2,3,4] → 3.0
        // i=3: [3,4,5] → 4.0
        // i=4: [4,5] → 4.5 (edge)
        Assert.AreEqual(1.5, result[0], 1e-10);
        Assert.AreEqual(2.0, result[1], 1e-10);
        Assert.AreEqual(3.0, result[2], 1e-10);
        Assert.AreEqual(4.0, result[3], 1e-10);
        Assert.AreEqual(4.5, result[4], 1e-10);
    }

    [TestMethod]
    public void RunningStd_ConstantData_ShouldReturnZero()
    {
        var data = new double[] { 3.0, 3.0, 3.0, 3.0, 3.0 };

        var result = SlidingWindowStatistics.RunningStd(data, 3);

        for (int i = 0; i < data.Length; i++)
            Assert.AreEqual(0.0, result[i], 1e-10);
    }

    [TestMethod]
    public void RunningStd_VariedData_ShouldBePositive()
    {
        var data = new double[] { 1.0, 3.0, 1.0, 3.0, 1.0, 3.0 };

        var result = SlidingWindowStatistics.RunningStd(data, 3);

        for (int i = 1; i < data.Length - 1; i++)
            Assert.IsTrue(result[i] > 0, $"Std at index {i} should be positive.");
    }

    [TestMethod]
    public void RunningMAD_ConstantData_ShouldReturnZero()
    {
        var data = new double[] { 7.0, 7.0, 7.0, 7.0, 7.0 };

        var result = SlidingWindowStatistics.RunningMAD(data, 3);

        for (int i = 0; i < data.Length; i++)
            Assert.AreEqual(0.0, result[i], 1e-10);
    }

    [TestMethod]
    public void RunningMAD_WithOutlier_ShouldBeRobust()
    {
        // Most data is 1.0, one outlier at 100.0
        var data = new double[] { 1.0, 1.0, 1.0, 100.0, 1.0, 1.0, 1.0 };

        var mad = SlidingWindowStatistics.RunningMAD(data, 5);
        var std = SlidingWindowStatistics.RunningStd(data, 5);

        // MAD should be much smaller than std at the outlier position
        // because MAD is robust to single outliers
        Assert.IsTrue(mad[0] < std[3],
            $"MAD ({mad[0]}) should be much smaller than std at outlier ({std[3]})");
    }

    [TestMethod]
    public void RunningMedian_KnownValues()
    {
        var data = new double[] { 1.0, 5.0, 3.0, 2.0, 4.0 };

        var result = SlidingWindowStatistics.RunningMedian(data, 3);

        // i=1: [1,5,3] sorted=[1,3,5] → median=3
        // i=2: [5,3,2] sorted=[2,3,5] → median=3
        // i=3: [3,2,4] sorted=[2,3,4] → median=3
        Assert.AreEqual(3.0, result[1], 1e-10);
        Assert.AreEqual(3.0, result[2], 1e-10);
        Assert.AreEqual(3.0, result[3], 1e-10);
    }

    [TestMethod]
    public void SlidingWindow_OutputLength_ShouldMatchInput()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        Assert.AreEqual(data.Length, SlidingWindowStatistics.RunningMean(data, 5).Length);
        Assert.AreEqual(data.Length, SlidingWindowStatistics.RunningStd(data, 5).Length);
        Assert.AreEqual(data.Length, SlidingWindowStatistics.RunningMAD(data, 5).Length);
        Assert.AreEqual(data.Length, SlidingWindowStatistics.RunningMedian(data, 5).Length);
    }

    [TestMethod]
    public void SlidingWindow_EvenWindowSize_ShouldIncrementToOdd()
    {
        // Even window size 4 → treated as 5
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7 };
        var mean4 = SlidingWindowStatistics.RunningMean(data, 4);
        var mean5 = SlidingWindowStatistics.RunningMean(data, 5);

        // Should be identical since 4 is incremented to 5
        for (int i = 0; i < data.Length; i++)
            Assert.AreEqual(mean5[i], mean4[i], 1e-10);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void RunningMean_NullData_ShouldThrow()
    {
        SlidingWindowStatistics.RunningMean(null, 3);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void RunningMean_ZeroWindowSize_ShouldThrow()
    {
        SlidingWindowStatistics.RunningMean(new double[] { 1 }, 0);
    }

    #endregion

    #region FalseAlarmProbability Tests

    [TestMethod]
    public void AnalyticalFAP_HighPower_ShouldBeLow()
    {
        // Very strong signal → FAP should be near 0
        double fap = FalseAlarmProbability.AnalyticalFAP(peakPower: 20.0, numFrequencies: 1000, numDataPoints: 100);
        Assert.IsTrue(fap < 0.01, $"FAP should be very low for strong signal, got {fap}");
    }

    [TestMethod]
    public void AnalyticalFAP_LowPower_ShouldBeHigh()
    {
        // Weak signal → FAP should be near 1
        double fap = FalseAlarmProbability.AnalyticalFAP(peakPower: 0.5, numFrequencies: 1000, numDataPoints: 100);
        Assert.IsTrue(fap > 0.5, $"FAP should be high for weak signal, got {fap}");
    }

    [TestMethod]
    public void AnalyticalFAP_MoreFrequencies_ShouldIncraseFAP()
    {
        double fap100 = FalseAlarmProbability.AnalyticalFAP(5.0, 100, 100);
        double fap10000 = FalseAlarmProbability.AnalyticalFAP(5.0, 10000, 100);

        Assert.IsTrue(fap10000 > fap100,
            $"More frequencies should increase FAP: {fap100} vs {fap10000}");
    }

    [TestMethod]
    public void AnalyticalFAP_ShouldBeInZeroOneRange()
    {
        for (double power = 0; power <= 30; power += 2)
        {
            double fap = FalseAlarmProbability.AnalyticalFAP(power, 500, 200);
            Assert.IsTrue(fap >= 0 && fap <= 1.0,
                $"FAP should be in [0,1], got {fap} for power={power}");
        }
    }

    [TestMethod]
    public void AnalyticalFAP_NegativePower_ShouldReturnOne()
    {
        double fap = FalseAlarmProbability.AnalyticalFAP(-1.0, 100, 100);
        Assert.AreEqual(1.0, fap);
    }

    [TestMethod]
    public void BootstrapFAP_StrongSignal_ShouldBeLow()
    {
        // Create a time series with a clear periodic signal
        int n = 200;
        double period = 5.0;
        double[] times = new double[n];
        double[] values = new double[n];
        var rng = new Random(42);

        for (int i = 0; i < n; i++)
        {
            times[i] = i * 0.1;
            values[i] = Math.Sin(2.0 * Math.PI * times[i] / period)
                + 0.1 * (rng.NextDouble() * 2 - 1);
        }

        // Compute periodogram to get peak power
        var timesVec = new VectorN(times);
        var valuesVec = new VectorN(values);
        var pgram = LombScarglePeriodogram.ComputeByPeriod(timesVec, valuesVec, 2.0, 10.0, 200);

        double fap = FalseAlarmProbability.BootstrapFAP(
            times, values, pgram.BestPower,
            nBootstrap: 100, minPeriod: 2.0, maxPeriod: 10.0,
            numFrequencies: 200, seed: 42);

        Assert.IsTrue(fap < 0.1,
            $"Bootstrap FAP should be low for strong sinusoidal signal, got {fap}");
    }

    [TestMethod]
    public void BootstrapFAP_PureNoise_ShouldBeHigh()
    {
        int n = 100;
        double[] times = new double[n];
        double[] values = new double[n];
        var rng = new Random(42);

        for (int i = 0; i < n; i++)
        {
            times[i] = i * 0.1;
            values[i] = rng.NextDouble() * 2 - 1;
        }

        var timesVec = new VectorN(times);
        var valuesVec = new VectorN(values);
        var pgram = LombScarglePeriodogram.ComputeByPeriod(timesVec, valuesVec, 1.0, 5.0, 100);

        double fap = FalseAlarmProbability.BootstrapFAP(
            times, values, pgram.BestPower,
            nBootstrap: 100, minPeriod: 1.0, maxPeriod: 5.0,
            numFrequencies: 100, seed: 42);

        Assert.IsTrue(fap > 0.05,
            $"Bootstrap FAP should be high for pure noise, got {fap}");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentNullException))]
    public void BootstrapFAP_NullTimes_ShouldThrow()
    {
        FalseAlarmProbability.BootstrapFAP(null, new double[] { 1 }, 1.0);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void BootstrapFAP_MismatchedLengths_ShouldThrow()
    {
        FalseAlarmProbability.BootstrapFAP(
            new double[] { 1, 2, 3, 4 }, new double[] { 1 }, 1.0);
    }

    #endregion

    #region MedianAbsoluteDeviation Tests

    [TestMethod]
    public void MAD_SymmetricData_ShouldReturnCorrectValues()
    {
        var data = new double[] { 1, 2, 3, 4, 5 };
        var result = MedianAbsoluteDeviation.Compute(data);
        Assert.AreEqual(3.0, result.Median, 1e-10);
        Assert.AreEqual(1.0, result.Mad, 1e-10);
        Assert.AreEqual(1.4826, result.ScaledMad, 1e-4);
    }

    [TestMethod]
    public void MAD_SingleElement_ShouldReturnZeroMad()
    {
        var data = new double[] { 42.0 };
        var result = MedianAbsoluteDeviation.Compute(data);
        Assert.AreEqual(42.0, result.Median, 1e-10);
        Assert.AreEqual(0.0, result.Mad, 1e-10);
    }

    [TestMethod]
    public void MAD_WithOutliers_ShouldBeRobust()
    {
        // 98 values near 0, 2 extreme outliers
        var data = new double[100];
        for (int i = 0; i < 98; i++) data[i] = i * 0.01;
        data[98] = 1000.0;
        data[99] = -1000.0;

        var result = MedianAbsoluteDeviation.Compute(data);
        // MAD should not be strongly influenced by the 2 outliers
        Assert.IsTrue(result.Mad < 1.0, $"MAD should be small, got {result.Mad}");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void MAD_EmptyArray_ShouldThrow()
    {
        MedianAbsoluteDeviation.Compute(Array.Empty<double>());
    }

    #endregion

    #region TrimmedMean Tests

    [TestMethod]
    public void TrimmedMean_NoTrim_ShouldReturnArithmeticMean()
    {
        var data = new double[] { 1, 2, 3, 4, 5 };
        double result = TrimmedMean.Compute(data, trimRatio: 0.0);
        Assert.AreEqual(3.0, result, 1e-10);
    }

    [TestMethod]
    public void TrimmedMean_WithOutliers_ShouldResist()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 100 };
        double ordinaryMean = 14.5;
        double trimmed = TrimmedMean.Compute(data, trimRatio: 0.1);
        // Trimmed mean should be much closer to 5 than ordinary mean
        Assert.IsTrue(Math.Abs(trimmed - 5.0) < Math.Abs(ordinaryMean - 5.0),
            $"Trimmed mean {trimmed} should be closer to 5 than ordinary mean {ordinaryMean}");
    }

    [TestMethod]
    public void TrimmedMean_HighTrim_ApproachesMedian()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
        double trimmed = TrimmedMean.Compute(data, trimRatio: 0.45);
        double median = 10.5;
        Assert.AreEqual(median, trimmed, 0.5);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void TrimmedMean_InvalidRatio_ShouldThrow()
    {
        TrimmedMean.Compute(new double[] { 1, 2, 3 }, trimRatio: 0.5);
    }

    #endregion

    #region WinsorizedMean Tests

    [TestMethod]
    public void WinsorizedMean_NoWinsorize_ShouldReturnArithmeticMean()
    {
        var data = new double[] { 1, 2, 3, 4, 5 };
        double result = WinsorizedMean.Compute(data, trimRatio: 0.0);
        Assert.AreEqual(3.0, result, 1e-10);
    }

    [TestMethod]
    public void WinsorizedMean_WithOutliers_ShouldResist()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 100 };
        double ordinaryMean = 14.5;
        double winsorized = WinsorizedMean.Compute(data, trimRatio: 0.1);
        Assert.IsTrue(Math.Abs(winsorized - 5.0) < Math.Abs(ordinaryMean - 5.0),
            $"Winsorized mean {winsorized} should be closer to 5 than ordinary mean {ordinaryMean}");
    }

    [TestMethod]
    public void WinsorizedMean_PreservesSampleSize()
    {
        // Winsorizing replaces extremes rather than removing them,
        // so the effective sample size is unchanged.
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double winsorized = WinsorizedMean.Compute(data, trimRatio: 0.2);
        // With 20% winsorization the two lowest become 3, two highest become 8
        // Expected: (3+3+3+4+5+6+7+8+8+8)/10 = 5.5
        Assert.AreEqual(5.5, winsorized, 1e-10);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentOutOfRangeException))]
    public void WinsorizedMean_InvalidRatio_ShouldThrow()
    {
        WinsorizedMean.Compute(new double[] { 1, 2, 3 }, trimRatio: -0.1);
    }

    #endregion

    #region HuberLoss Tests

    [TestMethod]
    public void HuberLoss_SmallResidual_ShouldBeQuadratic()
    {
        double delta = 1.35;
        double r = 0.5;
        double expected = 0.5 * r * r;
        Assert.AreEqual(expected, HuberLoss.Loss(r, delta), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_LargeResidual_ShouldBeLinear()
    {
        double delta = 1.35;
        double r = 5.0;
        double expected = delta * (Math.Abs(r) - 0.5 * delta);
        Assert.AreEqual(expected, HuberLoss.Loss(r, delta), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_AtThreshold_ShouldBeContinuous()
    {
        double delta = 2.0;
        double quadratic = 0.5 * delta * delta;
        double linear = delta * (delta - 0.5 * delta);
        Assert.AreEqual(quadratic, linear, 1e-10, "Huber loss should be continuous at threshold");
        Assert.AreEqual(quadratic, HuberLoss.Loss(delta, delta), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_Gradient_SmallResidual_ShouldEqualResidual()
    {
        Assert.AreEqual(0.5, HuberLoss.Gradient(0.5, 1.35), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_Gradient_LargeResidual_ShouldBeClipped()
    {
        Assert.AreEqual(1.35, HuberLoss.Gradient(10.0, 1.35), 1e-10);
        Assert.AreEqual(-1.35, HuberLoss.Gradient(-10.0, 1.35), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_MeanLoss_ShouldBeAverage()
    {
        double[] residuals = { 0.5, -0.5, 3.0, -3.0 };
        double delta = 1.35;
        double expected = (HuberLoss.Loss(0.5, delta) + HuberLoss.Loss(-0.5, delta)
                         + HuberLoss.Loss(3.0, delta) + HuberLoss.Loss(-3.0, delta)) / 4.0;
        Assert.AreEqual(expected, HuberLoss.MeanLoss(residuals, delta), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_Weight_SmallResidual_ShouldBeOne()
    {
        Assert.AreEqual(1.0, HuberLoss.Weight(0.5, 1.35), 1e-10);
    }

    [TestMethod]
    public void HuberLoss_Weight_LargeResidual_ShouldBeLessThanOne()
    {
        double w = HuberLoss.Weight(5.0, 1.35);
        Assert.IsTrue(w > 0 && w < 1.0, $"Weight should be in (0,1), got {w}");
        Assert.AreEqual(1.35 / 5.0, w, 1e-10);
    }

    #endregion

    #region Ransac Tests

    [TestMethod]
    public void Ransac_PerfectLine_ShouldFitExactly()
    {
        var x = new double[] { 1, 2, 3, 4, 5 };
        var y = new double[] { 3, 5, 7, 9, 11 }; // y = 1 + 2x

        var result = Ransac.FitLine(x, y, residualThreshold: 0.1);
        Assert.AreEqual(1.0, result.BestModel[0], 1e-6, "Intercept should be 1");
        Assert.AreEqual(2.0, result.BestModel[1], 1e-6, "Slope should be 2");
        Assert.AreEqual(5, result.InlierCount);
    }

    [TestMethod]
    public void Ransac_WithOutliers_ShouldResist()
    {
        // True line: y = 2x
        var x = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        var y = new double[] { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };

        // Add 3 outliers
        y[7] = 100;
        y[8] = -50;
        y[9] = 200;

        var result = Ransac.FitLine(x, y, residualThreshold: 1.0, maxIterations: 2000, seed: 42);

        // Slope should be close to 2
        Assert.AreEqual(2.0, result.BestModel[1], 0.5, $"Slope should be near 2, got {result.BestModel[1]}");
        Assert.IsTrue(result.InlierCount >= 7, $"Should have at least 7 inliers, got {result.InlierCount}");
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Ransac_TooFewPoints_ShouldThrow()
    {
        Ransac.FitLine(new double[] { 1 }, new double[] { 1 });
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Ransac_MismatchedLengths_ShouldThrow()
    {
        Ransac.FitLine(new double[] { 1, 2, 3 }, new double[] { 1, 2 });
    }

    #endregion

    #region OutlierDetection Tests

    [TestMethod]
    public void OutlierDetection_IQR_ShouldDetectExtremes()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100 };
        var result = OutlierDetection.Iqr(data, k: 1.5);
        Assert.IsTrue(result.OutlierMask[10], "100 should be detected as outlier");
        Assert.IsTrue(result.OutlierCount >= 1);
    }

    [TestMethod]
    public void OutlierDetection_IQR_NoOutliers_ShouldReturnNone()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        var result = OutlierDetection.Iqr(data, k: 1.5);
        Assert.AreEqual(0, result.OutlierCount);
    }

    [TestMethod]
    public void OutlierDetection_ZScore_ShouldDetectExtremes()
    {
        var rng = new Random(42);
        var data = new double[100];
        for (int i = 0; i < 98; i++) data[i] = rng.NextDouble() * 2 - 1;
        data[98] = 50.0;
        data[99] = -50.0;

        var result = OutlierDetection.ZScore(data, threshold: 3.0);
        Assert.IsTrue(result.OutlierMask[98], "50 should be an outlier");
        Assert.IsTrue(result.OutlierMask[99], "-50 should be an outlier");
    }

    [TestMethod]
    public void OutlierDetection_ModifiedZScore_ShouldDetectExtremes()
    {
        var data = new double[50];
        for (int i = 0; i < 48; i++) data[i] = i * 0.1;
        data[48] = 500.0;
        data[49] = -500.0;

        var result = OutlierDetection.ModifiedZScore(data, threshold: 3.5);
        Assert.IsTrue(result.OutlierMask[48], "500 should be an outlier");
        Assert.IsTrue(result.OutlierMask[49], "-500 should be an outlier");
    }

    [TestMethod]
    public void OutlierDetection_IndicesMatchMask()
    {
        var data = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 100 };
        var result = OutlierDetection.Iqr(data, k: 1.5);

        int maskCount = 0;
        for (int i = 0; i < result.OutlierMask.Length; i++)
            if (result.OutlierMask[i]) maskCount++;

        Assert.AreEqual(result.OutlierCount, maskCount);
        Assert.AreEqual(result.OutlierCount, result.OutlierIndices.Length);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void OutlierDetection_IQR_TooFewPoints_ShouldThrow()
    {
        OutlierDetection.Iqr(new double[] { 1, 2, 3 });
    }

    #endregion
}
