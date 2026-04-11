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
}
