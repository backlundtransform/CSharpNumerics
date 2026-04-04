using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;

namespace NumericTest;

[TestClass]
public class TimeSeriesAnalysisTests
{
    // ──────────────────────────────────────────────────
    // Lomb–Scargle Periodogram
    // ──────────────────────────────────────────────────

    [TestMethod]
    public void LombScargle_DetectsSingleFrequency()
    {
        // Sine wave at 0.1 Hz sampled unevenly
        var rng = new Random(42);
        int n = 200;
        double freq = 0.1; // Hz
        double[] t = new double[n];
        double[] y = new double[n];
        double time = 0;
        for (int i = 0; i < n; i++)
        {
            time += 0.8 + 0.4 * rng.NextDouble(); // ~1.0 ± 0.2 spacing
            t[i] = time;
            y[i] = 3.0 * Math.Sin(2 * Math.PI * freq * t[i]);
        }

        var result = LombScarglePeriodogram.Compute(
            new VectorN(t), new VectorN(y),
            minFrequency: 0.01, maxFrequency: 0.5, numFrequencies: 500);

        Assert.AreEqual(500, result.Frequencies.Length);
        // Best frequency should be near 0.1 Hz
        Assert.AreEqual(freq, result.BestFrequency, 0.01);
        Assert.IsTrue(result.BestPower > 0);
    }

    [TestMethod]
    public void LombScargle_ByPeriod_DetectsPeriod()
    {
        double period = 5.0;
        int n = 150;
        var rng = new Random(123);
        double[] t = new double[n];
        double[] y = new double[n];
        double time = 0;
        for (int i = 0; i < n; i++)
        {
            time += 0.3 + 0.2 * rng.NextDouble();
            t[i] = time;
            y[i] = 2.0 * Math.Sin(2 * Math.PI * t[i] / period);
        }

        var result = LombScarglePeriodogram.ComputeByPeriod(
            new VectorN(t), new VectorN(y),
            minPeriod: 1, maxPeriod: 20, numPeriods: 500);

        Assert.AreEqual(period, result.BestPeriod, 0.3);
    }

    [TestMethod]
    public void LombScargle_FalseAlarmProbability()
    {
        double fap = LombScarglePeriodogram.FalseAlarmProbability(10.0, 100, 50);
        Assert.IsTrue(fap < 0.01); // high power → low FAP

        double fap2 = LombScarglePeriodogram.FalseAlarmProbability(1.0, 100, 50);
        Assert.IsTrue(fap2 > fap); // lower power → higher FAP
    }

    [TestMethod]
    public void LombScargle_AngularFreqOverload()
    {
        double[] t = { 0, 1, 2, 3, 4, 5 };
        double[] y = { 0, 1, 0, -1, 0, 1 };
        double[] omega = new double[10];
        for (int i = 0; i < 10; i++)
            omega[i] = 0.5 + i * 0.5;

        var result = LombScarglePeriodogram.Compute(
            new VectorN(t), new VectorN(y), new VectorN(omega));

        Assert.AreEqual(10, result.Power.Length);
        Assert.IsTrue(result.BestPower > 0);
    }

    // ──────────────────────────────────────────────────
    // Box-fitting Least Squares
    // ──────────────────────────────────────────────────

    [TestMethod]
    public void BLS_DetectsTransitPeriod()
    {
        // Simulate a transit with period=10, duration fraction ~0.05, depth 0.02
        double truePeriod = 10.0;
        double depth = 0.02;
        double durFrac = 0.05;
        int n = 500;
        var rng = new Random(99);

        double[] t = new double[n];
        double[] y = new double[n];
        double time = 0;
        for (int i = 0; i < n; i++)
        {
            time += 0.18 + 0.04 * rng.NextDouble();
            t[i] = time;
            double phase = ((t[i] / truePeriod) % 1.0 + 1.0) % 1.0;
            y[i] = 1.0 + 0.001 * rng.NextDouble(); // baseline ~ 1
            if (phase < durFrac)
                y[i] -= depth; // transit dip
        }

        var result = BoxFittingLeastSquares.Compute(
            new VectorN(t), new VectorN(y),
            minPeriod: 5, maxPeriod: 20, numPeriods: 300,
            minDurationFraction: 0.01, maxDurationFraction: 0.15);

        Assert.AreEqual(300, result.Periods.Length);
        Assert.AreEqual(truePeriod, result.BestPeriod, 0.5);
        Assert.IsTrue(result.BestDepth > 0);
    }

    [TestMethod]
    public void BLS_WithExplicitPeriods()
    {
        double truePeriod = 5.0;
        int n = 300;
        var rng = new Random(55);
        double[] t = new double[n];
        double[] y = new double[n];
        double time = 0;
        for (int i = 0; i < n; i++)
        {
            time += 0.09 + 0.02 * rng.NextDouble();
            t[i] = time;
            double phase = ((t[i] / truePeriod) % 1.0 + 1.0) % 1.0;
            y[i] = 1.0;
            if (phase < 0.08) y[i] -= 0.03;
        }

        double[] periods = new double[100];
        for (int i = 0; i < 100; i++)
            periods[i] = 2.0 + i * 0.1;

        var result = BoxFittingLeastSquares.Compute(
            new VectorN(t), new VectorN(y), new VectorN(periods));

        Assert.AreEqual(truePeriod, result.BestPeriod, 0.5);
    }

    // ──────────────────────────────────────────────────
    // Phase Folding
    // ──────────────────────────────────────────────────

    [TestMethod]
    public void PhaseFold_ProducesCorrectPhases()
    {
        var times = new VectorN(new[] { 0.0, 2.5, 5.0, 7.5, 10.0, 12.5 });
        var values = new VectorN(new[] { 1.0, 2.0, 1.0, 2.0, 1.0, 2.0 });

        var (phases, vals) = PhaseFolding.Fold(times, values, period: 5.0);

        Assert.AreEqual(6, phases.Length);
        // All phases should be in [0, 1)
        for (int i = 0; i < phases.Length; i++)
        {
            Assert.IsTrue(phases[i] >= 0 && phases[i] < 1.0);
        }
        // Sorted
        for (int i = 1; i < phases.Length; i++)
            Assert.IsTrue(phases[i] >= phases[i - 1]);
    }

    [TestMethod]
    public void PhaseFold_WithEpoch()
    {
        var times = new VectorN(new[] { 1.0, 3.0, 5.0, 7.0 });
        var values = new VectorN(new[] { 10.0, 20.0, 10.0, 20.0 });

        var (phases, vals) = PhaseFolding.Fold(times, values, period: 4.0, epoch: 1.0);

        Assert.AreEqual(4, phases.Length);
        // t=1 with epoch=1 → phase=0, t=3 → phase=0.5
        Assert.AreEqual(0.0, phases[0], 1e-10);
    }

    [TestMethod]
    public void PhaseFoldAndBin_ProducesCorrectBins()
    {
        // Period=2, 100 points
        int n = 100;
        double[] t = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            t[i] = i * 0.1; // 0 to 9.9
            y[i] = Math.Sin(2 * Math.PI * t[i] / 2.0); // period=2
        }

        var (centers, means) = PhaseFolding.FoldAndBin(
            new VectorN(t), new VectorN(y), period: 2.0, numBins: 10);

        Assert.AreEqual(10, centers.Length);
        Assert.AreEqual(0.05, centers[0], 1e-10); // first bin center
        Assert.AreEqual(0.95, centers[9], 1e-10); // last bin center

        // No NaN bins
        for (int i = 0; i < 10; i++)
            Assert.IsFalse(double.IsNaN(means[i]));
    }

    [TestMethod]
    public void PhaseFoldAndBinWithErrors_ReturnsStdDevs()
    {
        int n = 200;
        double[] t = new double[n];
        double[] y = new double[n];
        var rng = new Random(7);
        for (int i = 0; i < n; i++)
        {
            t[i] = i * 0.05;
            y[i] = Math.Cos(2 * Math.PI * t[i] / 3.0) + 0.1 * rng.NextDouble();
        }

        var (centers, means, stds) = PhaseFolding.FoldAndBinWithErrors(
            new VectorN(t), new VectorN(y), period: 3.0, numBins: 15);

        Assert.AreEqual(15, stds.Length);
        for (int i = 0; i < 15; i++)
        {
            Assert.IsFalse(double.IsNaN(stds[i]));
            Assert.IsTrue(stds[i] >= 0);
        }
    }

    // ──────────────────────────────────────────────────
    // Detrending
    // ──────────────────────────────────────────────────

    [TestMethod]
    public void Detrend_Polynomial_RemovesLinearTrend()
    {
        int n = 50;
        double[] t = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            t[i] = i;
            y[i] = 2.0 * i + 5.0 + Math.Sin(2 * Math.PI * i / 10.0);
        }

        var (detrended, trend) = TimeSeriesDetrending.PolynomialDetrend(
            new VectorN(t), new VectorN(y), degree: 1);

        Assert.AreEqual(n, detrended.Length);
        // After removing the linear trend, the detrended signal should oscillate around 0
        double mean = 0;
        for (int i = 0; i < n; i++) mean += detrended[i];
        mean /= n;
        Assert.AreEqual(0, mean, 0.5);
    }

    [TestMethod]
    public void Detrend_Polynomial_QuadraticTrend()
    {
        int n = 80;
        double[] t = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            t[i] = i * 0.1;
            y[i] = 0.5 * t[i] * t[i] - 2 * t[i] + 3; // pure quadratic
        }

        var (detrended, trend) = TimeSeriesDetrending.PolynomialDetrend(
            new VectorN(t), new VectorN(y), degree: 2);

        // Should remove the quadratic perfectly
        for (int i = 0; i < n; i++)
            Assert.AreEqual(0, detrended[i], 1e-6);
    }

    [TestMethod]
    public void Detrend_MovingAverage_SmoothsSignal()
    {
        int n = 100;
        double[] t = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            t[i] = i;
            y[i] = 10.0 + 0.1 * i + Math.Sin(2 * Math.PI * i / 5.0);
        }

        var (detrended, trend) = TimeSeriesDetrending.MovingAverageDetrend(
            new VectorN(y), windowSize: 5);

        Assert.AreEqual(n, detrended.Length);
        Assert.AreEqual(n, trend.Length);
        // Trend should be smoother than original
    }

    [TestMethod]
    public void Detrend_MedianFilter_RobustToOutliers()
    {
        int n = 50;
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
            y[i] = 5.0; // flat signal
        y[25] = 100.0; // single outlier

        var (detrended, trend) = TimeSeriesDetrending.MedianFilterDetrend(
            new VectorN(y), windowSize: 7);

        // Median filter should not be much affected by the single outlier
        Assert.AreEqual(5.0, trend[20], 1e-10);
        Assert.AreEqual(5.0, trend[30], 1e-10);
    }

    [TestMethod]
    public void Detrend_SavitzkyGolay_PreservesShape()
    {
        int n = 100;
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
            y[i] = 3.0 * i + 1.0; // perfect linear → SG should reproduce exactly

        var (detrended, trend) = TimeSeriesDetrending.SavitzkyGolayDetrend(
            new VectorN(y), windowSize: 11, polyOrder: 1);

        // For linear data with polyOrder≥1, SG should reconstruct perfectly in interior
        for (int i = 10; i < n - 10; i++)
            Assert.AreEqual(0, detrended[i], 1e-8);
    }

    [TestMethod]
    public void Detrend_Differencing_ReducesOrder()
    {
        // Linear trend: first difference should be constant
        int n = 20;
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
            y[i] = 3.0 * i + 2.0;

        VectorN diff1 = TimeSeriesDetrending.Difference(new VectorN(y), order: 1);
        Assert.AreEqual(n - 1, diff1.Length);
        for (int i = 0; i < diff1.Length; i++)
            Assert.AreEqual(3.0, diff1[i], 1e-10);

        // Second difference of linear → all zeros
        VectorN diff2 = TimeSeriesDetrending.Difference(new VectorN(y), order: 2);
        Assert.AreEqual(n - 2, diff2.Length);
        for (int i = 0; i < diff2.Length; i++)
            Assert.AreEqual(0, diff2[i], 1e-10);
    }

    // ──────────────────────────────────────────────────
    // Peak Fitting
    // ──────────────────────────────────────────────────

    [TestMethod]
    public void FindPeaks_DetectsLocalMaxima()
    {
        // Signal with two peaks
        double[] y = new double[50];
        for (int i = 0; i < 50; i++)
        {
            y[i] = Math.Exp(-((i - 15.0) / 3.0) * ((i - 15.0) / 3.0)) +
                   0.7 * Math.Exp(-((i - 35.0) / 3.0) * ((i - 35.0) / 3.0));
        }

        var peaks = PeakFitting.FindPeaks(new VectorN(y));

        Assert.IsTrue(peaks.Count >= 2);
        Assert.IsTrue(peaks.Contains(15));
        Assert.IsTrue(peaks.Contains(35));
    }

    [TestMethod]
    public void FindPeaks_ProminenceFilter()
    {
        double[] y = new double[50];
        for (int i = 0; i < 50; i++)
        {
            y[i] = Math.Exp(-((i - 15.0) / 3.0) * ((i - 15.0) / 3.0)) +
                   0.05 * Math.Exp(-((i - 35.0) / 1.0) * ((i - 35.0) / 1.0));
        }

        // Without prominence filter
        var allPeaks = PeakFitting.FindPeaks(new VectorN(y));
        // With high prominence filter — should keep only the big peak
        var bigPeaks = PeakFitting.FindPeaks(new VectorN(y), prominence: 0.5);

        Assert.IsTrue(allPeaks.Count > bigPeaks.Count);
        Assert.AreEqual(1, bigPeaks.Count);
        Assert.AreEqual(15, bigPeaks[0]);
    }

    [TestMethod]
    public void FindValleys_DetectsLocalMinima()
    {
        double[] y = new double[30];
        for (int i = 0; i < 30; i++)
            y[i] = -Math.Exp(-((i - 15.0) / 4.0) * ((i - 15.0) / 4.0));

        var valleys = PeakFitting.FindValleys(new VectorN(y));

        Assert.IsTrue(valleys.Count >= 1);
        Assert.AreEqual(15, valleys[0]);
    }

    [TestMethod]
    public void FitGaussian_RecoversParameters()
    {
        // Generate Gaussian: A=5, μ=10, σ=2, baseline=1
        double A = 5.0, mu = 10.0, sigma = 2.0, baseline = 1.0;
        int n = 50;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i * 0.5;
            y[i] = A * Math.Exp(-0.5 * ((x[i] - mu) / sigma) * ((x[i] - mu) / sigma)) + baseline;
        }

        PeakResult result = PeakFitting.FitGaussian(
            new VectorN(x), new VectorN(y),
            centerEstimate: 10.0, widthEstimate: 2.5);

        Assert.AreEqual(mu, result.Center, 0.2);
        Assert.AreEqual(A, result.Amplitude, 0.5);
        Assert.AreEqual(sigma, result.Width, 0.5);
        Assert.AreEqual(baseline, result.Baseline, 0.5);
        Assert.AreEqual(PeakShape.Gaussian, result.Shape);
        Assert.IsTrue(result.RSquared > 0.99);
    }

    [TestMethod]
    public void FitLorentzian_RecoversParameters()
    {
        // Generate Lorentzian: A=4, x₀=8, γ=1.5, baseline=0.5
        double A = 4.0, x0 = 8.0, gamma = 1.5, baseline = 0.5;
        int n = 80;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i * 0.25;
            y[i] = A * (gamma * gamma) / ((x[i] - x0) * (x[i] - x0) + gamma * gamma) + baseline;
        }

        PeakResult result = PeakFitting.FitLorentzian(
            new VectorN(x), new VectorN(y),
            centerEstimate: 8.0, widthEstimate: 2.0);

        Assert.AreEqual(x0, result.Center, 0.3);
        Assert.AreEqual(A, result.Amplitude, 0.5);
        Assert.AreEqual(gamma, result.Width, 0.5);
        Assert.AreEqual(PeakShape.Lorentzian, result.Shape);
        Assert.IsTrue(result.RSquared > 0.99);
    }

    [TestMethod]
    public void FindAndFitPeaks_GaussianShape()
    {
        // Two Gaussian peaks
        int n = 100;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i;
            y[i] = 3.0 * Math.Exp(-0.5 * ((x[i] - 30) / 4.0) * ((x[i] - 30) / 4.0)) +
                   2.0 * Math.Exp(-0.5 * ((x[i] - 70) / 3.0) * ((x[i] - 70) / 3.0));
        }

        var results = PeakFitting.FindAndFitPeaks(
            new VectorN(x), new VectorN(y),
            prominence: 0.5,
            shape: PeakShape.Gaussian,
            fitRadius: 15);

        Assert.AreEqual(2, results.Count);

        // First peak near x=30
        Assert.AreEqual(30, results[0].Center, 2.0);
        Assert.AreEqual(3.0, results[0].Amplitude, 1.0);

        // Second peak near x=70
        Assert.AreEqual(70, results[1].Center, 2.0);
    }

    [TestMethod]
    public void FindAndFitPeaks_MinDistance()
    {
        // Two peaks close together
        int n = 50;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i;
            y[i] = 3.0 * Math.Exp(-0.5 * ((x[i] - 20) / 2.0) * ((x[i] - 20) / 2.0)) +
                   2.0 * Math.Exp(-0.5 * ((x[i] - 24) / 2.0) * ((x[i] - 24) / 2.0));
        }

        var resultsClose = PeakFitting.FindAndFitPeaks(
            new VectorN(x), new VectorN(y),
            prominence: 0.1, shape: PeakShape.Gaussian, fitRadius: 8, minDistance: 1);

        var resultsFar = PeakFitting.FindAndFitPeaks(
            new VectorN(x), new VectorN(y),
            prominence: 0.1, shape: PeakShape.Gaussian, fitRadius: 8, minDistance: 10);

        Assert.IsTrue(resultsClose.Count >= resultsFar.Count);
    }
}
