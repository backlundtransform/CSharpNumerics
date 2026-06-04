using System;
using System.Collections.Generic;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.StateEstimation;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace NumericsTests;

[TestClass]
public class StateEstimationTests
{
    // Box–Muller gaussian noise from a seeded RNG (deterministic tests).
    private static double NextGaussian(Random rng, double mean, double std)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        double z = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        return mean + std * z;
    }

    #region KalmanFilter

    [TestMethod]
    public void KalmanFilter_ConstantSignalWithNoise_ConvergesToTrueLevel()
    {
        const double trueLevel = 10.0;
        var rng = new Random(42);

        var kf = new KalmanFilter(
            new VectorN(new[] { 0.0 }),          // start far from the truth
            new Matrix(new[,] { { 1.0 } }));

        var F = new Matrix(new[,] { { 1.0 } });
        var Q = new Matrix(new[,] { { 1e-5 } }); // level is (almost) constant
        var H = new Matrix(new[,] { { 1.0 } });
        var R = new Matrix(new[,] { { 0.25 } }); // measurement variance

        for (int i = 0; i < 500; i++)
        {
            double z = NextGaussian(rng, trueLevel, 0.5);
            kf.Predict(F, Q);
            kf.Update(H, R, new VectorN(new[] { z }));
        }

        Assert.AreEqual(trueLevel, kf.State[0], 0.2,
            "Filter should converge to the true constant level.");
        Assert.IsTrue(kf.Covariance.values[0, 0] < 0.05,
            "Covariance should shrink as measurements accumulate.");
    }

    [TestMethod]
    public void KalmanFilter_TracksLinearRamp_WithCorrectVelocity()
    {
        const double trueVelocity = 2.0;
        var rng = new Random(7);

        // State = [position, velocity], constant-velocity model with dt = 1.
        var kf = new KalmanFilter(
            new VectorN(new[] { 0.0, 0.0 }),
            new Matrix(new[,] { { 10.0, 0.0 }, { 0.0, 10.0 } }));

        var F = new Matrix(new[,] { { 1.0, 1.0 }, { 0.0, 1.0 } });
        var Q = new Matrix(new[,] { { 1e-4, 0.0 }, { 0.0, 1e-4 } });
        var H = new Matrix(new[,] { { 1.0, 0.0 } });   // observe position only
        var R = new Matrix(new[,] { { 1.0 } });

        for (int t = 0; t < 300; t++)
        {
            double z = NextGaussian(rng, trueVelocity * t, 1.0);
            kf.Predict(F, Q);
            kf.Update(H, R, new VectorN(new[] { z }));
        }

        Assert.AreEqual(trueVelocity, kf.State[1], 0.1,
            "Estimated velocity should match the ramp slope.");
    }

    #endregion

    #region ExtendedKalmanFilter

    [TestMethod]
    public void ExtendedKalmanFilter_EstimatesSinusoid_WithKnownModel()
    {
        // State = phase θ. Process: θ_{k+1} = θ_k + ω (known ω) — linear.
        // Measurement: z = sin(θ) — non-linear, so the Jacobian H = cos(θ) is needed.
        const double omega = 0.1;
        const double measNoiseStd = 0.1;
        var rng = new Random(123);

        var ekf = new ExtendedKalmanFilter(
            new VectorN(new[] { 0.0 }),
            new Matrix(new[,] { { 0.1 } }));

        Func<VectorN, VectorN> f = x => new VectorN(new[] { x[0] + omega });
        Func<VectorN, Matrix> jacF = _ => new Matrix(new[,] { { 1.0 } });
        Func<VectorN, VectorN> h = x => new VectorN(new[] { Math.Sin(x[0]) });
        Func<VectorN, Matrix> jacH = x => new Matrix(new[,] { { Math.Cos(x[0]) } });

        var Q = new Matrix(new[,] { { 1e-5 } });
        var R = new Matrix(new[,] { { measNoiseStd * measNoiseStd } });

        int n = 200;
        double sseModel = 0.0, sseRaw = 0.0;

        for (int k = 0; k < n; k++)
        {
            double thetaTrue = omega * k;
            double clean = Math.Sin(thetaTrue);
            double z = clean + NextGaussian(rng, 0.0, measNoiseStd);

            ekf.Predict(f, jacF, Q);
            ekf.Update(h, jacH, R, new VectorN(new[] { z }));

            double reconstructed = Math.Sin(ekf.State[0]);
            sseModel += (reconstructed - clean) * (reconstructed - clean);
            sseRaw += (z - clean) * (z - clean);
        }

        double rmseModel = Math.Sqrt(sseModel / n);
        double rmseRaw = Math.Sqrt(sseRaw / n);

        Assert.IsTrue(rmseModel < measNoiseStd,
            $"EKF reconstruction RMSE ({rmseModel:F4}) should be below the measurement noise.");
        Assert.IsTrue(rmseModel < rmseRaw,
            "EKF should track the clean sinusoid better than the raw noisy measurement.");
    }

    #endregion

    #region KalmanSmoother

    [TestMethod]
    public void KalmanSmoother_HasLowerVariance_ThanForwardOnlyKalman()
    {
        const double trueLevel = 5.0;
        var rng = new Random(99);

        var F = new Matrix(new[,] { { 1.0 } });
        var Q = new Matrix(new[,] { { 1e-4 } });
        var H = new Matrix(new[,] { { 1.0 } });
        var R = new Matrix(new[,] { { 0.5 } });

        var measurements = new List<VectorN>();
        for (int i = 0; i < 100; i++)
            measurements.Add(new VectorN(new[] { NextGaussian(rng, trueLevel, Math.Sqrt(0.5)) }));

        var smoother = new KalmanSmoother();
        var result = smoother.Smooth(
            measurements,
            new VectorN(new[] { 0.0 }),
            new Matrix(new[,] { { 1.0 } }),
            F, Q, H, R);

        // At an interior index the smoothed covariance must not exceed the filtered one.
        int mid = measurements.Count / 2;
        double filteredVar = result.FilteredCovariances[mid].values[0, 0];
        double smoothedVar = result.SmoothedCovariances[mid].values[0, 0];

        Assert.IsTrue(smoothedVar < filteredVar,
            $"Smoothed variance ({smoothedVar:E3}) should be below the forward-only variance ({filteredVar:E3}).");

        // The smoothed track should also be closer to the truth on average.
        double sseFilt = 0.0, sseSmooth = 0.0;
        for (int k = 0; k < measurements.Count; k++)
        {
            sseFilt += Math.Pow(result.FilteredStates[k][0] - trueLevel, 2);
            sseSmooth += Math.Pow(result.SmoothedStates[k][0] - trueLevel, 2);
        }
        Assert.IsTrue(sseSmooth <= sseFilt,
            "Smoothed estimates should be at least as accurate as the filtered ones.");
    }

    #endregion

    #region HoltWinters

    [TestMethod]
    public void HoltWinters_CapturesDailySeason_ForecastErrorBelow5Percent()
    {
        const int season = 24;       // hourly samples, daily cycle
        const double baseline = 100.0;
        var rng = new Random(2024);

        // A clear daily pattern superimposed on a constant baseline.
        double Pattern(int hour) => 15.0 * Math.Sin(2.0 * Math.PI * hour / season);

        int days = 7;
        var data = new double[days * season];
        for (int t = 0; t < data.Length; t++)
            data[t] = baseline + Pattern(t % season) + NextGaussian(rng, 0.0, 1.0);

        var hw = new HoltWintersSmoothing(season, alpha: 0.3, beta: 0.05, gamma: 0.4, SeasonalType.Additive);
        hw.Fit(data);

        double[] forecast = hw.Forecast(season);

        // MAPE against the clean generating signal for the next day.
        double mape = 0.0;
        for (int h = 0; h < season; h++)
        {
            double expected = baseline + Pattern(h);
            mape += Math.Abs((forecast[h] - expected) / expected);
        }
        mape /= season;

        Assert.IsTrue(mape < 0.05,
            $"Daily-season forecast MAPE ({mape:P2}) should be below 5%.");
    }

    [TestMethod]
    public void HoltWinters_RequiresTwoFullSeasons()
    {
        var hw = new HoltWintersSmoothing(seasonLength: 24, alpha: 0.3, beta: 0.1, gamma: 0.3);
        var tooShort = new double[30];  // < 2 × 24

        Assert.ThrowsException<ArgumentException>(() => hw.Fit(tooShort));
    }

    #endregion
}
