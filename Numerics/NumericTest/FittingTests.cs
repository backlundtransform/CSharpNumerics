using Microsoft.VisualStudio.TestTools.UnitTesting;
using CSharpNumerics.Statistics.Fitting;
using CSharpNumerics.Numerics.Objects;
using System;

namespace NumericTest;

[TestClass]
public class FittingTests
{
    private const double Tol = 1e-6;

    #region Least Squares — Linear

    [TestMethod]
    public void LinearFit_ExactLine()
    {
        // y = 2 + 3x
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var y = new VectorN(new[] { 5.0, 8, 11, 14, 17 });

        var result = LeastSquaresFitter.Fit(x, y, degree: 1);

        Assert.AreEqual(2.0, result.Coefficients[0], Tol, "Intercept");
        Assert.AreEqual(3.0, result.Coefficients[1], Tol, "Slope");
        Assert.AreEqual(1.0, result.RSquared, Tol, "R²");
        Assert.AreEqual(0.0, result.RMSE, Tol, "RMSE");
        Assert.AreEqual(5, result.ObservationCount);
        Assert.AreEqual(2, result.ParameterCount);
        Assert.AreEqual(3, result.DegreesOfFreedom);
    }

    [TestMethod]
    public void LinearFit_NoisyData()
    {
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
        var y = new VectorN(new[] { 2.1, 4.0, 5.9, 8.1, 10.0, 11.9, 14.1, 16.0, 17.9, 20.1 });

        var result = LeastSquaresFitter.Fit(x, y);

        Assert.IsTrue(result.RSquared > 0.999, "R² should be near 1");
        Assert.AreEqual(10, result.Residuals.Length);
        Assert.IsTrue(result.RMSE < 0.1);
    }

    #endregion

    #region Least Squares — Polynomial

    [TestMethod]
    public void QuadraticFit_ExactParabola()
    {
        // y = 1 + 2x + 3x²
        double[] xArr = { -2, -1, 0, 1, 2, 3 };
        double[] yArr = new double[6];
        for (int i = 0; i < xArr.Length; i++)
            yArr[i] = 1 + 2 * xArr[i] + 3 * xArr[i] * xArr[i];

        var result = LeastSquaresFitter.Fit(new VectorN(xArr), new VectorN(yArr), degree: 2);

        Assert.AreEqual(1.0, result.Coefficients[0], Tol, "a0");
        Assert.AreEqual(2.0, result.Coefficients[1], Tol, "a1");
        Assert.AreEqual(3.0, result.Coefficients[2], Tol, "a2");
        Assert.AreEqual(1.0, result.RSquared, Tol);
    }

    [TestMethod]
    public void CubicFit()
    {
        // y = x³ − 2x + 1
        double[] xArr = { -3, -2, -1, 0, 1, 2, 3, 4 };
        double[] yArr = new double[xArr.Length];
        for (int i = 0; i < xArr.Length; i++)
            yArr[i] = xArr[i] * xArr[i] * xArr[i] - 2 * xArr[i] + 1;

        var result = LeastSquaresFitter.Fit(new VectorN(xArr), new VectorN(yArr), degree: 3);

        Assert.AreEqual(1.0, result.Coefficients[0], Tol, "a0");
        Assert.AreEqual(-2.0, result.Coefficients[1], Tol, "a1");
        Assert.AreEqual(0.0, result.Coefficients[2], Tol, "a2");
        Assert.AreEqual(1.0, result.Coefficients[3], Tol, "a3");
    }

    #endregion

    #region Least Squares — Multiple Regression

    [TestMethod]
    public void MultipleRegression_TwoPredictors()
    {
        // y = 1 + 2*x1 + 3*x2
        int n = 10;
        double[,] X = new double[n, 3];
        double[] yArr = new double[n];
        var rng = new Random(42);
        for (int i = 0; i < n; i++)
        {
            double x1 = rng.NextDouble() * 10;
            double x2 = rng.NextDouble() * 10;
            X[i, 0] = 1.0;
            X[i, 1] = x1;
            X[i, 2] = x2;
            yArr[i] = 1.0 + 2.0 * x1 + 3.0 * x2;
        }

        var result = LeastSquaresFitter.Fit(X, new VectorN(yArr));

        Assert.AreEqual(1.0, result.Coefficients[0], 1e-4, "Intercept");
        Assert.AreEqual(2.0, result.Coefficients[1], 1e-4, "β1");
        Assert.AreEqual(3.0, result.Coefficients[2], 1e-4, "β2");
        Assert.AreEqual(1.0, result.RSquared, 1e-4);
    }

    #endregion

    #region Weighted Least Squares

    [TestMethod]
    public void WeightedFit_HighWeightsDominant()
    {
        // Points along y = 2x with one outlier at (3, 100)
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var y = new VectorN(new[] { 2.0, 4, 100, 8, 10 });

        // Give low weight to the outlier
        var weights = new VectorN(new[] { 1.0, 1, 0.001, 1, 1 });

        var result = WeightedLeastSquaresFitter.Fit(x, y, weights);

        // The fit should be close to y ≈ 2x
        Assert.IsTrue(Math.Abs(result.Coefficients[1] - 2.0) < 0.5, "Slope near 2");
    }

    [TestMethod]
    public void WeightedFit_UniformWeights_MatchesOLS()
    {
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var y = new VectorN(new[] { 2.1, 3.9, 6.2, 7.8, 10.1 });
        var w = new VectorN(new[] { 1.0, 1, 1, 1, 1 });

        var wls = WeightedLeastSquaresFitter.Fit(x, y, w);
        var ols = LeastSquaresFitter.Fit(x, y);

        for (int i = 0; i < ols.Coefficients.Length; i++)
            Assert.AreEqual(ols.Coefficients[i], wls.Coefficients[i], 1e-10);
    }

    #endregion

    #region Nonlinear Least Squares

    [TestMethod]
    public void NonlinearFit_Exponential()
    {
        // y = 2 * exp(0.5 * x)
        double[] xArr = { 0, 1, 2, 3, 4, 5 };
        double[] yArr = new double[6];
        for (int i = 0; i < xArr.Length; i++)
            yArr[i] = 2.0 * Math.Exp(0.5 * xArr[i]);

        Func<double, VectorN, double> model = (xi, p) => p[0] * Math.Exp(p[1] * xi);
        var init = new VectorN(new[] { 1.0, 1.0 });

        var result = NonlinearLeastSquaresFitter.Fit(model, new VectorN(xArr), new VectorN(yArr), init);

        Assert.AreEqual(2.0, result.Coefficients[0], 0.01, "Amplitude");
        Assert.AreEqual(0.5, result.Coefficients[1], 0.01, "Rate");
        Assert.IsTrue(result.RSquared > 0.999);
    }

    [TestMethod]
    public void NonlinearFit_Logistic()
    {
        // y = L / (1 + exp(-k*(x - x0)))  with L=10, k=1, x0=5
        double L = 10, k = 1, x0 = 5;
        double[] xArr = new double[20];
        double[] yArr = new double[20];
        for (int i = 0; i < 20; i++)
        {
            xArr[i] = i * 0.5;
            yArr[i] = L / (1.0 + Math.Exp(-k * (xArr[i] - x0)));
        }

        Func<double, VectorN, double> model = (xi, p) =>
            p[0] / (1.0 + Math.Exp(-p[1] * (xi - p[2])));

        var init = new VectorN(new[] { 8.0, 0.5, 4.0 });
        var result = NonlinearLeastSquaresFitter.Fit(model,
            new VectorN(xArr), new VectorN(yArr), init);

        Assert.AreEqual(10.0, result.Coefficients[0], 0.1, "L");
        Assert.AreEqual(1.0, result.Coefficients[1], 0.1, "k");
        Assert.AreEqual(5.0, result.Coefficients[2], 0.1, "x0");
    }

    #endregion

    #region Robust Fitting

    [TestMethod]
    public void RobustFit_Huber_ResistsOutliers()
    {
        // True line: y = 1 + 2x with outliers
        double[] yArr = { 3, 5, 7, 9, 11, 13, 15, 17, 19, 21 };
        yArr[2] = 50;
        yArr[7] = -30;
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
        var y = new VectorN(yArr);

        var robust = RobustFitter.Fit(x, y, degree: 1, weightFunction: RobustWeightFunction.Huber);
        var ols = LeastSquaresFitter.Fit(x, y);

        // Robust slope should be closer to 2 than OLS slope
        Assert.IsTrue(
            Math.Abs(robust.Coefficients[1] - 2.0) < Math.Abs(ols.Coefficients[1] - 2.0),
            "Robust fit should be closer to true slope than OLS");
    }

    [TestMethod]
    public void RobustFit_TukeyBisquare()
    {
        double[] yArr = { 3, 5, 7, 9, 11, 13, 15, 17 };
        yArr[3] = 100; // outlier
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5, 6, 7, 8 });
        var y = new VectorN(yArr);

        var result = RobustFitter.Fit(x, y, degree: 1, weightFunction: RobustWeightFunction.TukeyBisquare);

        Assert.IsTrue(Math.Abs(result.Coefficients[1] - 2.0) < 1.0,
            "Tukey bisquare should resist the outlier");
    }

    #endregion

    #region Goodness of Fit

    [TestMethod]
    public void GoodnessOfFit_PerfectFit()
    {
        var obs = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var pred = new VectorN(new[] { 1.0, 2, 3, 4, 5 });

        Assert.AreEqual(1.0, GoodnessOfFit.RSquared(obs, pred), Tol);
        Assert.AreEqual(0.0, GoodnessOfFit.RMSE(obs, pred), Tol);
        Assert.AreEqual(0.0, GoodnessOfFit.MAE(obs, pred), Tol);
        Assert.AreEqual(0.0, GoodnessOfFit.SSE(obs, pred), Tol);
    }

    [TestMethod]
    public void GoodnessOfFit_Metrics()
    {
        var obs = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var pred = new VectorN(new[] { 1.1, 2.2, 2.8, 4.1, 4.8 });

        double r2 = GoodnessOfFit.RSquared(obs, pred);
        double rmse = GoodnessOfFit.RMSE(obs, pred);
        double mae = GoodnessOfFit.MAE(obs, pred);
        double aic = GoodnessOfFit.AIC(obs, pred, 2);
        double bic = GoodnessOfFit.BIC(obs, pred, 2);

        Assert.IsTrue(r2 > 0.95 && r2 < 1.0);
        Assert.IsTrue(rmse > 0);
        Assert.IsTrue(mae > 0);
        Assert.IsTrue(!double.IsNaN(aic));
        Assert.IsTrue(!double.IsNaN(bic));
    }

    [TestMethod]
    public void GoodnessOfFit_AdjustedRSquared()
    {
        double adjR2 = GoodnessOfFit.AdjustedRSquared(0.95, 100, 5);
        Assert.IsTrue(adjR2 < 0.95, "Adjusted R² should be less than R²");
        Assert.IsTrue(adjR2 > 0.94);
    }

    #endregion

    #region Residual Analysis

    [TestMethod]
    public void ResidualAnalysis_BasicProperties()
    {
        var obs = new VectorN(new[] { 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
        var pred = new VectorN(new[] { 1.1, 1.9, 3.2, 3.8, 5.1, 5.9, 7.2, 7.8, 9.1, 9.9 });

        var analysis = GoodnessOfFit.AnalyzeResiduals(obs, pred);

        Assert.AreEqual(10, analysis.StandardizedResiduals.Length);
        Assert.IsTrue(Math.Abs(analysis.MeanResidual) < 0.5);
        Assert.IsTrue(analysis.ResidualStandardDeviation > 0);
        // DW near 2 means no autocorrelation at this scale
        Assert.IsTrue(analysis.DurbinWatsonStatistic > 0 && analysis.DurbinWatsonStatistic < 4);
    }

    #endregion

    #region Parameter Estimation

    [TestMethod]
    public void ConfidenceIntervals_ContainTrueValues()
    {
        // y = 1 + 2x (exact)
        double[] xArr = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double[] yArr = new double[10];
        var rng = new Random(123);
        for (int i = 0; i < 10; i++)
            yArr[i] = 1.0 + 2.0 * xArr[i] + (rng.NextDouble() - 0.5) * 0.5;

        var result = LeastSquaresFitter.Fit(new VectorN(xArr), new VectorN(yArr));
        var (lower, upper) = ParameterEstimation.ConfidenceIntervals(result, 0.95);

        // True intercept ≈ 1 should be within CI
        Assert.IsTrue(lower[0] < 1.5 && upper[0] > 0.5,
            $"Intercept CI [{lower[0]:F3}, {upper[0]:F3}] should contain ~1");
        // True slope ≈ 2 should be within CI
        Assert.IsTrue(lower[1] < 2.2 && upper[1] > 1.8,
            $"Slope CI [{lower[1]:F3}, {upper[1]:F3}] should contain ~2");
    }

    [TestMethod]
    public void PredictionInterval_ContainsObservation()
    {
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var y = new VectorN(new[] { 2.1, 3.9, 6.2, 7.8, 10.1 });

        var result = LeastSquaresFitter.Fit(x, y);
        double[,] designMatrix = FitDesignMatrix(x.Values);

        // Predict at x = 3 → design row = [1, 3]
        var xNew = new VectorN(new[] { 1.0, 3.0 });
        var (lo, hi) = ParameterEstimation.PredictionInterval(result, designMatrix, xNew, 0.95);

        // Observed y[2] = 6.2 should be within the PI
        Assert.IsTrue(lo < 6.2 && hi > 6.2,
            $"PI [{lo:F3}, {hi:F3}] should contain observed 6.2");
    }

    [TestMethod]
    public void MLE_NormalMean()
    {
        // MLE for normal mean with known variance=1: minimise Σ(x_i − μ)² / 2
        double[] data = { 4.8, 5.1, 5.3, 4.9, 5.0 };

        Func<VectorN, double> negLL = theta =>
        {
            double mu = theta[0];
            double sum = 0;
            for (int i = 0; i < data.Length; i++)
                sum += (data[i] - mu) * (data[i] - mu);
            return sum / 2.0;
        };

        VectorN mle = ParameterEstimation.MaximumLikelihoodEstimate(negLL, new VectorN(new[] { 0.0 }));
        double trueMean = 0;
        for (int i = 0; i < data.Length; i++) trueMean += data[i];
        trueMean /= data.Length;

        Assert.AreEqual(trueMean, mle[0], 0.01);
    }

    [TestMethod]
    public void MethodOfMoments_Normal()
    {
        var data = new VectorN(new[] { 1.0, 2, 3, 4, 5 });
        var (mean, variance) = ParameterEstimation.MethodOfMomentsNormal(data);

        Assert.AreEqual(3.0, mean, Tol);
        Assert.AreEqual(2.5, variance, Tol); // sample variance with n-1
    }

    [TestMethod]
    public void MethodOfMoments_Exponential()
    {
        var data = new VectorN(new[] { 0.5, 1.0, 1.5, 2.0, 2.5 });
        double rate = ParameterEstimation.MethodOfMomentsExponential(data);

        // Mean = 1.5 → rate = 1/1.5
        Assert.AreEqual(1.0 / 1.5, rate, Tol);
    }

    [TestMethod]
    public void MethodOfMoments_Gamma()
    {
        var data = new VectorN(new[] { 2.0, 4, 6, 8, 10 });
        var (shape, scale) = ParameterEstimation.MethodOfMomentsGamma(data);

        Assert.IsTrue(shape > 0);
        Assert.IsTrue(scale > 0);
        // Mean ≈ shape * scale
        Assert.AreEqual(6.0, shape * scale, 0.1);
    }

    #endregion

    #region Validation

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Fit_MismatchedLengths_Throws()
    {
        LeastSquaresFitter.Fit(new VectorN(new[] { 1.0, 2 }), new VectorN(new[] { 1.0, 2, 3 }));
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Fit_TooFewPoints_Throws()
    {
        LeastSquaresFitter.Fit(new VectorN(new[] { 1.0 }), new VectorN(new[] { 1.0 }), degree: 1);
    }

    #endregion

    #region FittingResult Properties

    [TestMethod]
    public void FittingResult_AllPropertiesPopulated()
    {
        var x = new VectorN(new[] { 1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
        var y = new VectorN(new[] { 2.0, 4, 6, 8, 10, 12, 14, 16, 18, 20 });
        var result = LeastSquaresFitter.Fit(x, y);

        Assert.AreEqual(10, result.ObservationCount);
        Assert.AreEqual(2, result.ParameterCount);
        Assert.AreEqual(8, result.DegreesOfFreedom);
        Assert.IsTrue(result.AdjustedRSquared <= result.RSquared);
    }

    #endregion

    private static double[,] FitDesignMatrix(double[] x)
    {
        int n = x.Length;
        double[,] X = new double[n, 2];
        for (int i = 0; i < n; i++)
        {
            X[i, 0] = 1.0;
            X[i, 1] = x[i];
        }
        return X;
    }
}
