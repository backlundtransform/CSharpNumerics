using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Enums;
using CSharpNumerics.Numerics.Interpolation;
using CSharpNumerics.Statistics.Data;

namespace NumericTest;

[TestClass]
public class AdvancedInterpolationTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Polynomial Interpolation
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Polynomial_Lagrange_ShouldPassThroughDataPoints()
    {
        double[] x = { 0, 1, 2, 3 };
        double[] y = { 1, 3, 2, 5 };
        var poly = new PolynomialInterpolation(x, y);

        for (int i = 0; i < x.Length; i++)
            Assert.AreEqual(y[i], poly.Evaluate(x[i]), 1e-10,
                $"Should pass through point ({x[i]}, {y[i]})");
    }

    [TestMethod]
    public void Polynomial_Newton_ShouldMatchLagrange()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 1, -1, 3, 0, 2 };
        var poly = new PolynomialInterpolation(x, y);

        double[] testPoints = { 0.5, 1.5, 2.7, 3.3 };
        foreach (double xi in testPoints)
            Assert.AreEqual(poly.Evaluate(xi), poly.EvaluateNewton(xi), 1e-10,
                $"Lagrange and Newton should agree at x={xi}");
    }

    [TestMethod]
    public void Polynomial_Neville_ShouldMatchLagrange()
    {
        double[] x = { 0, 1, 2, 3 };
        double[] y = { 1, 2, 0, 5 };
        var poly = new PolynomialInterpolation(x, y);

        var (val, _) = poly.EvaluateNeville(1.5);
        Assert.AreEqual(poly.Evaluate(1.5), val, 1e-10);
    }

    [TestMethod]
    public void Polynomial_QuadraticFunction_ShouldBeExact()
    {
        // f(x) = 2x² - 3x + 1 — 3 points should be enough
        Func<double, double> f = x => 2 * x * x - 3 * x + 1;
        double[] x = { 0, 1, 2 };
        double[] y = x.Select(f).ToArray();
        var poly = new PolynomialInterpolation(x, y);

        Assert.AreEqual(f(0.5), poly.Evaluate(0.5), 1e-10);
        Assert.AreEqual(f(1.5), poly.Evaluate(1.5), 1e-10);
    }

    [TestMethod]
    public void Polynomial_Degree_ShouldBeNMinus1()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 1, 2, 3, 4, 5 };
        var poly = new PolynomialInterpolation(x, y);
        Assert.AreEqual(4, poly.Degree);
    }

    [TestMethod]
    public void Polynomial_ExtensionMethod_ShouldWork()
    {
        var serie = new List<Serie>
        {
            new() { Index = 0, Value = 1 },
            new() { Index = 1, Value = 3 },
            new() { Index = 2, Value = 2 },
            new() { Index = 3, Value = 5 }
        };

        double val = serie.Interpolate(p => (p.Index, p.Value), 1.5, InterpolationType.Polynomial);

        // Also evaluate directly for comparison
        var poly = new PolynomialInterpolation(
            serie.Select(s => s.Index).ToArray(),
            serie.Select(s => s.Value).ToArray());

        Assert.AreEqual(poly.Evaluate(1.5), val, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Cubic Spline Interpolation
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CubicSpline_Natural_ShouldPassThroughDataPoints()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 0, 1, 0, 1, 0 };
        var spline = new CubicSplineInterpolation(x, y);

        for (int i = 0; i < x.Length; i++)
            Assert.AreEqual(y[i], spline.Evaluate(x[i]), 1e-10,
                $"Spline should pass through ({x[i]}, {y[i]})");
    }

    [TestMethod]
    public void CubicSpline_Natural_EndpointSecondDerivative_ShouldBeZero()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 1, 2, 1, 3, 2 };
        var spline = new CubicSplineInterpolation(x, y, SplineBoundary.Natural);

        Assert.AreEqual(0, spline.SecondDerivative(x[0]), 1e-8, "S''(x0) should be 0");
        Assert.AreEqual(0, spline.SecondDerivative(x[4]), 1e-8, "S''(xn) should be 0");
    }

    [TestMethod]
    public void CubicSpline_Clamped_ShouldMatchEndpointDerivatives()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 0, 1, 0, 1, 0 };
        double leftDeriv = 2.0, rightDeriv = -1.0;

        var spline = new CubicSplineInterpolation(x, y, SplineBoundary.Clamped, leftDeriv, rightDeriv);

        Assert.AreEqual(leftDeriv, spline.Derivative(x[0]), 1e-8, "Left derivative should match");
        Assert.AreEqual(rightDeriv, spline.Derivative(x[4]), 1e-8, "Right derivative should match");
    }

    [TestMethod]
    public void CubicSpline_NotAKnot_ShouldPassThroughPoints()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 1, 2, 0, 3, 1 };
        var spline = new CubicSplineInterpolation(x, y, SplineBoundary.NotAKnot);

        for (int i = 0; i < x.Length; i++)
            Assert.AreEqual(y[i], spline.Evaluate(x[i]), 1e-10);
    }

    [TestMethod]
    public void CubicSpline_CubicPolynomial_NotAKnot_ShouldBeExact()
    {
        // Not-a-knot cubic spline should exactly reproduce any cubic polynomial
        Func<double, double> f = x => x * x * x - 2 * x * x + x - 3;
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = x.Select(f).ToArray();
        var spline = new CubicSplineInterpolation(x, y, SplineBoundary.NotAKnot);

        Assert.AreEqual(f(0.5), spline.Evaluate(0.5), 1e-6);
        Assert.AreEqual(f(2.5), spline.Evaluate(2.5), 1e-6);
        Assert.AreEqual(f(3.7), spline.Evaluate(3.7), 1e-6);
    }

    [TestMethod]
    public void CubicSpline_Derivative_ShouldApproximateAnalytic()
    {
        // f(x) = sin(x), f'(x) = cos(x)
        int n = 10;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i * Math.PI / (n - 1);
            y[i] = Math.Sin(x[i]);
        }

        var spline = new CubicSplineInterpolation(x, y);
        double mid = Math.PI / 4;
        double deriv = spline.Derivative(mid);

        Assert.AreEqual(Math.Cos(mid), deriv, 0.01, "Derivative should approximate cos(x)");
    }

    [TestMethod]
    public void CubicSpline_EvaluateArray_ShouldMatchSingleEvaluations()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 1, 0, 1, 0, 1 };
        var spline = new CubicSplineInterpolation(x, y);

        double[] xs = { 0.5, 1.5, 2.5, 3.5 };
        double[] results = spline.Evaluate(xs);

        for (int i = 0; i < xs.Length; i++)
            Assert.AreEqual(spline.Evaluate(xs[i]), results[i], 1e-15);
    }

    [TestMethod]
    public void CubicSpline_ExtensionMethod_ShouldWork()
    {
        var serie = new List<Serie>
        {
            new() { Index = 0, Value = 0 },
            new() { Index = 1, Value = 1 },
            new() { Index = 2, Value = 0 },
            new() { Index = 3, Value = 1 },
            new() { Index = 4, Value = 0 }
        };

        double val = serie.Interpolate(p => (p.Index, p.Value), 1.5, InterpolationType.CubicSpline);

        var spline = new CubicSplineInterpolation(
            serie.Select(s => s.Index).ToArray(),
            serie.Select(s => s.Value).ToArray());

        Assert.AreEqual(spline.Evaluate(1.5), val, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Rational Interpolation
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Rational_BulirschStoer_ShouldPassThroughPoints()
    {
        double[] x = { 0, 1, 2, 3 };
        double[] y = { 1, 0.5, 0.333, 0.25 };
        var rat = new RationalInterpolation(x, y);

        for (int i = 0; i < x.Length; i++)
        {
            var (val, _) = rat.Evaluate(x[i]);
            Assert.AreEqual(y[i], val, 1e-8, $"Should pass through ({x[i]}, {y[i]})");
        }
    }

    [TestMethod]
    public void Rational_BulirschStoer_1OverX_ShouldBeAccurate()
    {
        // f(x) = 1/(1+x) — rational functions excel at this
        Func<double, double> f = x => 1.0 / (1.0 + x);
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = x.Select(f).ToArray();
        var rat = new RationalInterpolation(x, y);

        var (val, err) = rat.Evaluate(2.5);
        Assert.AreEqual(f(2.5), val, 1e-6, "Should approximate 1/(1+x) well");
    }

    [TestMethod]
    public void Rational_FloaterHormann_ShouldApproximate()
    {
        Func<double, double> f = x => 1.0 / (1.0 + x);
        double[] x = { 0, 1, 2, 3, 4, 5 };
        double[] y = x.Select(f).ToArray();
        var rat = new RationalInterpolation(x, y);

        double val = rat.EvaluateFloaterHormann(2.5, d: 3);
        Assert.AreEqual(f(2.5), val, 0.05, "Floater-Hormann should be close");
    }

    [TestMethod]
    public void Rational_FloaterHormann_ExactAtDataPoints()
    {
        double[] x = { 0, 1, 2, 3, 4 };
        double[] y = { 2, 5, 1, 3, 7 };
        var rat = new RationalInterpolation(x, y);

        for (int i = 0; i < x.Length; i++)
        {
            double val = rat.EvaluateFloaterHormann(x[i], d: 2);
            Assert.AreEqual(y[i], val, 1e-10, $"Should be exact at x={x[i]}");
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Trigonometric Interpolation
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Trigonometric_PureSine_ShouldInterpolateExactly()
    {
        // f(x) = sin(x) on [0, 2π) with equispaced points
        int N = 9; // odd → exact interpolation
        double period = 2 * Math.PI;
        double[] x = new double[N];
        double[] y = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = period * i / N;
            y[i] = Math.Sin(x[i]);
        }

        var trig = new TrigonometricInterpolation(x, y, period);

        // Should be very close at a non-data point
        double xi = Math.PI / 5;
        Assert.AreEqual(Math.Sin(xi), trig.Evaluate(xi), 1e-6,
            "Should closely approximate sin(x)");
    }

    [TestMethod]
    public void Trigonometric_Derivative_ShouldApproximateCosine()
    {
        int N = 17;
        double period = 2 * Math.PI;
        double[] x = new double[N];
        double[] y = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = period * i / N;
            y[i] = Math.Sin(x[i]);
        }

        var trig = new TrigonometricInterpolation(x, y, period);
        double xi = 1.0;
        Assert.AreEqual(Math.Cos(xi), trig.Derivative(xi), 0.01,
            "Derivative of sin(x) should ≈ cos(x)");
    }

    [TestMethod]
    public void Trigonometric_CosineCoefficients_ShouldReflectSignal()
    {
        // f(x) = 3 + 2*cos(x) → a_0 ≈ 6 (since a_0/2 = 3), a_1 ≈ 2
        int N = 17;
        double period = 2 * Math.PI;
        double[] x = new double[N];
        double[] y = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = period * i / N;
            y[i] = 3 + 2 * Math.Cos(x[i]);
        }

        var trig = new TrigonometricInterpolation(x, y, period);
        double[] a = trig.CosineCoefficients;

        Assert.AreEqual(6.0, a[0], 0.1, "a_0 should be ~6 (since mean = 3 = a_0/2)");
        Assert.AreEqual(2.0, a[1], 0.1, "a_1 should be ~2");
    }

    [TestMethod]
    public void Trigonometric_EvaluateArray_ShouldMatchSingle()
    {
        int N = 8;
        double period = 2 * Math.PI;
        double[] x = new double[N], y = new double[N];
        for (int i = 0; i < N; i++)
        {
            x[i] = period * i / N;
            y[i] = Math.Sin(x[i]) + 0.5 * Math.Cos(2 * x[i]);
        }

        var trig = new TrigonometricInterpolation(x, y, period);
        double[] xs = { 0.3, 1.0, 2.5, 4.0 };
        double[] results = trig.Evaluate(xs);

        for (int i = 0; i < xs.Length; i++)
            Assert.AreEqual(trig.Evaluate(xs[i]), results[i], 1e-15);
    }

    [TestMethod]
    public void Trigonometric_ExtensionMethod_ShouldWork()
    {
        var serie = new List<Serie>();
        int N = 9;
        double period = 2 * Math.PI;
        for (int i = 0; i < N; i++)
        {
            double xi = period * i / N;
            serie.Add(new Serie { Index = xi, Value = Math.Sin(xi) });
        }

        double val = serie.Interpolate(p => (p.Index, p.Value), 1.0, InterpolationType.Trigonometric);

        var trig = new TrigonometricInterpolation(
            serie.Select(s => s.Index).ToArray(),
            serie.Select(s => s.Value).ToArray());

        Assert.AreEqual(trig.Evaluate(1.0), val, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multivariate Interpolation — IDW
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void IDW_ShouldReturnExactValueAtDataPoint()
    {
        double[][] points = { new[] { 0.0, 0.0 }, new[] { 1.0, 0.0 }, new[] { 0.0, 1.0 } };
        double[] values = { 10, 20, 30 };
        var interp = new MultivariateInterpolation(points, values);

        Assert.AreEqual(10.0, interp.EvaluateIDW(new[] { 0.0, 0.0 }), 1e-10);
        Assert.AreEqual(20.0, interp.EvaluateIDW(new[] { 1.0, 0.0 }), 1e-10);
    }

    [TestMethod]
    public void IDW_MidpointOfTwoEqualWeight_ShouldBeAverage()
    {
        // Two points equidistant from query → should return average
        double[][] points = { new[] { 0.0 }, new[] { 2.0 } };
        double[] values = { 10, 20 };
        var interp = new MultivariateInterpolation(points, values);

        double val = interp.EvaluateIDW(new[] { 1.0 }, power: 2);
        Assert.AreEqual(15.0, val, 1e-10, "Midpoint should be average of equal-distance values");
    }

    [TestMethod]
    public void IDW_CloserPointShouldDominate()
    {
        double[][] points = { new[] { 0.0 }, new[] { 10.0 } };
        double[] values = { 100, 0 };
        var interp = new MultivariateInterpolation(points, values);

        double val = interp.EvaluateIDW(new[] { 1.0 }, power: 2);
        Assert.IsTrue(val > 50, $"Value near first point should be > 50, got {val:F2}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multivariate Interpolation — RBF
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void RBF_Gaussian_ShouldPassThroughDataPoints()
    {
        double[][] points = {
            new[] { 0.0, 0.0 }, new[] { 1.0, 0.0 },
            new[] { 0.0, 1.0 }, new[] { 1.0, 1.0 }
        };
        double[] values = { 1, 2, 3, 4 };
        var interp = new MultivariateInterpolation(points, values);

        for (int i = 0; i < points.Length; i++)
        {
            double val = interp.EvaluateRBF(points[i], RbfKernel.Gaussian);
            Assert.AreEqual(values[i], val, 1e-6,
                $"RBF should pass through point {i}");
        }
    }

    [TestMethod]
    public void RBF_Multiquadric_ShouldPassThroughDataPoints()
    {
        double[][] points = {
            new[] { 0.0, 0.0 }, new[] { 1.0, 0.0 },
            new[] { 0.0, 1.0 }, new[] { 1.0, 1.0 }
        };
        double[] values = { 5, 10, 15, 20 };
        var interp = new MultivariateInterpolation(points, values);

        for (int i = 0; i < points.Length; i++)
        {
            double val = interp.EvaluateRBF(points[i], RbfKernel.Multiquadric);
            Assert.AreEqual(values[i], val, 1e-5);
        }
    }

    [TestMethod]
    public void RBF_QuadraticFunction_ShouldApproximate()
    {
        // f(x,y) = x² + y² — sample on a 3×3 grid
        Func<double, double, double> f = (x, y) => x * x + y * y;
        var points = new List<double[]>();
        var values = new List<double>();
        for (int i = 0; i <= 2; i++)
            for (int j = 0; j <= 2; j++)
            {
                points.Add(new double[] { i, j });
                values.Add(f(i, j));
            }

        var interp = new MultivariateInterpolation(points.ToArray(), values.ToArray());
        double val = interp.EvaluateRBF(new[] { 0.5, 0.5 }, RbfKernel.Multiquadric);

        Assert.AreEqual(f(0.5, 0.5), val, 0.5,
            $"RBF should approximate x²+y² at (0.5,0.5), got {val:F4}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multivariate Interpolation — Bilinear
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Bilinear_ShouldInterpolateLinearFunction()
    {
        // f(x,y) = 2x + 3y — bilinear should be exact for linear functions
        double[] xGrid = { 0, 1, 2 };
        double[] yGrid = { 0, 1, 2 };
        var grid = new double[3, 3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                grid[i, j] = 2 * xGrid[i] + 3 * yGrid[j];

        double val = MultivariateInterpolation.Bilinear(xGrid, yGrid, grid, 0.5, 0.5);
        double expected = 2 * 0.5 + 3 * 0.5;
        Assert.AreEqual(expected, val, 1e-10);
    }

    [TestMethod]
    public void Bilinear_CornerValues_ShouldBeExact()
    {
        double[] xGrid = { 0, 1 };
        double[] yGrid = { 0, 1 };
        var grid = new double[2, 2] { { 10, 20 }, { 30, 40 } };

        Assert.AreEqual(10, MultivariateInterpolation.Bilinear(xGrid, yGrid, grid, 0, 0), 1e-10);
        Assert.AreEqual(20, MultivariateInterpolation.Bilinear(xGrid, yGrid, grid, 0, 1), 1e-10);
        Assert.AreEqual(30, MultivariateInterpolation.Bilinear(xGrid, yGrid, grid, 1, 0), 1e-10);
        Assert.AreEqual(40, MultivariateInterpolation.Bilinear(xGrid, yGrid, grid, 1, 1), 1e-10);
    }

    [TestMethod]
    public void Bilinear_Center_ShouldBeAverageOfCorners()
    {
        double[] xGrid = { 0, 1 };
        double[] yGrid = { 0, 1 };
        var grid = new double[2, 2] { { 0, 0 }, { 0, 4 } };

        double val = MultivariateInterpolation.Bilinear(xGrid, yGrid, grid, 0.5, 0.5);
        Assert.AreEqual(1.0, val, 1e-10, "Center of [0,0,0,4] should be 1.0");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multivariate Interpolation — Trilinear
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Trilinear_LinearFunction_ShouldBeExact()
    {
        // f(x,y,z) = x + 2y + 3z
        double[] xg = { 0, 1 };
        double[] yg = { 0, 1 };
        double[] zg = { 0, 1 };
        var grid = new double[2, 2, 2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    grid[i, j, k] = xg[i] + 2 * yg[j] + 3 * zg[k];

        double val = MultivariateInterpolation.Trilinear(xg, yg, zg, grid, 0.3, 0.6, 0.8);
        double expected = 0.3 + 2 * 0.6 + 3 * 0.8;
        Assert.AreEqual(expected, val, 1e-10);
    }

    [TestMethod]
    public void Trilinear_CornerValues_ShouldBeExact()
    {
        double[] xg = { 0, 1 };
        double[] yg = { 0, 1 };
        double[] zg = { 0, 1 };
        var grid = new double[2, 2, 2];
        grid[0, 0, 0] = 1;
        grid[1, 1, 1] = 8;

        Assert.AreEqual(1, MultivariateInterpolation.Trilinear(xg, yg, zg, grid, 0, 0, 0), 1e-10);
        Assert.AreEqual(8, MultivariateInterpolation.Trilinear(xg, yg, zg, grid, 1, 1, 1), 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Edge cases & validation
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Polynomial_TooFewPoints_ShouldThrow()
    {
        new PolynomialInterpolation(new[] { 1.0 }, new[] { 2.0 });
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void CubicSpline_TooFewPoints_ShouldThrow()
    {
        new CubicSplineInterpolation(new[] { 1.0, 2.0 }, new[] { 1.0, 2.0 });
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Polynomial_MismatchedLengths_ShouldThrow()
    {
        new PolynomialInterpolation(new[] { 1.0, 2.0 }, new[] { 1.0 });
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void Multivariate_MismatchedLengths_ShouldThrow()
    {
        new MultivariateInterpolation(
            new[] { new[] { 0.0 }, new[] { 1.0 } },
            new[] { 1.0 });
    }
}
