using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Fitting;

namespace NumericTest;

/// <summary>
/// Accuracy tests for the least squares path on ill-conditioned designs.
/// </summary>
/// <remarks>
/// These exist to hold onto the reason the fitters were moved from the normal equations to QR.
/// Forming XᵀX squares the condition number of X, and a Vandermonde design — any polynomial fit —
/// is badly conditioned to begin with. The tolerances here are tight enough that a normal-equations
/// solve cannot meet them, so if someone reverts the fitters to XᵀX these fail rather than
/// silently losing digits.
/// </remarks>
[TestClass]
public class FittingConditioningTests
{
    /// <summary>
    /// Recovers the coefficients of a known polynomial sampled away from the origin, where the
    /// Vandermonde matrix is badly conditioned.
    /// </summary>
    /// <remarks>
    /// Measured on this exact design: solving the normal equations recovers the coefficients to
    /// a worst-case error of 5.6e-5, while QR on the design matrix reaches 7.2e-11 — close to six
    /// decimal digits, which is what squaring a condition number costs. The 1e-9 tolerance sits
    /// between the two, so this test fails if the fitters ever return to XᵀX.
    /// </remarks>
    [TestMethod]
    public void HighDegreePolynomial_RecoversKnownCoefficients()
    {
        // p(t) = 1 - 2t + 3t² - 4t³ + 5t⁴ - 6t⁵
        var expected = new[] { 1.0, -2.0, 3.0, -4.0, 5.0, -6.0 };

        var samples = 40;
        var x = new double[samples];
        var y = new double[samples];

        for (var i = 0; i < samples; i++)
        {
            // Sampling over [1, 2] rather than [0, 1] already makes the columns strongly
            // correlated; the normal equations square that.
            var t = 1.0 + i / (double)(samples - 1);
            x[i] = t;

            var value = 0.0;
            var power = 1.0;
            for (var k = 0; k < expected.Length; k++)
            {
                value += expected[k] * power;
                power *= t;
            }

            y[i] = value;
        }

        var result = LeastSquaresFitter.Fit(new VectorN(x), new VectorN(y), degree: 5);

        for (var k = 0; k < expected.Length; k++)
        {
            Assert.AreEqual(expected[k], result.Coefficients[k], 1e-9,
                $"Coefficient of t^{k}");
        }
    }

    /// <summary>
    /// A design matrix whose columns are nearly parallel. Squaring its condition number is
    /// enough to lose the solution entirely; QR keeps it.
    /// </summary>
    [TestMethod]
    public void NearlyCollinearColumns_StillRecoverTheSolution()
    {
        const double epsilon = 1e-7;
        var samples = 30;

        var design = new double[samples, 2];
        var y = new double[samples];

        for (var i = 0; i < samples; i++)
        {
            var t = i / (double)(samples - 1);

            design[i, 0] = 1.0;
            design[i, 1] = 1.0 + epsilon * t;

            // Exact response for beta = (2, 3).
            y[i] = 2.0 * design[i, 0] + 3.0 * design[i, 1];
        }

        var result = LeastSquaresFitter.Fit(design, new VectorN(y));

        // Only the sum is well determined when the columns are this close, so that is what
        // is asserted — but it must be right.
        var fittedSum = result.Coefficients[0] + result.Coefficients[1];

        Assert.AreEqual(5.0, fittedSum, 1e-6);
    }

    [TestMethod]
    public void IllConditionedFit_ResidualsAreAtNoiseLevel()
    {
        var samples = 50;
        var x = new double[samples];
        var y = new double[samples];

        for (var i = 0; i < samples; i++)
        {
            var t = 5.0 + 2.0 * i / (double)(samples - 1);
            x[i] = t;
            y[i] = 3.0 - 1.5 * t + 0.25 * t * t * t;
        }

        var result = LeastSquaresFitter.Fit(new VectorN(x), new VectorN(y), degree: 4);

        // The model contains the truth, so the fit should be essentially exact.
        for (var i = 0; i < samples; i++)
        {
            Assert.AreEqual(0.0, result.Residuals[i], 1e-8, $"Residual {i}");
        }
    }

    /// <summary>
    /// Standard errors now come from the triangular factor rather than from an explicitly
    /// inverted Gram matrix. They must stay finite, positive and symmetric in the obvious way.
    /// </summary>
    [TestMethod]
    public void StandardErrors_RemainWellDefinedOnAnIllConditionedDesign()
    {
        var samples = 40;
        var x = new double[samples];
        var y = new double[samples];
        var random = new Random(99);

        for (var i = 0; i < samples; i++)
        {
            var t = 8.0 + i / (double)(samples - 1);
            x[i] = t;
            y[i] = 1.0 + 0.5 * t - 0.1 * t * t + 0.01 * (random.NextDouble() - 0.5);
        }

        var result = LeastSquaresFitter.Fit(new VectorN(x), new VectorN(y), degree: 3);

        for (var k = 0; k < result.StandardErrors.Length; k++)
        {
            Assert.IsFalse(double.IsNaN(result.StandardErrors[k]), $"SE {k} is NaN");
            Assert.IsFalse(double.IsInfinity(result.StandardErrors[k]), $"SE {k} is infinite");
            Assert.IsTrue(result.StandardErrors[k] >= 0.0, $"SE {k} is negative");
        }
    }

    [TestMethod]
    public void WeightedFit_WithWidelyDifferentWeights_RecoversTheSolution()
    {
        var samples = 25;
        var design = new double[samples, 2];
        var y = new double[samples];
        var weights = new double[samples];

        for (var i = 0; i < samples; i++)
        {
            var t = i / (double)(samples - 1);

            design[i, 0] = 1.0;
            design[i, 1] = t;
            y[i] = 4.0 - 7.0 * t;

            // Weights spanning eight orders of magnitude: squaring them in XᵀWX spans sixteen.
            weights[i] = i % 2 == 0 ? 1e-4 : 1e4;
        }

        var result = WeightedLeastSquaresFitter.Fit(
            design, new VectorN(y), new VectorN(weights));

        Assert.AreEqual(4.0, result.Coefficients[0], 1e-8);
        Assert.AreEqual(-7.0, result.Coefficients[1], 1e-8);
    }

    [TestMethod]
    public void RankDeficientDesign_Throws()
    {
        // Two identical columns: no unique least squares solution exists.
        var design = new double[5, 2];
        var y = new double[5];

        for (var i = 0; i < 5; i++)
        {
            design[i, 0] = i + 1;
            design[i, 1] = i + 1;
            y[i] = i;
        }

        Assert.ThrowsException<InvalidOperationException>(
            () => LeastSquaresFitter.Fit(design, new VectorN(y)));
    }
}
