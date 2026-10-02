using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.RootFinding;

namespace NumericTest;

[TestClass]
public class RootFindingTests
{
    // f(x) = x² − 4, root at x = 2.
    private static readonly Func<double, double> Quadratic = x => x * x - 4.0;

    // f(x) = cos(x) − x, root at the Dottie number.
    private static readonly Func<double, double> CosineMinusX = x => Math.Cos(x) - x;
    private const double DottieNumber = 0.7390851332151607;

    // f(x) = x³ − 2x − 5, the classic test case from Numerical Recipes.
    private static readonly Func<double, double> Cubic = x => x * x * x - 2.0 * x - 5.0;
    private const double CubicRoot = 2.0945514815423265;

    #region Bisection

    [TestMethod]
    public void Bisection_Quadratic_FindsRoot()
    {
        var result = RootFinder.Bisection(Quadratic, 0.0, 5.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(2.0, result.Value, 1e-9);
    }

    [TestMethod]
    public void Bisection_ReversedBracket_IsAccepted()
    {
        var result = RootFinder.Bisection(Quadratic, 5.0, 0.0);

        Assert.IsTrue(result.Converged);
        Assert.AreEqual(2.0, result.Value, 1e-9);
    }

    [TestMethod]
    public void Bisection_RootExactlyOnEndpoint_ReturnsImmediately()
    {
        var result = RootFinder.Bisection(Quadratic, 2.0, 5.0);

        Assert.IsTrue(result.Converged);
        Assert.AreEqual(2.0, result.Value, 0.0);
        Assert.AreEqual(0, result.Iterations);
    }

    [TestMethod]
    public void Bisection_NoSignChange_Throws()
    {
        Assert.ThrowsException<ArgumentException>(
            () => RootFinder.Bisection(Quadratic, 3.0, 5.0));
    }

    [TestMethod]
    public void Bisection_NonSmoothFunction_StillConverges()
    {
        // |x| − 1 has a corner at the origin but a clean sign change at x = 1.
        Func<double, double> kinked = x => Math.Abs(x) - 1.0;

        var result = RootFinder.Bisection(kinked, 0.5, 3.0);

        Assert.IsTrue(result.Converged);
        Assert.AreEqual(1.0, result.Value, 1e-9);
    }

    #endregion

    #region Secant

    [TestMethod]
    public void Secant_Quadratic_FindsRoot()
    {
        var result = RootFinder.Secant(Quadratic, 0.0, 5.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(2.0, result.Value, 1e-9);
    }

    [TestMethod]
    public void Secant_EqualStartingPoints_Throws()
    {
        Assert.ThrowsException<ArgumentException>(
            () => RootFinder.Secant(Quadratic, 1.0, 1.0));
    }

    [TestMethod]
    public void Secant_HorizontalSecant_ReportsFailureInsteadOfDividingByZero()
    {
        // A constant non-zero function gives f(x₁) − f(x₀) = 0: no step is defined.
        Func<double, double> constant = _ => 3.0;

        var result = RootFinder.Secant(constant, 0.0, 1.0);

        Assert.IsFalse(result.Converged);
        Assert.IsFalse(double.IsNaN(result.Value));
    }

    #endregion

    #region Brent

    [TestMethod]
    public void Brent_Quadratic_FindsRoot()
    {
        var result = RootFinder.Brent(Quadratic, 0.0, 5.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(2.0, result.Value, 1e-12);
    }

    [TestMethod]
    public void Brent_Transcendental_MatchesKnownRoot()
    {
        var result = RootFinder.Brent(CosineMinusX, 0.0, 1.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(DottieNumber, result.Value, 1e-12);
    }

    [TestMethod]
    public void Brent_Cubic_MatchesKnownRoot()
    {
        var result = RootFinder.Brent(Cubic, 2.0, 3.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(CubicRoot, result.Value, 1e-12);
    }

    [TestMethod]
    public void Brent_TripleRoot_ConvergesDespiteVanishingDerivative()
    {
        // (x − 1)³ is flat at its root, which defeats plain interpolation.
        Func<double, double> tripleRoot = x => Math.Pow(x - 1.0, 3);

        var result = RootFinder.Brent(tripleRoot, -2.0, 4.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(1.0, result.Value, 1e-4);
    }

    [TestMethod]
    public void Brent_CubeRoot_ConvergesWhereNewtonDiverges()
    {
        // ∛x has an infinite derivative at the root — Newton's method never converges on it.
        Func<double, double> cubeRoot = Math.Cbrt;

        var brent = RootFinder.Brent(cubeRoot, -1.0, 2.0);
        var newton = RootFinder.Newton(cubeRoot, 1.0);

        Assert.IsTrue(brent.Converged, brent.ToString());
        Assert.AreEqual(0.0, brent.Value, 1e-9);
        Assert.IsFalse(newton.Converged, "Newton is expected to fail on the cube root.");
    }

    [TestMethod]
    public void Brent_NoSignChange_Throws()
    {
        Assert.ThrowsException<ArgumentException>(
            () => RootFinder.Brent(Quadratic, 3.0, 5.0));
    }

    [TestMethod]
    public void Brent_ConvergesInFarFewerIterationsThanBisection()
    {
        var brent = RootFinder.Brent(CosineMinusX, 0.0, 1.0);
        var bisection = RootFinder.Bisection(CosineMinusX, 0.0, 1.0);

        Assert.IsTrue(brent.Converged);
        Assert.IsTrue(bisection.Converged);
        Assert.IsTrue(
            brent.Iterations < bisection.Iterations / 2,
            $"Brent took {brent.Iterations} iterations, bisection took {bisection.Iterations}.");
    }

    #endregion

    #region Newton

    [TestMethod]
    public void Newton_Quadratic_FindsRoot()
    {
        var result = RootFinder.Newton(Quadratic, 1.0);

        Assert.IsTrue(result.Converged, result.ToString());
        Assert.AreEqual(2.0, result.Value, 1e-6);
    }

    [TestMethod]
    public void Newton_AnalyticDerivative_IsMoreAccurateThanFiniteDifference()
    {
        var analytic = RootFinder.Newton(Quadratic, x => 2.0 * x, 1.0);
        var finiteDifference = RootFinder.Newton(Quadratic, 1.0);

        Assert.IsTrue(analytic.Converged, analytic.ToString());
        Assert.IsTrue(finiteDifference.Converged, finiteDifference.ToString());
        Assert.AreEqual(2.0, analytic.Value, 1e-12);
        Assert.IsTrue(
            analytic.Residual <= finiteDifference.Residual,
            $"Analytic residual {analytic.Residual:G6} should not exceed " +
            $"finite-difference residual {finiteDifference.Residual:G6}.");
    }

    [TestMethod]
    public void Newton_StopsEarlyInsteadOfRunningTheFullBudget()
    {
        var result = RootFinder.Newton(Quadratic, x => 2.0 * x, 1.0, maxIterations: 100);

        Assert.IsTrue(result.Converged);
        Assert.IsTrue(result.Iterations < 20, $"Took {result.Iterations} iterations.");
    }

    [TestMethod]
    public void Newton_VanishingDerivative_ReportsFailureInsteadOfNaN()
    {
        // x³ − 3x + 3 has f'(1) = 0 exactly, so no Newton step exists from x₀ = 1.
        Func<double, double> function = x => x * x * x - 3.0 * x + 3.0;
        Func<double, double> derivative = x => 3.0 * x * x - 3.0;

        var result = RootFinder.Newton(function, derivative, 1.0);

        Assert.IsFalse(result.Converged);
        Assert.IsFalse(double.IsNaN(result.Value));
        Assert.IsFalse(double.IsInfinity(result.Value));
    }

    [TestMethod]
    public void Newton_DivergentStart_ReportsFailure()
    {
        // arctan diverges under Newton's method for any |x₀| greater than ≈1.3917.
        var result = RootFinder.Newton(Math.Atan, x => 1.0 / (1.0 + x * x), 5.0);

        Assert.IsFalse(result.Converged);
    }

    [TestMethod]
    public void Newton_NoRealRoot_ReportsFailure()
    {
        Func<double, double> noRoot = x => x * x + 1.0;

        var result = RootFinder.Newton(noRoot, x => 2.0 * x, 3.0);

        Assert.IsFalse(result.Converged);
    }

    [TestMethod]
    public void Newton_NullDerivative_Throws()
    {
        Assert.ThrowsException<ArgumentNullException>(
            () => RootFinder.Newton(Quadratic, null!, 1.0));
    }

    #endregion

    #region Argument validation

    [TestMethod]
    public void NonPositiveTolerance_Throws()
    {
        Assert.ThrowsException<ArgumentOutOfRangeException>(
            () => RootFinder.Brent(Quadratic, 0.0, 5.0, tolerance: 0.0));
    }

    [TestMethod]
    public void ZeroIterationBudget_Throws()
    {
        Assert.ThrowsException<ArgumentOutOfRangeException>(
            () => RootFinder.Brent(Quadratic, 0.0, 5.0, maxIterations: 0));
    }

    [TestMethod]
    public void DegenerateBracket_Throws()
    {
        Assert.ThrowsException<ArgumentException>(
            () => RootFinder.Brent(Quadratic, 2.5, 2.5));
    }

    [TestMethod]
    public void NullFunction_Throws()
    {
        Assert.ThrowsException<ArgumentNullException>(
            () => RootFinder.Brent(null!, 0.0, 5.0));
    }

    #endregion

    #region RootResult

    [TestMethod]
    public void EnsureConverged_OnSuccess_ReturnsValue()
    {
        var result = RootFinder.Brent(Quadratic, 0.0, 5.0);

        Assert.AreEqual(2.0, result.EnsureConverged(), 1e-12);
    }

    [TestMethod]
    public void EnsureConverged_OnFailure_Throws()
    {
        var result = RootFinder.Newton(Math.Atan, x => 1.0 / (1.0 + x * x), 5.0);

        Assert.ThrowsException<InvalidOperationException>(() => result.EnsureConverged());
    }

    [TestMethod]
    public void Result_ReportsResidualAtReturnedPoint()
    {
        var result = RootFinder.Brent(Cubic, 2.0, 3.0);

        Assert.AreEqual(Math.Abs(Cubic(result.Value)), result.Residual, 1e-15);
    }

    #endregion

    #region FindRoot facade and NewtonRaphson compatibility

    [TestMethod]
    public void FindRoot_DelegatesToBrent()
    {
        var result = CosineMinusX.FindRoot(0.0, 1.0);

        Assert.IsTrue(result.Converged);
        Assert.AreEqual(DottieNumber, result.Value, 1e-12);
        Assert.AreEqual("Brent", result.Method);
    }

    [TestMethod]
    public void NewtonRaphson_KeepsWorkingThroughTheNewImplementation()
    {
        var root = Quadratic.NewtonRaphson();

        Assert.AreEqual(2.0, root, 1e-6);
    }

    [TestMethod]
    public void NewtonRaphson_RespectsInitialGuess()
    {
        var negativeRoot = Quadratic.NewtonRaphson(-1.0);

        Assert.AreEqual(-2.0, negativeRoot, 1e-6);
    }

    #endregion
}
