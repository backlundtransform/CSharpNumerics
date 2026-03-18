using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
using CSharpNumerics.Numerics.Optimization.Strategies;

namespace NumericTest;

[TestClass]
public class OptimizationTests
{
    // ── f(x) = x² ──────────────────────────────────────────────
    // Minimum at x = 0. Gradient = 2x.

    private class Quadratic1D : IObjectiveFunction
    {
        public int Dimension => 1;
        public double Evaluate(double[] x) => x[0] * x[0];
        public double[] Gradient(double[] x) => new[] { 2 * x[0] };
    }

    // ── f(x,y) = (x-3)² + (y+1)² ──────────────────────────────
    // Minimum at (3, -1).

    private class Quadratic2D : IObjectiveFunction
    {
        public int Dimension => 2;
        public double Evaluate(double[] x) => (x[0] - 3) * (x[0] - 3) + (x[1] + 1) * (x[1] + 1);
        public double[] Gradient(double[] x) => new[] { 2 * (x[0] - 3), 2 * (x[1] + 1) };
    }

    // ── Rosenbrock f(x,y) = (1-x)² + 100(y-x²)² ──────────────
    // Minimum at (1, 1).

    private class Rosenbrock : IObjectiveFunction
    {
        public int Dimension => 2;
        public double Evaluate(double[] x)
            => (1 - x[0]) * (1 - x[0]) + 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
        public double[] Gradient(double[] x)
        {
            double dx = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
            double dy = 200 * (x[1] - x[0] * x[0]);
            return new[] { dx, dy };
        }
    }

    // ═══════════════════════════════════════════════════════════
    // GradientDescent tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void GradientDescent_Minimises_Quadratic1D()
    {
        var gd = new GradientDescent(learningRate: 0.1);
        var conv = new MaxIterationsOrTolerance(maxIterations: 200, tolerance: 1e-10);
        var minimizer = new Minimizer(gd, conv);

        var x0 = new VectorN(new[] { 5.0 });
        var result = minimizer.Minimize(new Quadratic1D(), x0);

        Assert.AreEqual(0.0, result[0], 1e-6);
    }

    [TestMethod]
    public void GradientDescent_Minimises_Quadratic2D()
    {
        var gd = new GradientDescent(learningRate: 0.1);
        var conv = new MaxIterationsOrTolerance(maxIterations: 500, tolerance: 1e-10);
        var minimizer = new Minimizer(gd, conv);

        var x0 = new VectorN(new[] { 0.0, 0.0 });
        var result = minimizer.Minimize(new Quadratic2D(), x0);

        Assert.AreEqual(3.0, result[0], 1e-5);
        Assert.AreEqual(-1.0, result[1], 1e-5);
    }

    [TestMethod]
    public void GradientDescent_WithMomentum_ConvergesFaster()
    {
        var gdPlain = new GradientDescent(learningRate: 0.01);
        var gdMom = new GradientDescent(learningRate: 0.01, momentum: 0.9);
        var conv1 = new MaxIterationsOrTolerance(maxIterations: 100, tolerance: 1e-10);
        var conv2 = new MaxIterationsOrTolerance(maxIterations: 100, tolerance: 1e-10);

        var f = new Quadratic2D();
        var x0 = new VectorN(new[] { 10.0, -10.0 });

        var m1 = new Minimizer(gdPlain, conv1);
        var m2 = new Minimizer(gdMom, conv2);

        var r1 = m1.Minimize(f, x0);
        var r2 = m2.Minimize(f, x0);

        double loss1 = f.Evaluate(r1.Values);
        double loss2 = f.Evaluate(r2.Values);

        Assert.IsTrue(loss2 < loss1, $"Momentum ({loss2:E3}) should beat plain GD ({loss1:E3}) in same iterations");
    }

    [TestMethod]
    public void GradientDescent_WithL2_RegularisesTowardsZero()
    {
        var gd = new GradientDescent(learningRate: 0.1, l2: 0.5);
        var x = new VectorN(new[] { 3.0, -1.0 });
        // gradient is zero but L2 pulls towards 0
        var zeroGrad = new VectorN(new[] { 0.0, 0.0 });
        var updated = gd.Step(x, zeroGrad);

        Assert.IsTrue(Math.Abs(updated[0]) < Math.Abs(x[0]));
        Assert.IsTrue(Math.Abs(updated[1]) < Math.Abs(x[1]));
    }

    [TestMethod]
    public void GradientDescent_Reset_ClearsState()
    {
        var gd = new GradientDescent(learningRate: 0.1, momentum: 0.9);
        var x = new VectorN(new[] { 5.0 });
        var g = new VectorN(new[] { 2.0 });

        gd.Step(x, g); // builds velocity
        gd.Reset();

        // After reset, first step should behave like fresh optimizer
        var result = gd.Step(x, g);
        var expected = x[0] - 0.1 * 2.0; // lr * gradient (no momentum contribution)
        Assert.AreEqual(expected, result[0], 1e-12);
    }

    // ═══════════════════════════════════════════════════════════
    // Adam tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void Adam_Minimises_Quadratic2D()
    {
        var adam = new Adam(learningRate: 0.1);
        var conv = new MaxIterationsOrTolerance(maxIterations: 500, tolerance: 1e-10);
        var minimizer = new Minimizer(adam, conv);

        var x0 = new VectorN(new[] { 0.0, 0.0 });
        var result = minimizer.Minimize(new Quadratic2D(), x0);

        Assert.AreEqual(3.0, result[0], 1e-3);
        Assert.AreEqual(-1.0, result[1], 1e-3);
    }

    [TestMethod]
    public void Adam_Minimises_Rosenbrock()
    {
        var adam = new Adam(learningRate: 0.01);
        var conv = new MaxIterationsOrTolerance(maxIterations: 20000, tolerance: 1e-12);
        var minimizer = new Minimizer(adam, conv);

        var x0 = new VectorN(new[] { -1.0, 1.0 });
        var result = minimizer.Minimize(new Rosenbrock(), x0);

        Assert.AreEqual(1.0, result[0], 0.05);
        Assert.AreEqual(1.0, result[1], 0.05);
    }

    [TestMethod]
    public void Adam_Reset_ClearsState()
    {
        var adam = new Adam(learningRate: 0.01);
        var x = new VectorN(new[] { 5.0 });
        var g = new VectorN(new[] { 2.0 });

        var step1 = adam.Step(x, g);
        adam.Reset();
        var step2 = adam.Step(x, g);

        // After reset both should produce the same result
        Assert.AreEqual(step1[0], step2[0], 1e-12);
    }

    // ═══════════════════════════════════════════════════════════
    // CoordinateDescent tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void CoordinateDescent_SolvesLasso()
    {
        // y = 2*x1 + 3*x2 (no noise)
        var X = new double[,]
        {
            { 1, 0 },
            { 0, 1 },
            { 1, 1 },
            { 2, 1 }
        };
        var y = new double[] { 2, 3, 5, 7 };

        var cd = new CoordinateDescent(maxIterations: 1000);
        var w = cd.Solve(X, y, l1: 0.0, l2: 0.0);

        Assert.AreEqual(2.0, w[0], 0.1);
        Assert.AreEqual(3.0, w[1], 0.1);
    }

    [TestMethod]
    public void CoordinateDescent_L1_SparsifiesCoefficients()
    {
        // y = 3*x1 + 0*x2, with moderate L1 the irrelevant x2 should stay zero
        var X = new double[,]
        {
            { 1, 2 },
            { 2, -1 },
            { 3, 0 },
            { 4, 3 },
            { 5, -2 }
        };
        var y = new double[] { 3, 6, 9, 12, 15 }; // exactly 3*x1

        var cd = new CoordinateDescent(maxIterations: 1000);
        var w = cd.Solve(X, y, l1: 0.5, l2: 0.0);

        // x2 is uncorrelated with residual → should be zero under L1
        Assert.AreEqual(0.0, w[1], 0.01);
        Assert.IsTrue(Math.Abs(w[0]) > 1.0, "x1 coefficient should remain significant");
    }

    [TestMethod]
    public void CoordinateDescent_WithBiasSkip_DoesNotPenaliseBias()
    {
        var X = new double[,]
        {
            { 1, 1 },
            { 1, 2 },
            { 1, 3 },
            { 1, 4 }
        };
        var y = new double[] { 10, 12, 14, 16 }; // y = 8 + 2*x

        var cd = new CoordinateDescent(maxIterations: 1000);
        var w = cd.Solve(X, y, l1: 0.0, l2: 0.0, skipBiasRegularisation: true);

        Assert.AreEqual(8.0, w[0], 0.5);  // bias
        Assert.AreEqual(2.0, w[1], 0.5);  // slope
    }

    // ═══════════════════════════════════════════════════════════
    // EarlyStopping tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void EarlyStopping_TriggersAfterPatience()
    {
        var es = new EarlyStopping(patience: 3, minDelta: 0.01);

        Assert.IsFalse(es.Check(1.0));  // new best
        Assert.IsFalse(es.Check(0.5));  // new best
        Assert.IsFalse(es.Check(0.5));  // no improvement #1
        Assert.IsFalse(es.Check(0.5));  // no improvement #2
        Assert.IsTrue(es.Check(0.5));   // no improvement #3 → trigger
    }

    [TestMethod]
    public void EarlyStopping_ResetsOnImprovement()
    {
        var es = new EarlyStopping(patience: 3, minDelta: 0.01);

        Assert.IsFalse(es.Check(1.0));
        Assert.IsFalse(es.Check(1.0));  // stale #1
        Assert.IsFalse(es.Check(1.0));  // stale #2
        Assert.IsFalse(es.Check(0.5));  // improvement! resets counter
        Assert.IsFalse(es.Check(0.5));  // stale #1
        Assert.IsFalse(es.Check(0.5));  // stale #2
        Assert.IsTrue(es.Check(0.5));   // stale #3 → trigger
    }

    [TestMethod]
    public void EarlyStopping_RespectsMinDelta()
    {
        var es = new EarlyStopping(patience: 2, minDelta: 0.1);

        Assert.IsFalse(es.Check(1.0));      // new best
        Assert.IsFalse(es.Check(0.95));     // improvement < minDelta → stale #1
        Assert.IsTrue(es.Check(0.92));      // still < minDelta → stale #2 → trigger
    }

    [TestMethod]
    public void EarlyStopping_BestValue_TrackedCorrectly()
    {
        var es = new EarlyStopping(patience: 5);

        es.Check(1.0);
        es.Check(0.8);
        es.Check(0.5);
        es.Check(0.7); // worse

        Assert.AreEqual(0.5, es.BestValue, 1e-12);
    }

    [TestMethod]
    public void EarlyStopping_Reset_ClearsState()
    {
        var es = new EarlyStopping(patience: 2);

        es.Check(0.5);
        es.Check(0.5);
        Assert.IsTrue(es.Check(0.5)); // triggered

        es.Reset();
        Assert.IsFalse(es.ShouldStop);
        Assert.AreEqual(double.MaxValue, es.BestValue);
    }

    // ═══════════════════════════════════════════════════════════
    // MaxIterationsOrTolerance tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void MaxIterationsOrTolerance_StopsAtMaxIterations()
    {
        var c = new MaxIterationsOrTolerance(maxIterations: 10, tolerance: 1e-10);

        Assert.IsFalse(c.HasConverged(5, 1.0, 1.0));
        Assert.IsTrue(c.HasConverged(10, 1.0, 1.0));
    }

    [TestMethod]
    public void MaxIterationsOrTolerance_StopsAtGradientTolerance()
    {
        var c = new MaxIterationsOrTolerance(maxIterations: 10000, tolerance: 1e-5);

        Assert.IsFalse(c.HasConverged(0, 1.0, 0.1));
        Assert.IsTrue(c.HasConverged(0, 1.0, 1e-6));
    }

    // ═══════════════════════════════════════════════════════════
    // LearningRateSchedule tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void LearningRateSchedule_Constant_ReturnsInitialLr()
    {
        var sched = new LearningRateSchedule(0.01, LearningRateSchedule.ScheduleType.Constant);

        Assert.AreEqual(0.01, sched.GetLearningRate(0), 1e-12);
        Assert.AreEqual(0.01, sched.GetLearningRate(100), 1e-12);
    }

    [TestMethod]
    public void LearningRateSchedule_ExponentialDecay_Decreases()
    {
        var sched = new LearningRateSchedule(0.1,
            LearningRateSchedule.ScheduleType.ExponentialDecay,
            decayRate: 0.5, decaySteps: 100);

        double lr0 = sched.GetLearningRate(0);
        double lr50 = sched.GetLearningRate(50);
        double lr100 = sched.GetLearningRate(100);

        Assert.AreEqual(0.1, lr0, 1e-6);
        Assert.IsTrue(lr50 < lr0);
        Assert.IsTrue(lr100 < lr50);
        Assert.AreEqual(0.05, lr100, 1e-6); // 0.1 * 0.5^1
    }

    [TestMethod]
    public void LearningRateSchedule_CosineAnnealing_ReachesMinAtEnd()
    {
        double minLr = 1e-5;
        var sched = new LearningRateSchedule(0.1,
            LearningRateSchedule.ScheduleType.CosineAnnealing,
            decaySteps: 100, minLr: minLr);

        double lrEnd = sched.GetLearningRate(100);
        Assert.AreEqual(minLr, lrEnd, 1e-8);
    }
}
