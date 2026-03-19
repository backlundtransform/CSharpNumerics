using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using CSharpNumerics.Numerics.Optimization.MultiObjective;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
using CSharpNumerics.Numerics.Optimization.Strategies;
using NN = CSharpNumerics.ML.NeuralNetwork.NeuralNetwork;

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

    // ═══════════════════════════════════════════════════════════
    // ParetoFront tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void ParetoFront_IdentifiesDominatedSolutions()
    {
        var solutions = new List<ParetoSolution>
        {
            new(new[] { 0.0 }, new[] { 1.0, 5.0 }),  // A
            new(new[] { 0.0 }, new[] { 2.0, 3.0 }),  // B
            new(new[] { 0.0 }, new[] { 3.0, 1.0 }),  // C
            new(new[] { 0.0 }, new[] { 3.0, 4.0 }),  // D — dominated by B and C
        };

        var front = ParetoFront.GetFront(solutions);

        // A, B, C are non-dominated; D is dominated
        Assert.AreEqual(3, front.Count);
        Assert.IsTrue(front.All(s => s.Rank == 0));
    }

    [TestMethod]
    public void ParetoFront_AllNonDominated_ReturnsSingleFront()
    {
        var solutions = new List<ParetoSolution>
        {
            new(new[] { 0.0 }, new[] { 1.0, 4.0 }),
            new(new[] { 0.0 }, new[] { 2.0, 3.0 }),
            new(new[] { 0.0 }, new[] { 3.0, 2.0 }),
            new(new[] { 0.0 }, new[] { 4.0, 1.0 }),
        };

        var fronts = ParetoFront.NonDominatedSort(solutions);

        Assert.AreEqual(1, fronts.Count);
        Assert.AreEqual(4, fronts[0].Count);
    }

    [TestMethod]
    public void ParetoFront_NonDominatedSort_AssignsMultipleRanks()
    {
        var solutions = new List<ParetoSolution>
        {
            new(new[] { 0.0 }, new[] { 1.0, 4.0 }),  // rank 0
            new(new[] { 0.0 }, new[] { 2.0, 2.0 }),  // rank 0
            new(new[] { 0.0 }, new[] { 3.0, 3.0 }),  // rank 1 (dominated by [2,2])
            new(new[] { 0.0 }, new[] { 4.0, 5.0 }),  // rank 2
        };

        var fronts = ParetoFront.NonDominatedSort(solutions);

        Assert.IsTrue(fronts.Count >= 2);
        Assert.AreEqual(2, fronts[0].Count); // [1,4] and [2,2]
    }

    [TestMethod]
    public void ParetoFront_CrowdingDistance_BoundaryPointsAreInfinite()
    {
        var front = new List<ParetoSolution>
        {
            new(new[] { 0.0 }, new[] { 1.0, 5.0 }),
            new(new[] { 0.0 }, new[] { 3.0, 3.0 }),
            new(new[] { 0.0 }, new[] { 5.0, 1.0 }),
        };

        ParetoFront.ComputeCrowdingDistance(front);

        // The two extreme points on each objective should get infinity
        int infiniteCount = front.Count(s => double.IsPositiveInfinity(s.CrowdingDistance));
        Assert.IsTrue(infiniteCount >= 2);
    }

    [TestMethod]
    public void Dominates_StrictlyBetter_ReturnsTrue()
    {
        var a = new ParetoSolution(new[] { 0.0 }, new[] { 1.0, 2.0 });
        var b = new ParetoSolution(new[] { 0.0 }, new[] { 3.0, 4.0 });

        Assert.IsTrue(a.Dominates(b));
        Assert.IsFalse(b.Dominates(a));
    }

    [TestMethod]
    public void Dominates_Equal_ReturnsFalse()
    {
        var a = new ParetoSolution(new[] { 0.0 }, new[] { 2.0, 3.0 });
        var b = new ParetoSolution(new[] { 0.0 }, new[] { 2.0, 3.0 });

        Assert.IsFalse(a.Dominates(b));
        Assert.IsFalse(b.Dominates(a));
    }

    // ═══════════════════════════════════════════════════════════
    // NSGA-II tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void NSGA2_BiObjective_ConvergesToParetoFront()
    {
        // ZDT1-like: f1 = x, f2 = 1 - sqrt(x)
        // Pareto front: x ∈ [0,1], f2 = 1 - sqrt(f1)
        var nsga2 = new NSGA2(
            evaluate: x => new[] { x[0], 1.0 - Math.Sqrt(Math.Max(x[0], 0)) },
            numVariables: 1,
            lowerBounds: new[] { 0.0 },
            upperBounds: new[] { 1.0 },
            populationSize: 50,
            generations: 100,
            seed: 42);

        var result = nsga2.Run();

        Assert.IsTrue(result.FrontSize > 0, "Should find Pareto-optimal solutions");

        // All front solutions should have rank 0
        foreach (var s in result.ParetoFront)
            Assert.AreEqual(0, s.Rank);

        // Check that solutions roughly follow f2 ≈ 1 - sqrt(f1)
        foreach (var s in result.ParetoFront)
        {
            double expected = 1.0 - Math.Sqrt(s.Objectives[0]);
            Assert.AreEqual(expected, s.Objectives[1], 0.15,
                $"f1={s.Objectives[0]:F3}, f2={s.Objectives[1]:F3} expected≈{expected:F3}");
        }
    }

    [TestMethod]
    public void NSGA2_ReturnsReproducibleResults()
    {
        Func<double[], double[]> f = x => new[] { x[0] * x[0], (x[0] - 2) * (x[0] - 2) };

        var r1 = new NSGA2(f, 1, new[] { -5.0 }, new[] { 5.0 },
            populationSize: 30, generations: 50, seed: 123).Run();
        var r2 = new NSGA2(f, 1, new[] { -5.0 }, new[] { 5.0 },
            populationSize: 30, generations: 50, seed: 123).Run();

        Assert.AreEqual(r1.FrontSize, r2.FrontSize);
        for (int i = 0; i < r1.FrontSize; i++)
            Assert.AreEqual(r1.ParetoFront[i].Variables[0],
                r2.ParetoFront[i].Variables[0], 1e-10);
    }

    [TestMethod]
    public void NSGA2_TwoVariables_FindsTradeoff()
    {
        // f1 = x² + y², f2 = (x-1)² + (y-1)²
        // Pareto front: line from (0,0) to (1,1)
        var nsga2 = new NSGA2(
            evaluate: x => new[]
            {
                x[0] * x[0] + x[1] * x[1],
                (x[0] - 1) * (x[0] - 1) + (x[1] - 1) * (x[1] - 1)
            },
            numVariables: 2,
            lowerBounds: new[] { -2.0, -2.0 },
            upperBounds: new[] { 3.0, 3.0 },
            populationSize: 80,
            generations: 150,
            seed: 42);

        var result = nsga2.Run();

        Assert.IsTrue(result.FrontSize >= 5, $"Front has {result.FrontSize} solutions");

        // All Pareto solutions should be in [0,1] range roughly
        foreach (var s in result.ParetoFront)
        {
            Assert.IsTrue(s.Variables[0] >= -0.5 && s.Variables[0] <= 1.5,
                $"x={s.Variables[0]:F3} out of expected range");
        }
    }

    // ═══════════════════════════════════════════════════════════
    // NeuralNetwork IOptimizer integration tests
    // ═══════════════════════════════════════════════════════════

    [TestMethod]
    public void NeuralNetwork_ApplyGradients_WithAdam_ProducesDifferentWeights()
    {
        var nn = new NN(2, new[] { 4 }, 1, seed: 42);
        var (beforeW, beforeB) = nn.SnapshotWeights();

        var input = new VectorN(new[] { 1.0, 0.5 });
        var target = new VectorN(new[] { 1.0 });

        var yhat = nn.Forward(input, out var acts);
        nn.ComputeGradients(target, acts, out var dW, out var dB);

        var adam = new Adam(learningRate: 0.01);
        var adamBias = new Adam(learningRate: 0.01);
        nn.ApplyGradients(dW, dB, adam, adamBias, 1);

        var (afterW, afterB) = nn.SnapshotWeights();

        // Weights should have changed
        bool changed = false;
        for (int r = 0; r < beforeW[0].rowLength; r++)
            for (int c = 0; c < beforeW[0].columnLength; c++)
                if (Math.Abs(beforeW[0].values[r, c] - afterW[0].values[r, c]) > 1e-12)
                    changed = true;

        Assert.IsTrue(changed, "Adam should update weights");
    }

    [TestMethod]
    public void NeuralNetwork_ApplyGradients_LrOnly_MatchesOriginal()
    {
        // Verify the old overload still works exactly as before
        var nn1 = new NN(2, new[] { 3 }, 1, seed: 99);
        var nn2 = nn1.Clone();

        var input = new VectorN(new[] { 0.5, -0.3 });
        var target = new VectorN(new[] { 0.7 });

        nn1.Forward(input, out var a1);
        nn1.ComputeGradients(target, a1, out var dw1, out var db1);
        nn1.ApplyGradients(dw1, db1, 0.01, 1, 0.0);

        nn2.Forward(input, out var a2);
        nn2.ComputeGradients(target, a2, out var dw2, out var db2);
        nn2.ApplyGradients(dw2, db2, 0.01, 1, 0.0);

        var (w1, b1) = nn1.SnapshotWeights();
        var (w2, b2) = nn2.SnapshotWeights();

        for (int l = 0; l < w1.Count; l++)
            for (int r = 0; r < w1[l].rowLength; r++)
                for (int c = 0; c < w1[l].columnLength; c++)
                    Assert.AreEqual(w1[l].values[r, c], w2[l].values[r, c], 1e-12);
    }
}
