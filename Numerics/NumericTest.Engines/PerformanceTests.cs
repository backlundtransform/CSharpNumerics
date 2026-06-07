using CSharpNumerics.Engines.Game;
using CSharpNumerics.Engines.Game.BroadPhase;
using CSharpNumerics.Engines.Game.Fluids;
using CSharpNumerics.Engines.Game.Performance;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics.Objects;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace NumericTest;

[TestClass]
public class PerformanceTests
{
    // ═══════════════════════════════════════════════════════════════
    //  SIMD Math
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void SimdMath_Add_ProducesCorrectResults()
    {
        var a = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };
        var b = new double[] { 10, 20, 30, 40, 50, 60, 70, 80 };
        var result = new double[8];

        SimdMath.Add(a, b, result);

        Assert.AreEqual(11, result[0]);
        Assert.AreEqual(22, result[1]);
        Assert.AreEqual(88, result[7]);
    }

    [TestMethod]
    public void SimdMath_Scale_ProducesCorrectResults()
    {
        var a = new double[] { 1, 2, 3, 4, 5 };
        var result = new double[5];

        SimdMath.Scale(a, 3.0, result);

        Assert.AreEqual(3.0, result[0]);
        Assert.AreEqual(6.0, result[1]);
        Assert.AreEqual(15.0, result[4]);
    }

    [TestMethod]
    public void SimdMath_ScaleAdd_AccumulatesCorrectly()
    {
        var a = new double[] { 1, 2, 3, 4 };
        var result = new double[] { 10, 20, 30, 40 };

        SimdMath.ScaleAdd(a, 2.0, result);

        Assert.AreEqual(12.0, result[0]);
        Assert.AreEqual(24.0, result[1]);
        Assert.AreEqual(36.0, result[2]);
        Assert.AreEqual(48.0, result[3]);
    }

    [TestMethod]
    public void SimdMath_Dot_ProducesCorrectResult()
    {
        var a = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };
        var b = new double[] { 8, 7, 6, 5, 4, 3, 2, 1 };

        double dot = SimdMath.Dot(a, b);

        // 1*8 + 2*7 + 3*6 + 4*5 + 5*4 + 6*3 + 7*2 + 8*1 = 120
        Assert.AreEqual(120.0, dot);
    }

    [TestMethod]
    public void SimdMath_Decay_AppliesFactor()
    {
        var a = new double[] { 100, 200, 300, 400, 500, 600, 700, 800 };

        SimdMath.Decay(a, 0.5);

        Assert.AreEqual(50.0, a[0]);
        Assert.AreEqual(100.0, a[1]);
        Assert.AreEqual(400.0, a[7]);
    }

    [TestMethod]
    public void SimdMath_Subtract_ProducesCorrectResults()
    {
        var a = new double[] { 10, 20, 30, 40 };
        var b = new double[] { 1, 2, 3, 4 };
        var result = new double[4];

        SimdMath.Subtract(a, b, result);

        Assert.AreEqual(9.0, result[0]);
        Assert.AreEqual(18.0, result[1]);
        Assert.AreEqual(36.0, result[3]);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Span / stackalloc — cached Project arrays
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidSolver3D_CachedProject_StillProducesCorrectResults()
    {
        var config = new FluidConfig { GridX = 16, GridY = 16, GridZ = 16 };
        var solver = new GameFluidSolver3D(config);

        // Add emitter and step multiple times — should not crash
        solver.AddEmitter(new FluidEmitter(new Vector(8, 8, 8), densityRate: 5));
        for (int i = 0; i < 20; i++)
            solver.Step(0.016);

        // Density should have accumulated near emitter
        double d = solver.GetDensity(8, 8, 8);
        Assert.IsTrue(d > 0, "Density should be positive near emitter after steps with cached arrays");
    }

    [TestMethod]
    public void FluidSolver2D_CachedProject_StillProducesCorrectResults()
    {
        var config = new FluidConfig { GridX = 32, GridY = 32 };
        var solver = new GameFluidSolver2D(config);

        solver.AddEmitter(new FluidEmitter(new Vector(16, 16, 0), densityRate: 5));
        for (int i = 0; i < 20; i++)
            solver.Step(0.016);

        double d = solver.GetDensity(16, 16);
        Assert.IsTrue(d > 0, "Density should be positive near emitter");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multi-threaded broad phase
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ParallelBroadPhase_FindsSamePairsAsBruteForce()
    {
        var rng = new Random(42);
        int count = 200;
        var bodies = new RigidBody[count];
        var radii = new double[count];

        for (int i = 0; i < count; i++)
        {
            bodies[i] = RigidBody.CreateSolidSphere(mass: 1, radius: 1);
            bodies[i].Position = new Vector(
                rng.NextDouble() * 20 - 10,
                rng.NextDouble() * 20 - 10,
                rng.NextDouble() * 20 - 10);
            radii[i] = 1.0;
        }

        var bruteForce = new BruteForceBroadPhase();
        var parallel = new ParallelBroadPhase(threadCount: 4, singleThreadThreshold: 10);

        var bfResults = new List<(int a, int b)>();
        var parResults = new List<(int a, int b)>();

        bruteForce.FindPairs(bodies, radii, count, bfResults);
        parallel.FindPairs(bodies, radii, count, parResults);

        // Normalize pairs (smaller index first) and sort for comparison
        var bfSet = NormalizePairs(bfResults);
        var parSet = NormalizePairs(parResults);

        Assert.AreEqual(bfSet.Count, parSet.Count,
            $"Pair count mismatch: brute force={bfSet.Count}, parallel={parSet.Count}");

        foreach (var pair in bfSet)
            Assert.IsTrue(parSet.Contains(pair), $"Missing pair ({pair.a}, {pair.b})");
    }

    [TestMethod]
    public void ParallelBroadPhase_FallsBackSingleThread_SmallScene()
    {
        var bodies = new RigidBody[10];
        var radii = new double[10];
        for (int i = 0; i < 10; i++)
        {
            bodies[i] = RigidBody.CreateSolidSphere(mass: 1, radius: 1);
            bodies[i].Position = new Vector(i * 0.5, 0, 0);
            radii[i] = 1.0;
        }

        var parallel = new ParallelBroadPhase(singleThreadThreshold: 128);
        var results = new List<(int a, int b)>();

        parallel.FindPairs(bodies, radii, 10, results);

        Assert.IsTrue(results.Count > 0, "Should find overlapping pairs in close group");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluid solver thread
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidSolverThread_RunsOnBackground()
    {
        var config = new FluidConfig { GridX = 16, GridY = 16, GridZ = 16 };
        var solver = new GameFluidSolver3D(config);
        solver.AddEmitter(new FluidEmitter(new Vector(8, 8, 8), densityRate: 10));

        using var thread = new FluidSolverThread(solver);

        // Request a few steps
        for (int i = 0; i < 5; i++)
        {
            thread.RequestStep(0.016);
            thread.WaitForStep(5000);
        }

        Assert.AreEqual(5, thread.StepCount);

        // Check density snapshot
        int size = 18 * 18 * 18; // nx+2 * ny+2 * nz+2
        var density = new double[size];
        thread.GetDensitySnapshot(density);

        // Should have some non-zero density
        bool hasDensity = false;
        for (int i = 0; i < density.Length; i++)
            if (density[i] > 0) { hasDensity = true; break; }
        Assert.IsTrue(hasDensity, "Background solver should produce density");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Memory pooling
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void GameArrayPool_RentReturn_DoesNotThrow()
    {
        var arr = GameArrayPool.RentDouble(1000);
        Assert.IsTrue(arr.Length >= 1000);
        GameArrayPool.Return(arr);

        var arr2 = GameArrayPool.RentInt(500);
        Assert.IsTrue(arr2.Length >= 500);
        GameArrayPool.Return(arr2);
    }

    [TestMethod]
    public void GameArrayPool_RentDoubleCleared_IsZeroed()
    {
        var arr = GameArrayPool.RentDoubleCleared(100);
        for (int i = 0; i < 100; i++)
            Assert.AreEqual(0.0, arr[i]);
        GameArrayPool.Return(arr);
    }

    [TestMethod]
    public void GameArrayPool_RentScope_DisposesCorrectly()
    {
        using (var scope = GameArrayPool.RentScope(256))
        {
            Assert.IsTrue(scope.Array.Length >= 256);
            Assert.AreEqual(256, scope.Length);
            scope.Span[0] = 42;
            Assert.AreEqual(42.0, scope.Array[0]);
        }
        // Array returned to pool — no way to verify, but it shouldn't throw
    }

    // ═══════════════════════════════════════════════════════════════
    //  Profile timer
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ProfileTimer_MeasuresTiming()
    {
        var profiler = new ProfileTimer();

        profiler.Begin("Test");
        // Do some work
        double sum = 0;
        for (int i = 0; i < 100000; i++) sum += Math.Sqrt(i);
        profiler.End("Test");

        var stats = profiler.GetStats("Test");
        Assert.AreEqual(1, stats.Count);
        Assert.IsTrue(stats.TotalMs > 0, "Should measure non-zero time");
        Assert.IsTrue(stats.AverageMs > 0);
    }

    [TestMethod]
    public void ProfileTimer_Scope_AutoEndsOnDispose()
    {
        var profiler = new ProfileTimer();

        using (profiler.Scope("AutoTest"))
        {
            double sum = 0;
            for (int i = 0; i < 10000; i++) sum += Math.Sqrt(i);
        }

        var stats = profiler.GetStats("AutoTest");
        Assert.AreEqual(1, stats.Count);
        Assert.IsTrue(stats.TotalMs >= 0);
    }

    [TestMethod]
    public void ProfileTimer_TracksMinMax()
    {
        var profiler = new ProfileTimer();

        for (int t = 0; t < 10; t++)
        {
            profiler.Begin("Multi");
            double sum = 0;
            for (int i = 0; i < 1000; i++) sum += Math.Sqrt(i);
            profiler.End("Multi");
        }

        var stats = profiler.GetStats("Multi");
        Assert.AreEqual(10, stats.Count);
        Assert.IsTrue(stats.MinMs <= stats.AverageMs);
        Assert.IsTrue(stats.MaxMs >= stats.AverageMs);
    }

    // ═══════════════════════════════════════════════════════════════
    //  LOD fluid system
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidLOD_DualResolution_Steps()
    {
        var coarseConfig = new FluidConfig { GridX = 16, GridY = 16, GridZ = 16, PoissonIterations = 5 };
        var fineConfig = new FluidConfig { GridX = 16, GridY = 16, GridZ = 16, PoissonIterations = 5 };

        var lod = new FluidLOD(coarseConfig, fineConfig, scale: 2);

        // Add emitter to coarse solver
        lod.CoarseSolver.AddEmitter(new FluidEmitter(new Vector(8, 8, 8), densityRate: 10));

        // Step without crashing
        for (int i = 0; i < 5; i++)
            lod.Step(0.016);

        // Sample density
        double d = lod.SampleDensity(8, 8, 8);
        Assert.IsTrue(d >= 0, "Density should be non-negative");
    }

    [TestMethod]
    public void FluidLOD_SetFineRegion_DoesNotThrow()
    {
        var coarseConfig = new FluidConfig { GridX = 32, GridY = 32, GridZ = 32, PoissonIterations = 5 };
        var fineConfig = new FluidConfig { GridX = 16, GridY = 16, GridZ = 16, PoissonIterations = 5 };
        var lod = new FluidLOD(coarseConfig, fineConfig, scale: 2);

        lod.SetFineRegion(4, 4, 4);
        lod.Step(0.016);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Deterministic replay
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void SimulationRecorder_RecordAndPlayback()
    {
        var world = new PhysicsWorld();
        var ball = RigidBody.CreateSolidSphere(mass: 1, radius: 0.5);
        ball.Position = new Vector(0, 0, 10);
        world.AddBody(ball, 0.5);

        var recorder = new SimulationRecorder();

        // Record 10 frames
        double dt = 0.01;
        for (int i = 0; i < 10; i++)
        {
            world.Step(dt);
            recorder.RecordFrame(world, i * dt);
        }

        Assert.AreEqual(10, recorder.FrameCount);

        // Playback frame 0 — position should be close to initial
        recorder.ApplyFrame(world, 0);
        var pos0 = world.Body(0).Position;

        // Playback frame 9 — should have fallen
        recorder.ApplyFrame(world, 9);
        var pos9 = world.Body(0).Position;

        Assert.IsTrue(pos9.z < pos0.z, "Ball should have fallen during recorded frames");
    }

    [TestMethod]
    public void SimulationRecorder_BinarySerialization_RoundTrip()
    {
        var world = new PhysicsWorld();
        var ball = RigidBody.CreateSolidSphere(mass: 1, radius: 0.5);
        ball.Position = new Vector(5, 10, 15);
        ball.Velocity = new Vector(1, 2, 3);
        world.AddBody(ball, 0.5);

        var recorder = new SimulationRecorder();
        recorder.RecordFrame(world, 0);
        recorder.RecordFrame(world, 0.01);

        // Write to stream
        var stream = new MemoryStream();
        recorder.WriteTo(stream);

        // Read back
        stream.Position = 0;
        var recorder2 = new SimulationRecorder();
        recorder2.ReadFrom(stream);

        Assert.AreEqual(2, recorder2.FrameCount);

        var frame = recorder2.GetFrame(0);
        Assert.AreEqual(5.0, frame.Bodies[0].Position.x, 1e-10);
        Assert.AreEqual(10.0, frame.Bodies[0].Position.y, 1e-10);
        Assert.AreEqual(15.0, frame.Bodies[0].Position.z, 1e-10);
        Assert.AreEqual(1.0, frame.Bodies[0].Velocity.x, 1e-10);
    }

    [TestMethod]
    public void SimulationRecorder_NextFrame_AdvancesPlayback()
    {
        var world = new PhysicsWorld();
        var ball = RigidBody.CreateSolidSphere(mass: 1, radius: 0.5);
        ball.Position = new Vector(0, 0, 10);
        world.AddBody(ball, 0.5);

        var recorder = new SimulationRecorder();
        for (int i = 0; i < 5; i++)
        {
            world.Step(0.01);
            recorder.RecordFrame(world, i * 0.01);
        }

        recorder.Rewind();
        recorder.ApplyFrame(world, 0);
        Assert.AreEqual(0, recorder.CurrentFrame);

        bool more = recorder.NextFrame(world);
        Assert.IsTrue(more);
        Assert.AreEqual(1, recorder.CurrentFrame);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Profiling benchmarks — measure ms/frame for key scenarios
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Benchmark_Fluid3D_64Cube_PerformanceBaseline()
    {
        var config = new FluidConfig { GridX = 32, GridY = 32, GridZ = 32, PoissonIterations = 10 };
        var solver = new GameFluidSolver3D(config);
        solver.AddEmitter(new FluidEmitter(new Vector(16, 16, 16), densityRate: 10));

        // Warm up
        solver.Step(0.016);

        // Measure
        var sw = Stopwatch.StartNew();
        int steps = 10;
        for (int i = 0; i < steps; i++)
            solver.Step(0.016);
        sw.Stop();

        double msPerStep = sw.Elapsed.TotalMilliseconds / steps;
        // Just log it — no hard assertion on speed (varies by hardware)
        Console.WriteLine($"3D Fluid (32³, 10 Poisson): {msPerStep:F2} ms/step");
        Assert.IsTrue(msPerStep < 10000, "Sanity check: step should not take 10+ seconds");
    }

    [TestMethod]
    public void Benchmark_PhysicsWorld_100Bodies()
    {
        var world = new PhysicsWorld();
        var rng = new Random(42);
        for (int i = 0; i < 100; i++)
        {
            var b = RigidBody.CreateSolidSphere(mass: 1, radius: 0.5);
            b.Position = new Vector(rng.NextDouble() * 10, rng.NextDouble() * 10, rng.NextDouble() * 10);
            b.Velocity = new Vector(rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1);
            world.AddBody(b, 0.5);
        }

        // Warm up
        world.Step(0.01);

        var sw = Stopwatch.StartNew();
        int steps = 100;
        for (int i = 0; i < steps; i++)
            world.Step(0.01);
        sw.Stop();

        double msPerStep = sw.Elapsed.TotalMilliseconds / steps;
        Console.WriteLine($"PhysicsWorld (100 bodies): {msPerStep:F3} ms/step");
        Assert.IsTrue(msPerStep < 1000, "100 bodies should not take 1+ second per step");
    }

    [TestMethod]
    public void Benchmark_BroadPhase_Comparison()
    {
        var rng = new Random(42);
        int count = 200;
        var bodies = new RigidBody[count];
        var radii = new double[count];
        for (int i = 0; i < count; i++)
        {
            bodies[i] = RigidBody.CreateSolidSphere(mass: 1, radius: 1);
            bodies[i].Position = new Vector(
                rng.NextDouble() * 50, rng.NextDouble() * 50, rng.NextDouble() * 50);
            radii[i] = 1.0;
        }

        var results = new List<(int a, int b)>();

        // Brute force
        var bf = new BruteForceBroadPhase();
        var sw1 = Stopwatch.StartNew();
        for (int i = 0; i < 100; i++)
            bf.FindPairs(bodies, radii, count, results);
        sw1.Stop();

        // Sweep and prune
        var sap = new SweepAndPruneBroadPhase();
        var sw2 = Stopwatch.StartNew();
        for (int i = 0; i < 100; i++)
            sap.FindPairs(bodies, radii, count, results);
        sw2.Stop();

        // BVH
        var bvh = new BVHBroadPhase();
        var sw3 = Stopwatch.StartNew();
        for (int i = 0; i < 100; i++)
            bvh.FindPairs(bodies, radii, count, results);
        sw3.Stop();

        Console.WriteLine($"BroadPhase (200 bodies, 100 iters):");
        Console.WriteLine($"  BruteForce:     {sw1.Elapsed.TotalMilliseconds:F2} ms");
        Console.WriteLine($"  SweepAndPrune:  {sw2.Elapsed.TotalMilliseconds:F2} ms");
        Console.WriteLine($"  BVH:            {sw3.Elapsed.TotalMilliseconds:F2} ms");

        Assert.IsTrue(sw2.Elapsed < sw1.Elapsed * 5, "SAP should not be 5× slower than brute force");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static HashSet<(int a, int b)> NormalizePairs(List<(int a, int b)> pairs)
    {
        var set = new HashSet<(int a, int b)>();
        foreach (var (a, b) in pairs)
            set.Add(a < b ? (a, b) : (b, a));
        return set;
    }
}
