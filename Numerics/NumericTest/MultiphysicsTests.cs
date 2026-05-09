using CSharpNumerics.Engines.Multiphysics;
using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Engines.Multiphysics.Export;
using CSharpNumerics.Engines.Multiphysics.MonteCarlo;
using CSharpNumerics.Engines.Multiphysics.Snapshots;
using CSharpNumerics.ML.Clustering;
using CSharpNumerics.ML.Clustering.Algorithms;
using CSharpNumerics.ML.Clustering.Evaluators;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.Statistics.Random;
using CSharpNumerics.Numerics.FiniteDifference;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Engineering;
using CSharpNumerics.Physics.SolidMechanics.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.IO;
using System.Linq;

namespace NumericTest;

[TestClass]
public class MultiphysicsTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Engineering Material
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void EngineeringMaterial_Steel_HasCorrectProperties()
    {
        var steel = EngineeringLibrary.Steel;
        Assert.AreEqual("Steel", steel.Name);
        Assert.AreEqual(200e9, steel.YoungsModulus, 1e6);
        Assert.AreEqual(0.30, steel.PoissonsRatio, 0.001);
        Assert.IsTrue(steel.ThermalDiffusivity > 0);
    }

    [TestMethod]
    public void EngineeringMaterial_ThermalDiffusivity_MatchesFormula()
    {
        var al = EngineeringLibrary.Aluminum;
        double expected = al.ThermalConductivity / (al.Density * al.SpecificHeat);
        Assert.AreEqual(expected, al.ThermalDiffusivity, 1e-12);
    }

    [TestMethod]
    public void EngineeringMaterial_Water_HasViscosity()
    {
        var water = EngineeringLibrary.Water;
        Assert.IsTrue(water.DynamicViscosity > 0);
        Assert.IsTrue(water.KinematicViscosity > 0);
        double expectedNu = water.DynamicViscosity / water.Density;
        Assert.AreEqual(expectedNu, water.KinematicViscosity, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Poisson Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PoissonSolver_Laplace_ConvergesToLinearSolution()
    {
        // ∇²u = 0 with u(left) = 0, u(right) = 100 → linear gradient
        int nx = 20, ny = 20;
        double dx = 1.0 / (nx - 1);
        var grid = new Grid2D(nx, ny, dx);

        var rhs = grid.Zeros();
        var mask = new bool[grid.Length];
        var bcVals = new double[grid.Length];

        // Left = 0, right = 100, top/bottom = linear interpolation
        for (int iy = 0; iy < ny; iy++)
        {
            mask[grid.Index(0, iy)] = true;
            bcVals[grid.Index(0, iy)] = 0;
            mask[grid.Index(nx - 1, iy)] = true;
            bcVals[grid.Index(nx - 1, iy)] = 100;
        }
        for (int ix = 0; ix < nx; ix++)
        {
            double expected = 100.0 * ix / (nx - 1);
            mask[grid.Index(ix, 0)] = true;
            bcVals[grid.Index(ix, 0)] = expected;
            mask[grid.Index(ix, ny - 1)] = true;
            bcVals[grid.Index(ix, ny - 1)] = expected;
        }

        var (solution, iters) = GridOperators.SolvePoisson2D(rhs, grid, mask, bcVals);

        // Interior should be close to linear gradient in x
        for (int ix = 1; ix < nx - 1; ix++)
        {
            double expected = 100.0 * ix / (nx - 1);
            double actual = solution[grid.Index(ix, ny / 2)];
            Assert.AreEqual(expected, actual, 0.5, $"At ix={ix}: expected ~{expected:F1}, got {actual:F1}");
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — HeatPlate
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void HeatPlate_FluentAPI_ProducesResult()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
            .WithBoundary(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 50);

        Assert.AreEqual(MultiphysicsType.HeatPlate, result.Type);
        Assert.IsNotNull(result.Field);
        Assert.IsNotNull(result.Timeline);
        Assert.AreEqual(51, result.Timeline.Count); // initial + 50 steps
        Assert.AreEqual(100.0, result.MaxValue, 0.1);
    }

    [TestMethod]
    public void HeatPlate_BoundaryConditions_AreApplied()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 1.0, height: 1.0, nx: 10, ny: 10)
            .WithBoundary(top: 200.0, bottom: 50.0, left: 50.0, right: 50.0)
            .WithInitialCondition(50.0)
            .Solve(dt: 0.001, steps: 10);

        // Top-right corner should be 200 (top BC) or 50 (right BC)
        // Corners are set by both edges — last write wins (right edge written after top)
        // Top edge value at interior
        double topMid = result.Field[5, 9]; // ix=5, iy=9 (top row)
        Assert.AreEqual(200.0, topMid, 0.01);
    }

    [TestMethod]
    public void HeatPlate_WithSource_IncreasesTemperature()
    {
        var noSource = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
            .WithBoundary(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 50);

        var withSource = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
            .WithBoundary(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0)
            .WithInitialCondition(0.0)
            .AddSource(10, 10, 1e6)
            .Solve(dt: 0.0001, steps: 50);

        Assert.IsTrue(withSource.MaxValue > noSource.MaxValue,
            "Source should increase temperature");
    }

    [TestMethod]
    public void HeatPlate_RunAlias_SameAsSolve()
    {
        var builder = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 10, ny: 10)
            .WithBoundary(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0)
            .WithInitialCondition(0.0);

        // Run() is alias for Solve() — should not throw
        var result = builder.Run(dt: 0.001, steps: 5);
        Assert.IsNotNull(result.Field);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — BeamStress
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void BeamStress_CantileverPointLoad_MatchesAnalytical()
    {
        double P = 1000;   // 1 kN
        double L = 2.0;    // 2 m
        double E = 200e9;  // Steel
        double b = 0.05, h = 0.1;
        double I = b * h * h * h / 12.0;
        double EI = E * I;

        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: L, nodes: 201)
            .WithCrossSection(width: b, height: h)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: L, value: P)
            .Solve();

        Assert.AreEqual(MultiphysicsType.BeamStress, result.Type);
        Assert.AreEqual(201, result.Values.Length);
        Assert.AreEqual(201, result.Positions.Length);

        // Max deflection at free end: δ = PL³/(3EI)
        double expectedMax = P * L * L * L / (3.0 * EI);
        double actualMax = result.Values[result.Values.Length - 1];
        Assert.AreEqual(expectedMax, actualMax, expectedMax * 0.01,
            $"Expected max deflection ~{expectedMax:E3}, got {actualMax:E3}");
    }

    [TestMethod]
    public void BeamStress_SimplySupportedUniformLoad_MatchesAnalytical()
    {
        double q = 500;    // 500 N/m
        double L = 3.0;    // 3 m
        double b = 0.05, h = 0.1;
        double I = b * h * h * h / 12.0;
        double EI = EngineeringLibrary.Steel.YoungsModulus * I;

        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: L, nodes: 101)
            .WithCrossSection(width: b, height: h)
            .WithBoundary(BeamSupport.SimplySupported)
            .WithSource(q)
            .Solve();

        // Max deflection at midspan: δ = 5qL⁴/(384EI)
        double expectedMax = 5.0 * q * L * L * L * L / (384.0 * EI);
        double midIdx = result.Values.Length / 2;
        double actualMid = result.Values[(int)midIdx];
        Assert.AreEqual(expectedMax, actualMid, expectedMax * 0.02,
            $"Expected midspan deflection ~{expectedMax:E3}, got {actualMid:E3}");
    }

    [TestMethod]
    public void BeamStress_HasMomentShearStress()
    {
        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 1.0, nodes: 51)
            .WithCrossSection(width: 0.02, height: 0.04)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: 1.0, value: 500)
            .Solve();

        Assert.IsNotNull(result.BendingMoment);
        Assert.IsNotNull(result.ShearForce);
        Assert.IsNotNull(result.Stress);
        Assert.AreEqual(51, result.BendingMoment.Length);

        // Moment at the fixed end should be M = P·L = 500·1 = 500 N·m
        Assert.AreEqual(500.0, result.BendingMoment[0], 1.0);
    }

    [TestMethod]
    public void BeamStress_CircularCrossSection_ComputesI()
    {
        double r = 0.02; // 20 mm radius
        double expectedI = Math.PI * Math.Pow(r, 4) / 4.0;

        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(length: 1.0, nodes: 51)
            .WithCrossSection(radius: r)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: 1.0, value: 100)
            .Solve();

        // Verify deflection is consistent with circular I
        double P = 100, L = 1.0, EI = EngineeringLibrary.Aluminum.YoungsModulus * expectedI;
        double expectedDelta = P * L * L * L / (3.0 * EI);
        Assert.AreEqual(expectedDelta, result.Values[result.Values.Length - 1], expectedDelta * 0.01);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — ElectricField
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ElectricField_ParallelPlateCapacitor_LinearPotential()
    {
        // Left = 0 V, Right = 100 V, top/bottom = Neumann-like
        // Expect roughly linear potential gradient in x
        var result = SimulationType.Create(MultiphysicsType.ElectricField)
            .WithMaterial(EngineeringLibrary.Air)
            .WithGeometry(width: 1.0, height: 1.0, nx: 30, ny: 30)
            .WithBoundary(top: 50.0, bottom: 50.0, left: 0.0, right: 100.0)
            .Solve(maxIterations: 20000, tolerance: 1e-6);

        Assert.AreEqual(MultiphysicsType.ElectricField, result.Type);
        Assert.IsNotNull(result.Field);
        Assert.IsNotNull(result.Ex);
        Assert.IsNotNull(result.Ey);

        // Middle of domain: potential should be ~50 V
        double midPotential = result.Field[15, 15];
        Assert.AreEqual(50.0, midPotential, 5.0,
            $"Mid-domain potential should be ~50 V, got {midPotential:F1}");
    }

    [TestMethod]
    public void ElectricField_WithCharge_CreatesField()
    {
        var result = SimulationType.Create(MultiphysicsType.ElectricField)
            .WithMaterial(EngineeringLibrary.Air)
            .WithGeometry(width: 1.0, height: 1.0, nx: 20, ny: 20)
            .WithBoundary(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0)
            .AddSource(10, 10, 1e-6)
            .Solve(maxIterations: 10000, tolerance: 1e-8);

        // With a charge inside, the potential should be non-zero somewhere
        bool hasNonZero = false;
        for (int ix = 1; ix < 19; ix++)
            for (int iy = 1; iy < 19; iy++)
                if (Math.Abs(result.Field[ix, iy]) > 1e-10) hasNonZero = true;

        Assert.IsTrue(hasNonZero, "Charge source should create non-zero potential");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — PipeFlow
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PipeFlow_HagenPoiseuille_ConvergesToParabolic()
    {
        // Analytical steady-state: v(r) = (R² − r²)·(-dP/dx)/(4μ)
        var water = EngineeringLibrary.Water;
        double R = 0.005;   // 5 mm radius — smaller R → shorter diffusion time
        double dPdx = -100; // Pa/m (drives flow in +x)
        int n = 21;

        // Diffusion time scale τ = R²/ν ≈ 25 s; run for ~2τ
        var result = SimulationType.Create(MultiphysicsType.PipeFlow)
            .WithMaterial(water)
            .WithGeometry(length: 1.0, radius: R, nodes: n)
            .WithBoundary(pressureGradient: dPdx)
            .Solve(dt: 0.02, steps: 3000);

        Assert.AreEqual(MultiphysicsType.PipeFlow, result.Type);
        Assert.AreEqual(n, result.Values.Length);

        // Centreline velocity should be the max
        double vCenter = result.Values[0];
        double vWall = result.Values[n - 1];
        Assert.AreEqual(0.0, vWall, 1e-10, "No-slip at wall");
        Assert.IsTrue(vCenter > 0, "Flow should be in +x direction");

        // Analytical centreline: v_max = R²·(-dP/dx)/(4μ)
        double expectedCenter = R * R * (-dPdx) / (4.0 * water.DynamicViscosity);
        Assert.AreEqual(expectedCenter, vCenter, expectedCenter * 0.15,
            $"Centreline velocity: expected {expectedCenter:E3}, got {vCenter:E3}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Validation — missing config throws
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void Solve_WithoutMaterial_Throws()
    {
        SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithGeometry(width: 1, height: 1, nx: 10, ny: 10)
            .Solve();
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void BeamStress_WithoutSupport_Throws()
    {
        SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 1.0, nodes: 10)
            .WithCrossSection(width: 0.05, height: 0.1)
            .Solve();
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void BeamStress_WithoutCrossSection_Throws()
    {
        SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 1.0, nodes: 10)
            .WithBoundary(BeamSupport.Cantilever)
            .Solve();
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 3 — Snapshots
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FieldSnapshot_RoundTrips2DArray()
    {
        var field = new double[5, 4];
        field[2, 3] = 42.0;
        field[0, 0] = -1.5;

        var snap = new FieldSnapshot(field, MultiphysicsType.HeatPlate, 1.5, 3, 0.01, 0.01);

        Assert.AreEqual(5, snap.Nx);
        Assert.AreEqual(4, snap.Ny);
        Assert.AreEqual(1.5, snap.Time, 1e-12);
        Assert.AreEqual(3, snap.StepIndex);
        Assert.AreEqual(20, snap.Count);
        Assert.AreEqual(42.0, snap[2, 3], 1e-12);
        Assert.AreEqual(-1.5, snap[0, 0], 1e-12);
        Assert.AreEqual(42.0, snap.Max(), 1e-12);
        Assert.AreEqual(-1.5, snap.Min(), 1e-12);

        var roundTrip = snap.ToArray();
        Assert.AreEqual(42.0, roundTrip[2, 3], 1e-12);
        Assert.AreEqual(-1.5, roundTrip[0, 0], 1e-12);
    }

    [TestMethod]
    public void BeamSnapshot_FromResult_CapturesAllCurves()
    {
        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 2.0, nodes: 51)
            .WithCrossSection(width: 0.05, height: 0.1)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: 2.0, value: 1000)
            .Solve();

        var snap = BeamSnapshot.FromResult(result, BeamSupport.Cantilever, 2.0);

        Assert.AreEqual(51, snap.NodeCount);
        Assert.AreEqual(2.0, snap.Length, 1e-12);
        Assert.AreEqual(BeamSupport.Cantilever, snap.Support);
        Assert.IsTrue(snap.MaxDeflection > 0);
        Assert.IsTrue(snap.MaxStress > 0);
        Assert.AreEqual(result.Values.Length, snap.Deflection.Length);
        Assert.AreEqual(result.BendingMoment.Length, snap.BendingMoment.Length);
    }

    [TestMethod]
    public void SimulationTimeline_FromResult_HasCorrectCount()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 10, ny: 10)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 20);

        double dx = 0.1 / 10, dy = 0.1 / 10;
        var timeline = SimulationTimeline.FromResult(result, dt: 0.0001, dx: dx, dy: dy);

        Assert.AreEqual(21, timeline.Count); // initial + 20 steps
        Assert.AreEqual(MultiphysicsType.HeatPlate, timeline.Type);
        Assert.AreEqual(0.0, timeline.StartTime, 1e-12);
        Assert.AreEqual(0.002, timeline.EndTime, 1e-7);
        Assert.AreEqual(10, timeline[0].Nx);
        Assert.AreEqual(10, timeline[0].Ny);
    }

    [TestMethod]
    public void SimulationTimeline_Interpolation_BlendsBetweenFrames()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(50.0)
            .Solve(dt: 0.0001, steps: 10);

        double dx = 0.1 / 5, dy = 0.1 / 5;
        var timeline = SimulationTimeline.FromResult(result, dt: 0.0001, dx: dx, dy: dy);

        // Interpolate at midpoint between step 0 and step 1
        double midTime = 0.00005;
        var interp = timeline.InterpolateAt(midTime);
        Assert.AreEqual(5, interp.GetLength(0));
        Assert.AreEqual(5, interp.GetLength(1));

        // Value should be between step 0 and step 1 for an interior cell
        double v0 = timeline[0][2, 2];
        double v1 = timeline[1][2, 2];
        double vi = interp[2, 2];
        double expectedMid = (v0 + v1) / 2.0;
        Assert.AreEqual(expectedMid, vi, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 3 — Binary Export
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void BinaryExporter_2DTimeline_RoundTrips()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 0.1, height: 0.1, nx: 8, ny: 8)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 5);

        double dx = 0.1 / 8, dy = 0.1 / 8;
        var timeline = SimulationTimeline.FromResult(result, dt: 0.0001, dx: dx, dy: dy);

        string path = Path.Combine(Path.GetTempPath(), "test_heat.mphy");
        try
        {
            MultiphysicsBinaryExporter.Save(timeline, path);
            Assert.IsTrue(File.Exists(path));

            var header = MultiphysicsBinaryExporter.ReadHeader(path);
            Assert.AreEqual(1, header.Version);
            Assert.AreEqual(MultiphysicsType.HeatPlate, header.Type);
            Assert.AreEqual(8, header.Nx);
            Assert.AreEqual(8, header.Ny);
            Assert.AreEqual(6, header.TimeStepCount);
            Assert.AreEqual(1, header.LayerCount);
            Assert.IsTrue(header.Is2D);

            var data = MultiphysicsBinaryExporter.Read(path);
            Assert.IsNull(data.Positions);
            Assert.AreEqual(6, data.Layers.Length); // 6 steps × 1 layer
            Assert.AreEqual(64, data.Layers[0].Length); // 8×8

            // Verify last frame top-row value is ~100
            float topMid = data.Layers[5][7 * 8 + 4]; // iy=7, ix=4
            Assert.AreEqual(100.0f, topMid, 0.1f);
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    [TestMethod]
    public void BinaryExporter_BeamSnapshot_RoundTrips()
    {
        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 1.0, nodes: 21)
            .WithCrossSection(width: 0.02, height: 0.04)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: 1.0, value: 500)
            .Solve();

        var snap = BeamSnapshot.FromResult(result, BeamSupport.Cantilever, 1.0);

        string path = Path.Combine(Path.GetTempPath(), "test_beam.mphy");
        try
        {
            MultiphysicsBinaryExporter.Save(snap, path);
            Assert.IsTrue(File.Exists(path));

            var header = MultiphysicsBinaryExporter.ReadHeader(path);
            Assert.AreEqual(MultiphysicsType.BeamStress, header.Type);
            Assert.AreEqual(21, header.Nx);
            Assert.AreEqual(0, header.Ny);
            Assert.AreEqual(4, header.LayerCount);
            Assert.IsFalse(header.Is2D);

            var data = MultiphysicsBinaryExporter.Read(path);
            Assert.IsNotNull(data.Positions);
            Assert.AreEqual(21, data.Positions.Length);
            Assert.AreEqual(4, data.Layers.Length); // 4 layers × 1 step
            Assert.AreEqual(21, data.Layers[0].Length);

            // Deflection at tip should match
            float tipDeflection = data.Layers[0][20];
            Assert.AreEqual((float)result.Values[20], tipDeflection, 1e-4f);
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 3 — JSON Export
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void JsonExporter_SimulationResult_ContainsType()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.001, steps: 2);

        string json = MultiphysicsJsonExporter.ToJson(result);
        Assert.IsTrue(json.Contains("\"type\":\"HeatPlate\""));
        Assert.IsTrue(json.Contains("\"field\":["));
        Assert.IsTrue(json.Contains("\"maxValue\""));
    }

    [TestMethod]
    public void JsonExporter_BeamResult_ContainsAllCurves()
    {
        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 1.0, nodes: 11)
            .WithCrossSection(width: 0.05, height: 0.1)
            .WithBoundary(BeamSupport.SimplySupported)
            .WithSource(500)
            .Solve();

        string json = MultiphysicsJsonExporter.ToJson(result);
        Assert.IsTrue(json.Contains("\"type\":\"BeamStress\""));
        Assert.IsTrue(json.Contains("\"positions\":["));
        Assert.IsTrue(json.Contains("\"values\":["));
        Assert.IsTrue(json.Contains("\"bendingMoment\":["));
        Assert.IsTrue(json.Contains("\"shearForce\":["));
        Assert.IsTrue(json.Contains("\"stress\":["));
    }

    [TestMethod]
    public void JsonExporter_BeamSnapshot_ContainsMetadata()
    {
        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(length: 2.0, nodes: 11)
            .WithCrossSection(radius: 0.01)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: 2.0, value: 100)
            .Solve();

        var snap = BeamSnapshot.FromResult(result, BeamSupport.Cantilever, 2.0);
        var meta = new ExportMetadata { Simulation = "Test beam", Unit = "m" };
        string json = MultiphysicsJsonExporter.ToJson(snap, meta);

        Assert.IsTrue(json.Contains("\"type\":\"BeamStress\""));
        Assert.IsTrue(json.Contains("\"support\":\"Cantilever\""));
        Assert.IsTrue(json.Contains("\"simulation\":\"Test beam\""));
        Assert.IsTrue(json.Contains("\"unit\":\"m\""));
        Assert.IsTrue(json.Contains("\"deflection\":["));
    }

    [TestMethod]
    public void JsonExporter_Timeline_HasAllSteps()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 3);

        double dx = 0.1 / 5, dy = 0.1 / 5;
        var timeline = SimulationTimeline.FromResult(result, dt: 0.0001, dx: dx, dy: dy);

        string json = MultiphysicsJsonExporter.ToJson(timeline);
        Assert.IsTrue(json.Contains("\"timeSteps\":["));
        Assert.IsTrue(json.Contains("\"stepIndex\":0"));
        Assert.IsTrue(json.Contains("\"stepIndex\":3"));
        Assert.IsTrue(json.Contains("\"dimensions\":{"));
    }

    [TestMethod]
    public void JsonExporter_ElectricField_ContainsExEy()
    {
        var result = SimulationType.Create(MultiphysicsType.ElectricField)
            .WithMaterial(EngineeringLibrary.Air)
            .WithGeometry(width: 1.0, height: 1.0, nx: 10, ny: 10)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 100)
            .Solve(maxIterations: 5000, tolerance: 1e-4);

        string json = MultiphysicsJsonExporter.ToJson(result);
        Assert.IsTrue(json.Contains("\"type\":\"ElectricField\""));
        Assert.IsTrue(json.Contains("\"field\":["));
        Assert.IsTrue(json.Contains("\"ex\":["));
        Assert.IsTrue(json.Contains("\"ey\":["));
    }

    [TestMethod]
    public void JsonExporter_SaveToFile_CreatesFile()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.001, steps: 1);

        string path = Path.Combine(Path.GetTempPath(), "test_heat.json");
        try
        {
            MultiphysicsJsonExporter.Save(result, path);
            Assert.IsTrue(File.Exists(path));
            string content = File.ReadAllText(path);
            Assert.IsTrue(content.StartsWith("{"));
            Assert.IsTrue(content.Contains("\"type\":\"HeatPlate\""));
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 4 — Monte Carlo
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ParameterVariation_FluentSetters_StoreRanges()
    {
        var v = new ParameterVariation()
            .ThermalConductivity(40, 60)
            .YoungsModulus(180e9, 220e9)
            .Load(800, 1200);

        Assert.AreEqual(40, v.ThermalConductivityMin);
        Assert.AreEqual(60, v.ThermalConductivityMax);
        Assert.AreEqual(180e9, v.YoungsModulusMin);
        Assert.AreEqual(220e9, v.YoungsModulusMax);
        Assert.AreEqual(800, v.LoadMin);
        Assert.AreEqual(1200, v.LoadMax);
        Assert.IsTrue(v.HasVariation);
    }

    [TestMethod]
    public void MonteCarlo_HeatPlate_RunBatch_ProducesScenarioMatrix()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(40, 60)
            .BoundaryTemperature(90, 110);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(10, seed: 42);

        Assert.AreEqual(10, result.Iterations);
        Assert.AreEqual(25, result.FeatureCount); // 5×5
        Assert.AreEqual(10, result.Results.Count);
        Assert.AreEqual(MultiphysicsType.HeatPlate, result.Type);
    }

    [TestMethod]
    public void MonteCarlo_Beam_RunBatch_VariesOutput()
    {
        var baseline = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: 1.0, nodes: 11)
            .WithCrossSection(width: 0.05, height: 0.1)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: 1.0, value: 1000);

        var variation = new ParameterVariation()
            .YoungsModulus(180e9, 220e9)
            .Load(800, 1200);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(8, seed: 123);

        Assert.AreEqual(8, result.Iterations);
        Assert.AreEqual(11, result.FeatureCount); // 11 nodes

        // Different iterations should give different max values due to stochastic sampling
        var maxVals = result.Results.Select(r => r.MaxValue).Distinct().ToList();
        Assert.IsTrue(maxVals.Count > 1, "Stochastic variation should produce different outputs");
    }

    [TestMethod]
    public void MonteCarlo_ScenarioResult_Percentile_ReturnsCorrectLength()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(100, 300);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(10, seed: 99);

        var p50 = result.ComputePercentile(50);
        var p95 = result.ComputePercentile(95);

        Assert.AreEqual(25, p50.Length);
        Assert.AreEqual(25, p95.Length);
    }

    [TestMethod]
    public void MonteCarlo_GetFeatureDistribution_ReturnsAllIterations()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(40, 60);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(5, seed: 7);

        var dist = result.GetFeatureDistribution(0);
        Assert.AreEqual(5, dist.Length);
    }

    [TestMethod]
    public void MonteCarlo_IMonteCarloModel_Evaluate_ReturnsScalar()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(40, 60);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var rng = new RandomGenerator(42);
        double val = mc.Evaluate(rng);

        Assert.IsTrue(val > 0, "Evaluate should return positive scalar (peak value)");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 4 — Surrogate
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void SurrogateTrainer_TrainAndPredict_ProducesOutput()
    {
        // Use 3×3 grid (9 features) with 15 samples so the system is overdetermined
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 3, ny: 3)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(40, 60);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(15, seed: 42);

        var surrogate = new SurrogateTrainer(result, new Ridge());
        Assert.IsFalse(surrogate.IsTrained);

        surrogate.Train(r => r.MaxValue);
        Assert.IsTrue(surrogate.IsTrained);

        // Predict on self — should return something reasonable
        var predictions = surrogate.Predict(result.ScenarioMatrix);
        Assert.AreEqual(15, predictions.Length);
    }

    [TestMethod]
    public void SurrogateTrainer_PredictOne_ReturnsSingleValue()
    {
        // Use 3×3 grid (9 features) with 15 samples for overdetermined system
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 3, ny: 3)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(40, 60);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(15, seed: 77);

        var surrogate = new SurrogateTrainer(result, new Ridge());
        surrogate.Train(r => r.MaxValue);

        var scenario = result.GetScenarioVector(0);
        double prediction = surrogate.PredictOne(scenario);
        Assert.IsFalse(double.IsNaN(prediction));
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void SurrogateTrainer_PredictBeforeTrain_Throws()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 3, ny: 3)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation();
        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var result = mc.RunBatch(3, seed: 1);

        var surrogate = new SurrogateTrainer(result, new Ridge());
        surrogate.Predict(result.ScenarioMatrix); // Should throw
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 4 — Cluster Analyzer
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ClusterAnalyzer_HeatPlate_ProducesLabels()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(20, 80);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var mcResult = mc.RunBatch(20, seed: 42);

        var analysis = MultiphysicsClusterAnalyzer
            .For(mcResult)
            .WithAlgorithm(new KMeans { K = 3 })
            .TryClusterCounts(2, 4)
            .WithEvaluator(new SilhouetteEvaluator())
            .Run();

        Assert.IsTrue(analysis.BestClusterCount >= 2);
        Assert.IsTrue(analysis.BestClusterCount <= 4);
        Assert.AreEqual(20, analysis.Labels.Length);

        // Dominant cluster should have at least some iterations
        var dominant = analysis.GetClusterIterations(analysis.DominantCluster);
        Assert.IsTrue(dominant.Length > 0);
    }

    [TestMethod]
    public void ClusterAnalyzer_GetClusterMeanOutput_CorrectDimension()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 5, ny: 5)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(50, 300);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var mcResult = mc.RunBatch(15, seed: 55);

        var analysis = MultiphysicsClusterAnalyzer
            .For(mcResult)
            .WithAlgorithm(new KMeans { K = 2 })
            .TryClusterCounts(2, 3)
            .WithEvaluator(new InertiaEvaluator())
            .Run();

        var mean = analysis.GetClusterMeanOutput(analysis.DominantCluster);
        Assert.AreEqual(25, mean.Length); // 5×5 features
    }

    // ═══════════════════════════════════════════════════════════════
    //  Phase 5 — Analytical Validation
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void HeatPlate_SteadyState_ConvergesToLaplaceSolution()
    {
        // ∇²T = 0 with T(top)=100, T(bottom,left,right)=0.
        // Run long enough to reach steady state on a small domain.
        // For Copper on 0.1m plate: α≈117e-6, τ=L²/α≈0.1²/117e-6≈85s
        // dt=0.01 × 50000 steps = 500s ≈ 6τ → sufficient for steady state
        int nx = 21, ny = 21;
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper) // high α → fast convergence
            .WithGeometry(width: 0.1, height: 0.1, nx: nx, ny: ny)
            .WithBoundary(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.01, steps: 50000);

        // Fourier series for T(x,y) on unit square, top=T0, others=0:
        // T(x,y) = sum_{n odd} (4*T0/(n*pi)) * sin(n*pi*x/L) * sinh(n*pi*y/L)/sinh(n*pi)
        // At domain centre (0.05, 0.05) on a 0.1×0.1 plate:
        double T0 = 100.0;
        double Lx = 0.1;
        double xc = 0.05, yc = 0.05;
        double fourierSum = 0;
        for (int n = 1; n <= 99; n += 2)
        {
            double np = n * Math.PI;
            fourierSum += (4.0 * T0 / (n * Math.PI))
                        * Math.Sin(np * xc / Lx)
                        * Math.Sinh(np * yc / Lx) / Math.Sinh(np);
        }

        double simulated = result.Field[nx / 2, ny / 2];
        Assert.AreEqual(fourierSum, simulated, fourierSum * 0.10,
            $"Centre T: Fourier={fourierSum:F2}, Simulated={simulated:F2}");
    }

    [TestMethod]
    public void HeatPlate_UniformIC_NoSource_ConvergesToBoundary()
    {
        // If all BCs are 50°C and initial condition is 50°C, field should stay at 50°C
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 10, ny: 10)
            .WithBoundary(top: 50.0, bottom: 50.0, left: 50.0, right: 50.0)
            .WithInitialCondition(50.0)
            .Solve(dt: 0.0001, steps: 100);

        // Every cell should be exactly 50
        for (int ix = 0; ix < 10; ix++)
            for (int iy = 0; iy < 10; iy++)
                Assert.AreEqual(50.0, result.Field[ix, iy], 1e-10,
                    $"Field[{ix},{iy}] should be 50.0");
    }

    [TestMethod]
    public void PipeFlow_VelocityProfile_IsParabolic()
    {
        // Hagen-Poiseuille: v(r) = (R² - r²)·(-dP/dx)/(4μ)
        // Verify the ratio v(r)/v(0) ≈ 1 - (r/R)² at several radial stations
        var water = EngineeringLibrary.Water;
        double R = 0.005;
        double dPdx = -100;
        int n = 21;

        var result = SimulationType.Create(MultiphysicsType.PipeFlow)
            .WithMaterial(water)
            .WithGeometry(length: 1.0, radius: R, nodes: n)
            .WithBoundary(pressureGradient: dPdx)
            .Solve(dt: 0.02, steps: 3000);

        double vCenter = result.Values[0];
        double dr = R / (n - 1);

        // Check interior stations (skip wall node)
        for (int i = 1; i < n - 1; i++)
        {
            double r = i * dr;
            double expectedRatio = 1.0 - (r / R) * (r / R);
            double actualRatio = result.Values[i] / vCenter;
            Assert.AreEqual(expectedRatio, actualRatio, 0.15,
                $"At r={r:F4}: expected ratio {expectedRatio:F3}, got {actualRatio:F3}");
        }
    }

    [TestMethod]
    public void ElectricField_ParallelPlate_UniformExField()
    {
        // Between parallel plates at x=0 (0V) and x=L (100V),
        // the E-field should be approximately uniform: Ex ≈ -100/L
        int nx = 40, ny = 40;
        double L = 1.0;

        var result = SimulationType.Create(MultiphysicsType.ElectricField)
            .WithMaterial(EngineeringLibrary.Air)
            .WithGeometry(width: L, height: L, nx: nx, ny: ny)
            .WithBoundary(top: 50.0, bottom: 50.0, left: 0.0, right: 100.0)
            .Solve(maxIterations: 30000, tolerance: 1e-8);

        double expectedEx = -100.0 / L;

        // Sample Ex at interior points along mid-height
        int midY = ny / 2;
        for (int ix = 5; ix < nx - 5; ix++)
        {
            Assert.AreEqual(expectedEx, result.Ex[ix, midY], Math.Abs(expectedEx) * 0.15,
                $"Ex[{ix},{midY}] should be ~{expectedEx:F1}");
        }

        // Ey should be near zero in the interior
        for (int ix = 5; ix < nx - 5; ix++)
        {
            Assert.AreEqual(0.0, result.Ey[ix, midY], 5.0,
                $"Ey[{ix},{midY}] should be ~0");
        }
    }

    [TestMethod]
    public void BeamStress_FixedFixed_UniformLoad_MatchesAnalytical()
    {
        // Fixed-fixed beam with uniform load q:
        // Max deflection at midspan: δ = qL⁴/(384EI)
        double q = 1000;  // 1000 N/m
        double L = 2.0;
        double b = 0.05, h = 0.1;
        double I = b * h * h * h / 12.0;
        double E = EngineeringLibrary.Steel.YoungsModulus;
        double EI = E * I;

        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: L, nodes: 201)
            .WithCrossSection(width: b, height: h)
            .WithBoundary(BeamSupport.FixedFixed)
            .WithSource(q)
            .Solve();

        double expectedMax = q * L * L * L * L / (384.0 * EI);
        int midIdx = result.Values.Length / 2;
        double actualMid = result.Values[midIdx];
        Assert.AreEqual(expectedMax, actualMid, expectedMax * 0.05,
            $"Fixed-fixed midspan δ: expected {expectedMax:E3}, got {actualMid:E3}");
    }

    [TestMethod]
    public void BeamStress_Cantilever_MomentDistribution_Linear()
    {
        // Cantilever with tip load P:
        // Moment at distance x from fixed end: M(x) = P·(L - x)
        // M(0) = P·L, M(L) = 0
        double P = 500;
        double L = 2.0;
        int n = 101;

        var result = SimulationType.Create(MultiphysicsType.BeamStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(length: L, nodes: n)
            .WithCrossSection(width: 0.05, height: 0.1)
            .WithBoundary(BeamSupport.Cantilever)
            .AddSource(position: L, value: P)
            .Solve();

        double dx = L / (n - 1);
        // Check at 25%, 50%, 75% along the beam
        foreach (double frac in new[] { 0.25, 0.50, 0.75 })
        {
            int idx = (int)(frac * (n - 1));
            double x = idx * dx;
            double expectedM = P * (L - x);
            Assert.AreEqual(expectedM, result.BendingMoment[idx], expectedM * 0.05,
                $"M at x={x:F2}: expected {expectedM:F1}, got {result.BendingMoment[idx]:F1}");
        }
    }

    [TestMethod]
    public void ElectricField_PoissonSolver_ConvergenceBelowTolerance()
    {
        // Verify the solver converges within the specified tolerance
        // by checking that re-solving with tighter tolerance gives
        // a result very close to the first
        var result1 = SimulationType.Create(MultiphysicsType.ElectricField)
            .WithMaterial(EngineeringLibrary.Air)
            .WithGeometry(width: 1.0, height: 1.0, nx: 20, ny: 20)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 100)
            .Solve(maxIterations: 20000, tolerance: 1e-4);

        var result2 = SimulationType.Create(MultiphysicsType.ElectricField)
            .WithMaterial(EngineeringLibrary.Air)
            .WithGeometry(width: 1.0, height: 1.0, nx: 20, ny: 20)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 100)
            .Solve(maxIterations: 20000, tolerance: 1e-8);

        // Interior values should agree to ~1e-3 (the coarser tolerance)
        for (int ix = 2; ix < 18; ix++)
            for (int iy = 2; iy < 18; iy++)
                Assert.AreEqual(result2.Field[ix, iy], result1.Field[ix, iy], 0.1,
                    $"Field[{ix},{iy}] drift between tolerances too large");
    }

    [TestMethod]
    public void HeatPlate_EnergyConservation_NoSource()
    {
        // With no source and insulating-like BCs (all = 0),
        // starting from Gaussian IC, total energy should decrease
        // (heat dissipates to cold boundaries)
        int nx = 20, ny = 20;
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 1.0, height: 1.0, nx: nx, ny: ny)
            .WithBoundary(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0)
            .WithInitialCondition((x, y) => 100.0 * Math.Exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 0.02))
            .Solve(dt: 0.0001, steps: 500);

        // Total "energy" at initial state should be > final state
        double initialSum = 0, finalSum = 0;
        var initial = result.Timeline[0];
        var final2 = result.Timeline[result.Timeline.Count - 1];
        for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
            {
                initialSum += initial[ix, iy] * initial[ix, iy];
                finalSum += final2[ix, iy] * final2[ix, iy];
            }

        Assert.IsTrue(initialSum > finalSum,
            $"Energy should decrease: initial {initialSum:F2} > final {finalSum:F2}");
    }

    [TestMethod]
    public void BinaryExporter_NonSquareGrid_HeaderDimensions()
    {
        // Use a non-square HeatPlate (15×12) to verify dimension encoding
        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 1.0, height: 0.8, nx: 15, ny: 12)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 3);

        double dx = 1.0 / 15, dy = 0.8 / 12;
        var timeline = SimulationTimeline.FromResult(result, dt: 0.0001, dx: dx, dy: dy);

        string path = Path.Combine(Path.GetTempPath(), "test_nonsquare.mphy");
        try
        {
            MultiphysicsBinaryExporter.Save(timeline, path);
            var header = MultiphysicsBinaryExporter.ReadHeader(path);

            Assert.AreEqual(MultiphysicsType.HeatPlate, header.Type);
            Assert.AreEqual(15, header.Nx);
            Assert.AreEqual(12, header.Ny);
            Assert.IsTrue(header.Is2D);
            Assert.AreEqual(4, header.TimeStepCount); // initial + 3 steps
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    [TestMethod]
    public void JsonExporter_MonteCarlo_ResultContainsIterations()
    {
        var baseline = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 3, ny: 3)
            .WithBoundary(top: 100.0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(0.0);

        var variation = new ParameterVariation()
            .ThermalConductivity(40, 60);

        var mc = new MultiphysicsMonteCarloModel(baseline, variation);
        var mcResult = mc.RunBatch(5, seed: 42);

        // Export each individual result and verify it has valid JSON
        foreach (var simResult in mcResult.Results)
        {
            string json = MultiphysicsJsonExporter.ToJson(simResult);
            Assert.IsTrue(json.StartsWith("{"));
            Assert.IsTrue(json.Contains("\"type\":\"HeatPlate\""));
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  CylinderFlow — Basic Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CylinderFlow_FluentAPI_ProducesResult()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.5, height: 0.2, nx: 50, ny: 20)
            .WithCylinder(centerX: 0.1, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 20);

        Assert.AreEqual(MultiphysicsType.CylinderFlow, result.Type);
        Assert.IsNotNull(result.Vx);
        Assert.IsNotNull(result.Vy);
        Assert.IsNotNull(result.Pressure);
        Assert.IsNotNull(result.Vorticity);
        Assert.IsNotNull(result.CylinderMask);
        Assert.AreEqual(50, result.Vx.GetLength(0));
        Assert.AreEqual(20, result.Vx.GetLength(1));
    }

    [TestMethod]
    public void CylinderFlow_CylinderMask_IsCorrect()
    {
        int nx = 40, ny = 20;
        double width = 0.4, height = 0.2;
        double cx = 0.1, cy = 0.1, cr = 0.02;

        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: width, height: height, nx: nx, ny: ny)
            .WithCylinder(centerX: cx, centerY: cy, radius: cr)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 5);

        // Verify cells inside cylinder are masked
        double dx = width / nx, dy = height / ny;
        int maskedCount = 0;
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                double x = ix * dx, y = iy * dy;
                double dist = Math.Sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                if (dist <= cr)
                {
                    Assert.IsTrue(result.CylinderMask[ix, iy],
                        $"Cell ({ix},{iy}) at dist {dist:F4} should be masked");
                    maskedCount++;
                }
            }
        Assert.IsTrue(maskedCount > 0, "Should have some masked cells");
    }

    [TestMethod]
    public void CylinderFlow_NoSlip_VelocityZeroInsideCylinder()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.4, height: 0.2, nx: 40, ny: 20)
            .WithCylinder(centerX: 0.1, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 20);

        // All cells inside the cylinder should have zero velocity
        int nx = result.Vx.GetLength(0), ny = result.Vx.GetLength(1);
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                if (result.CylinderMask[ix, iy])
                {
                    Assert.AreEqual(0.0, result.Vx[ix, iy], 1e-12,
                        $"Vx inside cylinder at ({ix},{iy}) should be 0");
                    Assert.AreEqual(0.0, result.Vy[ix, iy], 1e-12,
                        $"Vy inside cylinder at ({ix},{iy}) should be 0");
                }
            }
    }

    [TestMethod]
    public void CylinderFlow_InletBC_MaintainsVelocity()
    {
        double uInlet = 0.05;
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.5, height: 0.2, nx: 50, ny: 20)
            .WithCylinder(centerX: 0.15, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(uInlet)
            .Solve(dt: 0.0005, steps: 30);

        // Left column (inlet) should have Vx = uInlet, Vy = 0
        int ny = result.Vx.GetLength(1);
        for (int iy = 0; iy < ny; iy++)
        {
            Assert.AreEqual(uInlet, result.Vx[0, iy], 1e-10,
                $"Inlet Vx at iy={iy} should be {uInlet}");
            Assert.AreEqual(0.0, result.Vy[0, iy], 1e-10,
                $"Inlet Vy at iy={iy} should be 0");
        }
    }

    [TestMethod]
    public void CylinderFlow_MassConservation_DivergenceNearZero()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.4, height: 0.2, nx: 40, ny: 20)
            .WithCylinder(centerX: 0.1, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 30);

        // Reconstruct divergence and check it's small at interior fluid cells
        int nx = 40, ny = 20;
        double dx = 0.4 / nx, dy = 0.2 / ny;
        var grid = new Grid2D(nx, ny, dx, dy);

        var vxFlat = new double[nx * ny];
        var vyFlat = new double[nx * ny];
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++)
            {
                int idx = grid.Index(ix, iy);
                vxFlat[idx] = result.Vx[ix, iy];
                vyFlat[idx] = result.Vy[ix, iy];
            }

        var div = GridOperators.Divergence2D(
            new VectorN(vxFlat), new VectorN(vyFlat), grid, BoundaryCondition.Dirichlet);

        // Check divergence at interior fluid cells (exclude boundaries and cylinder)
        double maxDiv = 0;
        for (int iy = 2; iy < ny - 2; iy++)
            for (int ix = 2; ix < nx - 2; ix++)
            {
                if (result.CylinderMask[ix, iy]) continue;
                double d = Math.Abs(div[grid.Index(ix, iy)]);
                if (d > maxDiv) maxDiv = d;
            }

        Assert.IsTrue(maxDiv < 1.0,
            $"Max divergence at interior fluid cells should be small, got {maxDiv:E3}");
    }

    [TestMethod]
    public void CylinderFlow_Timeline_HasSnapshots()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.4, height: 0.2, nx: 30, ny: 15)
            .WithCylinder(centerX: 0.1, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 25);

        Assert.IsNotNull(result.Timeline);
        Assert.IsTrue(result.Timeline.Count > 0, "Should have vorticity timeline snapshots");
        Assert.AreEqual(30, result.Timeline[0].GetLength(0));
        Assert.AreEqual(15, result.Timeline[0].GetLength(1));
    }

    [TestMethod]
    public void CylinderFlow_DragCoefficient_IsPositive()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.5, height: 0.2, nx: 50, ny: 20)
            .WithCylinder(centerX: 0.12, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(0.02)
            .Solve(dt: 0.0005, steps: 100);

        // After some steps, drag coefficient should be positive
        Assert.IsTrue(result.DragCoefficient > 0,
            $"Drag coefficient should be positive, got {result.DragCoefficient:E3}");
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void CylinderFlow_WithoutCylinder_Throws()
    {
        SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.5, height: 0.2, nx: 50, ny: 20)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 10);
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void CylinderFlow_WithoutInletVelocity_Throws()
    {
        SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.5, height: 0.2, nx: 50, ny: 20)
            .WithCylinder(centerX: 0.1, centerY: 0.1, radius: 0.02)
            .Solve(dt: 0.001, steps: 10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  CylinderFlow — Export
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CylinderFlow_JsonExport_ContainsVelocityFields()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.3, height: 0.15, nx: 20, ny: 10)
            .WithCylinder(centerX: 0.07, centerY: 0.075, radius: 0.015)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 10);

        string json = MultiphysicsJsonExporter.ToJson(result);
        Assert.IsTrue(json.Contains("\"type\":\"CylinderFlow\""));
        Assert.IsTrue(json.Contains("\"vx\":["));
        Assert.IsTrue(json.Contains("\"vy\":["));
        Assert.IsTrue(json.Contains("\"pressure\":["));
        Assert.IsTrue(json.Contains("\"vorticity\":["));
        Assert.IsTrue(json.Contains("\"dragCoefficient\""));
        Assert.IsTrue(json.Contains("\"strouhalNumber\""));
    }

    [TestMethod]
    public void CylinderFlow_BinaryExport_RoundTrips()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.3, height: 0.15, nx: 20, ny: 10)
            .WithCylinder(centerX: 0.07, centerY: 0.075, radius: 0.015)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 10);

        string path = Path.Combine(Path.GetTempPath(), "test_cylinder.mphy");
        try
        {
            MultiphysicsBinaryExporter.Save(result, path);
            Assert.IsTrue(File.Exists(path));

            var header = MultiphysicsBinaryExporter.ReadHeader(path);
            Assert.AreEqual(MultiphysicsType.CylinderFlow, header.Type);
            Assert.AreEqual(20, header.Nx);
            Assert.AreEqual(10, header.Ny);
            Assert.AreEqual(4, header.LayerCount);
            Assert.IsTrue(header.Is2D);

            var data = MultiphysicsBinaryExporter.Read(path);
            Assert.AreEqual(4, data.Layers.Length); // vx, vy, p, vorticity
            Assert.AreEqual(200, data.Layers[0].Length); // 20×10
        }
        finally
        {
            if (File.Exists(path)) File.Delete(path);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — HeatBlock3D
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void HeatBlock3D_FluentAPI_ProducesResult()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 10, ny: 10, nz: 10)
            .WithBoundary3D(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 20);

        Assert.AreEqual(MultiphysicsType.HeatBlock3D, result.Type);
        Assert.IsNotNull(result.Field3D);
        Assert.IsNotNull(result.Timeline3D);
        Assert.AreEqual(21, result.Timeline3D.Count); // initial + 20 steps
        Assert.AreEqual(100.0, result.MaxValue, 0.1);
    }

    [TestMethod]
    public void HeatBlock3D_BoundaryConditions_AreApplied()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: 8, ny: 8, nz: 8)
            .WithBoundary3D(top: 200.0, bottom: 50.0, left: 50.0, right: 50.0, front: 50.0, back: 50.0)
            .WithInitialCondition(50.0)
            .Solve(dt: 0.001, steps: 5);

        // Top face interior point should be 200
        double topMid = result.Field3D[4, 7, 4]; // ix=4, iy=7 (top row), iz=4
        Assert.AreEqual(200.0, topMid, 0.01);

        // Bottom face interior point should be 50
        double botMid = result.Field3D[4, 0, 4];
        Assert.AreEqual(50.0, botMid, 0.01);
    }

    [TestMethod]
    public void HeatBlock3D_WithSource_IncreasesTemperature()
    {
        var noSource = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 10, ny: 10, nz: 10)
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 20);

        var withSource = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 10, ny: 10, nz: 10)
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(0.0)
            .AddSource3D(5, 5, 5, 1e6)
            .Solve(dt: 0.0001, steps: 20);

        Assert.IsTrue(withSource.MaxValue > noSource.MaxValue,
            "Source should increase temperature");
    }

    [TestMethod]
    public void HeatBlock3D_UniformIC_BoundaryDriven_Decays()
    {
        // Start at 100, cool to 0 on all faces — interior should decrease
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 10, ny: 10, nz: 10)
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(100.0)
            .Solve(dt: 0.0001, steps: 50);

        // Centre should have cooled from 100 towards 0
        double centre = result.Field3D[5, 5, 5];
        Assert.IsTrue(centre < 100.0, $"Centre should cool, got {centre}");
        Assert.IsTrue(centre > 0.0, $"Centre should still be above 0, got {centre}");
    }

    [TestMethod]
    public void HeatBlock3D_InitialCondition3D_IsApplied()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: 8, ny: 8, nz: 8)
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition3D((x, y, z) => 50.0)
            .Solve(dt: 0.0001, steps: 1);

        // First timeline entry (initial) should have 50 at interior
        double[,,] initial = result.Timeline3D[0];
        double interior = initial[4, 4, 4];
        Assert.AreEqual(50.0, interior, 1.0);
    }

    [TestMethod]
    public void HeatBlock3D_SliceXY_MatchesField()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 8, ny: 8, nz: 8)
            .WithBoundary3D(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(50.0)
            .Solve(dt: 0.0001, steps: 10);

        int iz = 4;
        var slice = result.SliceXY(iz);
        Assert.AreEqual(8, slice.GetLength(0));
        Assert.AreEqual(8, slice.GetLength(1));

        for (int iy = 0; iy < 8; iy++)
            for (int ix = 0; ix < 8; ix++)
                Assert.AreEqual(result.Field3D[ix, iy, iz], slice[ix, iy], 1e-12);
    }

    [TestMethod]
    public void HeatBlock3D_SliceXZ_MatchesField()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 8, ny: 8, nz: 8)
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 100.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(20.0)
            .Solve(dt: 0.0001, steps: 5);

        int iy = 3;
        var slice = result.SliceXZ(iy);
        Assert.AreEqual(8, slice.GetLength(0));
        Assert.AreEqual(8, slice.GetLength(1));

        for (int iz = 0; iz < 8; iz++)
            for (int ix = 0; ix < 8; ix++)
                Assert.AreEqual(result.Field3D[ix, iy, iz], slice[ix, iz], 1e-12);
    }

    [TestMethod]
    public void HeatBlock3D_SliceYZ_MatchesField()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 8, ny: 8, nz: 8)
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0, front: 100.0, back: 0.0)
            .WithInitialCondition(20.0)
            .Solve(dt: 0.0001, steps: 5);

        int ix = 3;
        var slice = result.SliceYZ(ix);
        Assert.AreEqual(8, slice.GetLength(0));
        Assert.AreEqual(8, slice.GetLength(1));

        for (int iz = 0; iz < 8; iz++)
            for (int iy = 0; iy < 8; iy++)
                Assert.AreEqual(result.Field3D[ix, iy, iz], slice[iy, iz], 1e-12);
    }

    [TestMethod]
    public void HeatBlock3D_Timeline_GrowsWithSteps()
    {
        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 6, ny: 6, nz: 6)
            .WithBoundary3D(top: 100.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 30);

        Assert.AreEqual(31, result.Timeline3D.Count);
        Assert.AreEqual(6, result.Timeline3D[0].GetLength(0));
        Assert.AreEqual(6, result.Timeline3D[0].GetLength(1));
        Assert.AreEqual(6, result.Timeline3D[0].GetLength(2));
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void HeatBlock3D_WithoutNz_Throws()
    {
        SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 10, ny: 10) // 2D geometry — missing Nz
            .WithBoundary3D(top: 0.0, bottom: 0.0, left: 0.0, right: 0.0, front: 0.0, back: 0.0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.0001, steps: 5);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — FluidDiffusion3D
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidDiffusion3D_PureDiffusion_SphericalSpread()
    {
        // Zero velocity, point source IC at centre → should spread radially
        int n = 12;
        var result = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithInitialCondition3D((x, y, z) =>
            {
                double cx = 0.5, cy = 0.5, cz = 0.5;
                double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz);
                return r2 < 0.04 ? 100.0 : 0.0;
            })
            .Solve(dt: 0.001, steps: 50);

        Assert.AreEqual(MultiphysicsType.FluidDiffusion3D, result.Type);
        Assert.IsNotNull(result.Field3D);

        // Centre should have decreased from initial
        double centre = result.Field3D[n / 2, n / 2, n / 2];
        Assert.IsTrue(centre < 100.0, $"Centre should diffuse away, got {centre}");
        Assert.IsTrue(centre > 0.0, $"Centre should still be positive, got {centre}");
    }

    [TestMethod]
    public void FluidDiffusion3D_UniformAdvection_ShiftsPlume()
    {
        // Uniform velocity in +x should shift mass centroid downstream
        // compared to pure diffusion (zero velocity).
        int n = 12;
        double vxVal = 1.0;

        Func<double, double, double, double> blob = (x, y, z) =>
        {
            double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5);
            return r2 < 0.04 ? 50.0 : 0.0;
        };

        var noDrift = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithInitialCondition3D(blob)
            .Solve(dt: 0.005, steps: 40);

        var withDrift = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithVelocityField3D((x, y, z) => (vxVal, 0, 0))
            .WithInitialCondition3D(blob)
            .Solve(dt: 0.005, steps: 40);

        // Compute x-centroid for both: ΣcᵢXᵢ / Σcᵢ
        double centroidNoDrift = XCentroid(noDrift.Field3D, n, 1.0 / n);
        double centroidDrift = XCentroid(withDrift.Field3D, n, 1.0 / n);

        Assert.IsTrue(centroidDrift > centroidNoDrift,
            $"Advected centroid ({centroidDrift:F3}) should be downstream of pure diffusion ({centroidNoDrift:F3})");
    }

    private static double XCentroid(double[,,] field, int n, double dx)
    {
        double sumCX = 0, sumC = 0;
        for (int iz = 0; iz < n; iz++)
            for (int iy = 0; iy < n; iy++)
                for (int ix = 0; ix < n; ix++)
                {
                    double c = field[ix, iy, iz];
                    double x = ix * dx;
                    sumCX += c * x;
                    sumC += c;
                }
        return sumC > 0 ? sumCX / sumC : 0;
    }

    [TestMethod]
    public void FluidDiffusion3D_WithSource_IncreasesConcentration()
    {
        int n = 10;
        var noSource = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.5, height: 0.5, depth: 0.5, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.001, steps: 20);

        var withSource = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.5, height: 0.5, depth: 0.5, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithInitialCondition(0.0)
            .AddSource3D(5, 5, 5, 100.0)
            .Solve(dt: 0.001, steps: 20);

        Assert.IsTrue(withSource.MaxValue > noSource.MaxValue,
            "Source should increase concentration");
    }

    [TestMethod]
    public void FluidDiffusion3D_CFL_Warning()
    {
        // Very large velocity should trigger CFL flag
        int n = 8;
        var result = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithVelocityField3D((x, y, z) => (100.0, 0, 0)) // very fast
            .WithInitialCondition(0.0)
            .Solve(dt: 0.1, steps: 1); // large dt

        Assert.IsTrue(result.CflClamped, "CFL flag should be set for large velocity");
    }

    [TestMethod]
    public void FluidDiffusion3D_Timeline_GrowsWithSteps()
    {
        int n = 6;
        var result = SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.5, height: 0.5, depth: 0.5, nx: n, ny: n, nz: n)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithDiffusionCoefficient(0.01)
            .WithInitialCondition(10.0)
            .Solve(dt: 0.001, steps: 15);

        Assert.AreEqual(16, result.Timeline3D.Count); // initial + 15
        Assert.AreEqual(n, result.Timeline3D[0].GetLength(0));
        Assert.AreEqual(n, result.Timeline3D[0].GetLength(1));
        Assert.AreEqual(n, result.Timeline3D[0].GetLength(2));
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void FluidDiffusion3D_WithoutDiffusionCoeff_Throws()
    {
        SimulationType.Create(MultiphysicsType.FluidDiffusion3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 1.0, height: 1.0, depth: 1.0, nx: 8, ny: 8, nz: 8)
            .WithBoundary3D(top: 0, bottom: 0, left: 0, right: 0, front: 0, back: 0)
            .WithInitialCondition(0.0)
            .Solve(dt: 0.001, steps: 5);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluent API — CylinderFlow3D
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void CylinderFlow3D_FluentAPI_ProducesResult()
    {
        int nx = 30, ny = 15, nz = 6;
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.3, height: 0.15, depth: 0.06, nx: nx, ny: ny, nz: nz)
            .WithCylinder(centerX: 0.08, centerY: 0.075, radius: 0.015)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 5);

        Assert.AreEqual(MultiphysicsType.CylinderFlow3D, result.Type);
        Assert.IsNotNull(result.Field3D);
        Assert.IsNotNull(result.Vx3D);
        Assert.IsNotNull(result.Vy3D);
        Assert.IsNotNull(result.Vz3D);
        Assert.IsNotNull(result.Pressure3D);
        Assert.IsNotNull(result.CylinderMask3D);

        Assert.AreEqual(nx, result.Vx3D.GetLength(0));
        Assert.AreEqual(ny, result.Vx3D.GetLength(1));
        Assert.AreEqual(nz, result.Vx3D.GetLength(2));
    }

    [TestMethod]
    public void CylinderFlow3D_NoSlip_VelocityZeroInsideCylinder()
    {
        int nx = 30, ny = 15, nz = 4;
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.3, height: 0.15, depth: 0.04, nx: nx, ny: ny, nz: nz)
            .WithCylinder(centerX: 0.08, centerY: 0.075, radius: 0.015)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 10);

        // All cells inside cylinder should have zero velocity
        for (int iz = 0; iz < nz; iz++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                {
                    if (result.CylinderMask3D[ix, iy])
                    {
                        Assert.AreEqual(0.0, result.Vx3D[ix, iy, iz], 1e-12,
                            $"Vx inside cylinder at ({ix},{iy},{iz})");
                        Assert.AreEqual(0.0, result.Vy3D[ix, iy, iz], 1e-12,
                            $"Vy inside cylinder at ({ix},{iy},{iz})");
                        Assert.AreEqual(0.0, result.Vz3D[ix, iy, iz], 1e-12,
                            $"Vz inside cylinder at ({ix},{iy},{iz})");
                    }
                }
    }

    [TestMethod]
    public void CylinderFlow3D_Drag_IsReasonable()
    {
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.5, height: 0.2, depth: 0.04, nx: 40, ny: 16, nz: 4)
            .WithCylinder(centerX: 0.12, centerY: 0.1, radius: 0.02)
            .WithInletVelocity(0.02)
            .Solve(dt: 0.0005, steps: 50);

        // After some steps, drag coefficient should have non-trivial magnitude
        double absCd = Math.Abs(result.DragCoefficient);
        Assert.IsTrue(absCd > 1e-6,
            $"Drag should have non-trivial magnitude, got {result.DragCoefficient:E3}");
    }

    [TestMethod]
    public void CylinderFlow3D_AxialSymmetry_ZSlicesAreSimilar()
    {
        // With periodic z BCs and uniform initial conditions, all z-slices should be similar
        int nx = 20, ny = 10, nz = 6;
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.2, height: 0.1, depth: 0.06, nx: nx, ny: ny, nz: nz)
            .WithCylinder(centerX: 0.05, centerY: 0.05, radius: 0.01)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 10);

        // Compare z-slices of vx: take interior z-slices (1 and nz-2)
        double maxDiff = 0;
        for (int iy = 1; iy < ny - 1; iy++)
            for (int ix = 1; ix < nx - 1; ix++)
            {
                double diff = Math.Abs(result.Vx3D[ix, iy, 1] - result.Vx3D[ix, iy, nz - 2]);
                if (diff > maxDiff) maxDiff = diff;
            }

        Assert.IsTrue(maxDiff < 0.01,
            $"Z-slices should be similar (periodic z), max diff = {maxDiff:E3}");
    }

    [TestMethod]
    public void CylinderFlow3D_VelocityMagnitude_PositiveEverywhere()
    {
        int nx = 20, ny = 10, nz = 4;
        var result = SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.2, height: 0.1, depth: 0.04, nx: nx, ny: ny, nz: nz)
            .WithCylinder(centerX: 0.05, centerY: 0.05, radius: 0.01)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 5);

        // Field3D is velocity magnitude, should be ≥ 0 everywhere
        for (int iz = 0; iz < nz; iz++)
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                    Assert.IsTrue(result.Field3D[ix, iy, iz] >= 0,
                        $"Velocity magnitude at ({ix},{iy},{iz}) should be ≥ 0");
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void CylinderFlow3D_WithoutCylinder_Throws()
    {
        SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.3, height: 0.15, depth: 0.06, nx: 20, ny: 10, nz: 4)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 5);
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void CylinderFlow3D_WithoutInletVelocity_Throws()
    {
        SimulationType.Create(MultiphysicsType.CylinderFlow3D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry3D(width: 0.3, height: 0.15, depth: 0.06, nx: 20, ny: 10, nz: 4)
            .WithCylinder(centerX: 0.08, centerY: 0.075, radius: 0.015)
            .Solve(dt: 0.001, steps: 5);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Convection (Robin) Boundary Conditions
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void HeatPlate_ConvectionBC_CoolsToAmbient()
    {
        // Hot plate (100 K) cooling to ambient (20 K) via convection on all faces.
        // Use a small domain so diffusion reaches the centre quickly.
        double ambient = 20.0;
        double h = 50.0; // W/(m²·K)

        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.02, height: 0.02, nx: 10, ny: 10)
            .WithConvectionBoundary(h, ambient)
            .WithInitialCondition(100.0)
            .Solve(dt: 0.00001, steps: 2000);

        Assert.IsNotNull(result.Field);

        // Interior temperature should have dropped from 100 towards ambient
        double centerT = result.Field[5, 5];
        Assert.IsTrue(centerT < 100.0, $"Centre ({centerT:F2}) should cool from 100 K");
        Assert.IsTrue(centerT > ambient, $"Centre ({centerT:F2}) should still be above ambient");
    }

    [TestMethod]
    public void HeatPlate_ConvectionBC_ConvergesLowerThanDirichlet()
    {
        // Convection BC with T∞=0 should cool slower than Dirichlet T=0
        // because convection provides finite resistance while Dirichlet instantly clamps.
        double h = 10.0;

        var dirichlet = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 0.05, height: 0.05, nx: 10, ny: 10)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 0)
            .WithInitialCondition(100.0)
            .Solve(dt: 0.00005, steps: 200);

        var convection = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 0.05, height: 0.05, nx: 10, ny: 10)
            .WithConvectionBoundary(h, 0.0)
            .WithInitialCondition(100.0)
            .Solve(dt: 0.00005, steps: 200);

        double dirichletCenter = dirichlet.Field[5, 5];
        double convectionCenter = convection.Field[5, 5];

        // Convection should retain more heat (cool slower)
        Assert.IsTrue(convectionCenter > dirichletCenter,
            $"Convection centre ({convectionCenter:F2}) should be warmer than Dirichlet ({dirichletCenter:F2})");
    }

    [TestMethod]
    public void HeatPlate_MixedBC_DirichletAndConvection()
    {
        // Top = Dirichlet 200 K, other faces = convection to 20 K.
        // Use small domain + aluminum (high α) so heat reaches centre quickly.
        double ambient = 20.0;
        double h = 25.0;

        var result = SimulationType.Create(MultiphysicsType.HeatPlate)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.02, height: 0.02, nx: 10, ny: 10)
            .WithConvectionBoundary(h, ambient)
            .WithFaceBoundary("top", FaceBoundaryCondition.Dirichlet(200.0))
            .WithInitialCondition(ambient)
            .Solve(dt: 0.000005, steps: 3000);

        // Top row should be fixed at 200
        double topMid = result.Field[5, 9];
        Assert.AreEqual(200.0, topMid, 0.01);

        // Interior should have heated above ambient due to top heating
        double center = result.Field[5, 5];
        Assert.IsTrue(center > ambient, $"Centre ({center:F2}) should heat up from top BC");
    }

    [TestMethod]
    public void HeatBlock3D_ConvectionBC_CoolsToAmbient()
    {
        // Hot block cooling convectively on all six faces
        double ambient = 25.0;
        double h = 40.0;

        var result = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry3D(width: 0.1, height: 0.1, depth: 0.1, nx: 8, ny: 8, nz: 8)
            .WithConvectionBoundary3D(h, ambient)
            .WithInitialCondition(200.0)
            .Solve(dt: 0.00005, steps: 300);

        Assert.IsNotNull(result.Field3D);

        double centerT = result.Field3D[4, 4, 4];
        Assert.IsTrue(centerT < 200.0, "Centre should cool from 200 K");
        Assert.IsTrue(centerT > ambient, "Centre should still be above ambient");
    }

    [TestMethod]
    public void HeatBlock3D_ConvectionBC_SlowerThanDirichlet()
    {
        double h = 10.0;
        int n = 6;

        var dirichlet = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry3D(width: 0.05, height: 0.05, depth: 0.05, nx: n, ny: n, nz: n)
            .WithBoundary3D(0, 0, 0, 0, 0, 0)
            .WithInitialCondition(100.0)
            .Solve(dt: 0.00005, steps: 200);

        var convection = SimulationType.Create(MultiphysicsType.HeatBlock3D)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry3D(width: 0.05, height: 0.05, depth: 0.05, nx: n, ny: n, nz: n)
            .WithConvectionBoundary3D(h, 0.0)
            .WithInitialCondition(100.0)
            .Solve(dt: 0.00005, steps: 200);

        double dirichletCenter = dirichlet.Field3D[n / 2, n / 2, n / 2];
        double convectionCenter = convection.Field3D[n / 2, n / 2, n / 2];

        Assert.IsTrue(convectionCenter > dirichletCenter,
            $"Convection centre ({convectionCenter:F2}) should be warmer than Dirichlet ({dirichletCenter:F2})");
    }

    [TestMethod]
    public void ConvectiveBoundaryRate_CorrectSign()
    {
        // When T_cell > T_ambient, the rate should be negative (cooling)
        var model = new CSharpNumerics.Physics.Thermodynamics.HeatTransferModel();
        double rate = model.ConvectiveBoundaryRate(
            h: 50, density: 2700, specificHeat: 900, dx: 0.01,
            cellTemperature: 100, ambientTemperature: 20);

        Assert.IsTrue(rate < 0, "Rate should be negative when cell is hotter than ambient");
    }

    [TestMethod]
    public void ConvectiveBoundaryRate_ZeroWhenEqual()
    {
        var model = new CSharpNumerics.Physics.Thermodynamics.HeatTransferModel();
        double rate = model.ConvectiveBoundaryRate(
            h: 50, density: 2700, specificHeat: 900, dx: 0.01,
            cellTemperature: 20, ambientTemperature: 20);

        Assert.AreEqual(0.0, rate, 1e-15);
    }

    [TestMethod]
    public void FaceBoundaryCondition_DirichletFactory()
    {
        var bc = FaceBoundaryCondition.Dirichlet(300.0);
        Assert.AreEqual(FaceBCType.Dirichlet, bc.Type);
        Assert.AreEqual(300.0, bc.Temperature, 1e-10);
        Assert.AreEqual(0.0, bc.HeatTransferCoefficient, 1e-10);
    }

    [TestMethod]
    public void FaceBoundaryCondition_ConvectionFactory()
    {
        var bc = FaceBoundaryCondition.Convection(heatTransferCoefficient: 25, ambientTemperature: 20);
        Assert.AreEqual(FaceBCType.Convection, bc.Type);
        Assert.AreEqual(20.0, bc.Temperature, 1e-10);
        Assert.AreEqual(25.0, bc.HeatTransferCoefficient, 1e-10);
    }
}
