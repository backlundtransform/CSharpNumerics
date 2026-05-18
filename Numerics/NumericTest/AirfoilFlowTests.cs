using CSharpNumerics.Engines.Multiphysics;
using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;
using CSharpNumerics.Physics.Materials.Engineering;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest;

[TestClass]
public class AirfoilFlowTests
{
    // ═══════════════════════════════════════════════════════════════
    //  NACA Geometry Tests
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void NACAGeometry_Symmetric0012_GeneratesClosedContour()
    {
        var (x, y) = NACAGeometry.Generate("0012", 60);

        // Should have 61 points (60 panels + 1 to close)
        Assert.AreEqual(61, x.Length);
        Assert.AreEqual(61, y.Length);

        // First and last points should be at the trailing edge (approximately same location)
        Assert.AreEqual(x[0], x[60], 0.02); // open TE for NACA 4-digit
        Assert.AreEqual(y[0], y[60], 0.02);
    }

    [TestMethod]
    public void NACAGeometry_Symmetric_IsSymmetricAboutChordLine()
    {
        var (x, y) = NACAGeometry.Generate("0012", 100);

        // For symmetric airfoil, upper and lower surfaces should be mirror images
        // Upper surface: indices 0 to 50, lower: indices 50 to 100
        for (int i = 1; i < 50; i++)
        {
            int upperIdx = i;
            int lowerIdx = 100 - i;
            Assert.AreEqual(x[upperIdx], x[lowerIdx], 1e-10,
                $"X-coordinates should match at index {i}");
            Assert.AreEqual(y[upperIdx], -y[lowerIdx], 1e-10,
                $"Y-coordinates should be mirrored at index {i}");
        }
    }

    [TestMethod]
    public void NACAGeometry_Cambered2412_HasPositiveCamber()
    {
        var (x, y) = NACAGeometry.Generate("2412", 100);

        // Camber line should be above zero for a 2412 airfoil
        // The midpoint of upper and lower should be positive (positive camber)
        double midY = (y[25] + y[75]) / 2.0; // roughly at 25% chord
        Assert.IsTrue(midY > 0, "NACA 2412 should have positive camber");
    }

    [TestMethod]
    public void NACAGeometry_ChordScaling_Works()
    {
        var (x1, _) = NACAGeometry.Generate("0012", 40, 1.0);
        var (x2, _) = NACAGeometry.Generate("0012", 40, 2.0);

        // Chord=2 should produce coordinates scaled by factor 2
        for (int i = 0; i < x1.Length; i++)
            Assert.AreEqual(x1[i] * 2.0, x2[i], 1e-12);
    }

    [TestMethod]
    public void NACAGeometry_Rotate_ZeroAlpha_NoChange()
    {
        var (x, y) = NACAGeometry.Generate("0012", 40);
        var (xr, yr) = NACAGeometry.Rotate(x, y, 0.0);

        for (int i = 0; i < x.Length; i++)
        {
            Assert.AreEqual(x[i], xr[i], 1e-12);
            Assert.AreEqual(y[i], yr[i], 1e-12);
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Panel Method Tests
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PanelMethod_SymmetricAirfoil_ZeroAlpha_ZeroLift()
    {
        var (x, y) = NACAGeometry.Generate("0012", 80);
        var result = PanelMethod.Solve(x, y, 0.0);

        // Symmetric airfoil at zero AoA should have Cl ≈ 0
        Assert.AreEqual(0.0, result.Cl, 0.05, "Cl should be ~0 for symmetric airfoil at α=0");
        Assert.AreEqual(0.0, result.Cm, 0.05, "Cm should be ~0 for symmetric airfoil at α=0");
    }

    [TestMethod]
    public void PanelMethod_SymmetricAirfoil_PositiveAlpha_PositiveLift()
    {
        var (x, y) = NACAGeometry.Generate("0012", 100);
        double alpha = 5.0 * Math.PI / 180.0; // 5 degrees
        var (xr, yr) = NACAGeometry.Rotate(x, y, alpha);
        var result = PanelMethod.Solve(xr, yr, alpha);

        // At 5° AoA, thin airfoil theory predicts Cl ≈ 2π·α ≈ 0.55
        Assert.IsTrue(result.Cl > 0.3, $"Cl should be positive at α=5°, got {result.Cl}");
        Assert.IsTrue(result.Cl < 1.0, $"Cl should be reasonable at α=5°, got {result.Cl}");
    }

    [TestMethod]
    public void PanelMethod_LiftIncreasesWithAlpha()
    {
        var (x, y) = NACAGeometry.Generate("0012", 80);

        double alpha1 = 2.0 * Math.PI / 180.0;
        double alpha2 = 6.0 * Math.PI / 180.0;

        var (xr1, yr1) = NACAGeometry.Rotate(x, y, alpha1);
        var (xr2, yr2) = NACAGeometry.Rotate(x, y, alpha2);

        var result1 = PanelMethod.Solve(xr1, yr1, alpha1);
        var result2 = PanelMethod.Solve(xr2, yr2, alpha2);

        Assert.IsTrue(result2.Cl > result1.Cl,
            $"Cl at 6° ({result2.Cl}) should exceed Cl at 2° ({result1.Cl})");
    }

    [TestMethod]
    public void PanelMethod_CpDistribution_HasStagnationPoint()
    {
        var (x, y) = NACAGeometry.Generate("0012", 100);
        var result = PanelMethod.Solve(x, y, 0.0);

        // At stagnation point, Cp = 1.0
        double maxCp = double.MinValue;
        for (int i = 0; i < result.Cp.Length; i++)
            if (result.Cp[i] > maxCp) maxCp = result.Cp[i];

        Assert.IsTrue(maxCp > 0.8, $"Max Cp should be near 1.0 (stagnation), got {maxCp}");
    }

    [TestMethod]
    public void PanelMethod_CpDistribution_HasSuction()
    {
        var (x, y) = NACAGeometry.Generate("0012", 100);
        double alpha = 5.0 * Math.PI / 180.0;
        var (xr, yr) = NACAGeometry.Rotate(x, y, alpha);
        var result = PanelMethod.Solve(xr, yr, alpha);

        // At positive alpha, upper surface should have Cp < 0 (suction)
        double minCp = double.MaxValue;
        for (int i = 0; i < result.Cp.Length; i++)
            if (result.Cp[i] < minCp) minCp = result.Cp[i];

        Assert.IsTrue(minCp < -0.5, $"Min Cp should be well below 0 (suction peak), got {minCp}");
    }

    [TestMethod]
    public void PanelMethod_VelocityField_FarawayMatchesFreestream()
    {
        var (x, y) = NACAGeometry.Generate("0012", 60);
        var result = PanelMethod.Solve(x, y, 0.0, 1.0);

        // Query point far upstream should have velocity ≈ freestream
        var queryX = new double[] { -10.0 };
        var queryY = new double[] { 0.0 };

        var (u, v) = PanelMethod.VelocityField(result, x, y, queryX, queryY, 1.0);

        Assert.AreEqual(1.0, u[0], 0.1, "u far upstream should ≈ freestream");
        Assert.AreEqual(0.0, v[0], 0.1, "v far upstream should ≈ 0");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Multiphysics Engine Integration Tests
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void AirfoilFlow_MultiphysicsScenario_BasicSolve()
    {
        var result = SimulationType.Create(MultiphysicsType.AirfoilFlow)
            .WithMaterial(EngineeringLibrary.Air)
            .WithAirfoil("0012")
            .WithAngleOfAttack(5.0 * Math.PI / 180.0)
            .WithFreestream(30.0) // 30 m/s
            .Solve();

        Assert.AreEqual(MultiphysicsType.AirfoilFlow, result.Type);
        Assert.IsNotNull(result.CpDistribution);
        Assert.IsNotNull(result.SurfaceVelocity);
        Assert.IsNotNull(result.AirfoilX);
        Assert.IsNotNull(result.AirfoilY);
        Assert.IsTrue(result.LiftCoefficient > 0, "Lift should be positive at α=5°");
    }

    [TestMethod]
    public void AirfoilFlow_WithVelocityField_ProducesGridData()
    {
        var result = SimulationType.Create(MultiphysicsType.AirfoilFlow)
            .WithMaterial(EngineeringLibrary.Air)
            .WithAirfoil("2412")
            .WithAngleOfAttack(3.0 * Math.PI / 180.0)
            .WithFreestream(20.0)
            .WithGeometry(4.0, 3.0, 40, 30) // 4m × 3m domain, 40×30 grid
            .Solve();

        Assert.IsNotNull(result.Vx, "Should produce Vx field");
        Assert.IsNotNull(result.Vy, "Should produce Vy field");
        Assert.IsNotNull(result.Pressure, "Should produce pressure field");
        Assert.AreEqual(40, result.Vx.GetLength(0));
        Assert.AreEqual(30, result.Vx.GetLength(1));
    }

    [TestMethod]
    public void AirfoilFlow_NACA2412_HasPositiveLiftAtZeroAlpha()
    {
        // Cambered airfoil produces lift even at zero AoA
        var result = SimulationType.Create(MultiphysicsType.AirfoilFlow)
            .WithMaterial(EngineeringLibrary.Air)
            .WithAirfoil("2412")
            .WithAngleOfAttack(0.0)
            .WithFreestream(1.0)
            .Solve();

        // Cambered airfoils have positive Cl at α=0
        Assert.IsTrue(result.LiftCoefficient > 0,
            $"NACA 2412 should produce positive lift at α=0, got Cl={result.LiftCoefficient}");
    }

    [TestMethod]
    public void AirfoilFlow_DifferentPanelCounts_ConvergesLift()
    {
        double alpha = 4.0 * Math.PI / 180.0;

        var r40 = SimulationType.Create(MultiphysicsType.AirfoilFlow)
            .WithMaterial(EngineeringLibrary.Air)
            .WithAirfoil("0012", numPanels: 40)
            .WithAngleOfAttack(alpha)
            .WithFreestream(1.0)
            .Solve();

        var r100 = SimulationType.Create(MultiphysicsType.AirfoilFlow)
            .WithMaterial(EngineeringLibrary.Air)
            .WithAirfoil("0012", numPanels: 100)
            .WithAngleOfAttack(alpha)
            .WithFreestream(1.0)
            .Solve();

        // Both should give similar Cl (panel convergence)
        Assert.AreEqual(r40.LiftCoefficient, r100.LiftCoefficient, 0.15,
            "Panel method should converge with increasing panels");
    }
}
