using CSharpNumerics.Engines.Multiphysics;
using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials.Engineering;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class NewMultiphysicsTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Issue 5: VectorField.EvaluateGrid2D
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void VectorField_EvaluateGrid2D_ReturnsCorrectGridSize()
    {
        var field = new VectorField(
            p => p.x,
            p => p.y,
            p => 0.0);

        var result = field.EvaluateGrid2D(0, 1, 0, 1, 5, 4);

        Assert.AreEqual(20, result.Count); // 5 * 4 = 20
    }

    [TestMethod]
    public void VectorField_EvaluateGrid2D_EvaluatesCorrectly()
    {
        var field = new VectorField(
            p => 2.0 * p.x,
            p => -p.y,
            p => 0.0);

        var result = field.EvaluateGrid2D(0, 2, 0, 3, 3, 4);

        // Check corner: (2, 3) → (4, -3, 0)
        var corner = result.First(kv =>
            Math.Abs(kv.Key.x - 2.0) < 1e-10 && Math.Abs(kv.Key.y - 3.0) < 1e-10);
        Assert.AreEqual(4.0, corner.Value.x, 1e-10);
        Assert.AreEqual(-3.0, corner.Value.y, 1e-10);
    }

    [TestMethod]
    public void VectorField_EvaluateGrid2D_SinglePoint()
    {
        var field = new VectorField(
            p => 1.0,
            p => 2.0,
            p => 3.0);

        var result = field.EvaluateGrid2D(5, 5, 7, 7, 1, 1);

        Assert.AreEqual(1, result.Count);
        var entry = result.First();
        Assert.AreEqual(5.0, entry.Key.x, 1e-10);
        Assert.AreEqual(7.0, entry.Key.y, 1e-10);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Issue 4: Extended Material Library
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void EngineeringLibrary_HasNewMaterials()
    {
        Assert.AreEqual("Wood", EngineeringLibrary.Wood.Name);
        Assert.AreEqual("Rubber", EngineeringLibrary.Rubber.Name);
        Assert.AreEqual("Titanium", EngineeringLibrary.Titanium.Name);
        Assert.AreEqual("Brass", EngineeringLibrary.Brass.Name);
        Assert.AreEqual("Stainless Steel", EngineeringLibrary.StainlessSteel.Name);
        Assert.AreEqual("Oil", EngineeringLibrary.Oil.Name);
        Assert.AreEqual("Glycerin", EngineeringLibrary.Glycerin.Name);
        Assert.AreEqual("Plastic (HDPE)", EngineeringLibrary.Plastic.Name);
    }

    [TestMethod]
    public void EngineeringMaterial_HasMagneticPermeability()
    {
        Assert.AreEqual(100.0, EngineeringLibrary.Steel.MagneticPermeability);
        Assert.IsTrue(EngineeringLibrary.Copper.MagneticPermeability > 0);
        Assert.AreEqual(1.0, EngineeringLibrary.Brass.MagneticPermeability);
    }

    [TestMethod]
    public void EngineeringMaterial_CustomWithMagneticPermeability()
    {
        var mat = new EngineeringMaterial(
            "Iron", thermalConductivity: 80, specificHeat: 450, density: 7874,
            dynamicViscosity: 0, electricPermittivity: 0,
            youngsModulus: 211e9, poissonsRatio: 0.29,
            magneticPermeability: 5000.0);

        Assert.AreEqual(5000.0, mat.MagneticPermeability);
    }

    [TestMethod]
    public void EngineeringMaterial_DefaultMagneticPermeabilityIsOne()
    {
        var mat = new EngineeringMaterial(
            "Test", thermalConductivity: 1, specificHeat: 1, density: 1,
            dynamicViscosity: 0, electricPermittivity: 0,
            youngsModulus: 1, poissonsRatio: 0);

        Assert.AreEqual(1.0, mat.MagneticPermeability);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Issue 1: FluidFlow2D Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidFlow2D_BasicLidDrivenCavity()
    {
        var result = SimulationType.Create(MultiphysicsType.FluidFlow2D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
            .WithBoundary(top: 0.1, bottom: 0, left: 0, right: 0)
            .WithInletVelocity(0.01)
            .Solve(dt: 0.001, steps: 50);

        Assert.AreEqual(MultiphysicsType.FluidFlow2D, result.Type);
        Assert.IsNotNull(result.Vx);
        Assert.IsNotNull(result.Vy);
        Assert.IsNotNull(result.Pressure);
        Assert.IsNotNull(result.Timeline);
        Assert.AreEqual(20, result.Vx.GetLength(0));
        Assert.AreEqual(20, result.Vx.GetLength(1));
        Assert.IsTrue(result.Timeline.Count == 51); // initial + 50 steps
    }

    [TestMethod]
    public void FluidFlow2D_VelocityFieldsProduced()
    {
        var result = SimulationType.Create(MultiphysicsType.FluidFlow2D)
            .WithMaterial(EngineeringLibrary.Water)
            .WithGeometry(width: 0.05, height: 0.05, nx: 10, ny: 10)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 0)
            .WithInletVelocity(0.1)
            .Solve(dt: 0.0005, steps: 20);

        // Inlet should drive flow; Vx at interior should be non-zero
        bool hasFlow = false;
        for (int ix = 1; ix < 9; ix++)
            for (int iy = 1; iy < 9; iy++)
                if (Math.Abs(result.Vx[ix, iy]) > 1e-10)
                    hasFlow = true;

        Assert.IsTrue(hasFlow, "Expected non-zero interior velocity.");
    }

    [TestMethod]
    [ExpectedException(typeof(InvalidOperationException))]
    public void FluidFlow2D_ThrowsWithoutGeometry()
    {
        SimulationType.Create(MultiphysicsType.FluidFlow2D)
            .WithMaterial(EngineeringLibrary.Water)
            .Solve();
    }

    // ═══════════════════════════════════════════════════════════════
    //  Issue 2: MagneticField Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void MagneticField_SingleWire()
    {
        var result = SimulationType.Create(MultiphysicsType.MagneticField)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 30, ny: 30)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 0)
            .AddSource(15, 15, 1e6) // Current density at center
            .Solve(maxIterations: 5000, tolerance: 1e-6);

        Assert.AreEqual(MultiphysicsType.MagneticField, result.Type);
        Assert.IsNotNull(result.VectorPotential);
        Assert.IsNotNull(result.Bx);
        Assert.IsNotNull(result.By);
        Assert.AreEqual(30, result.Bx.GetLength(0));
        Assert.AreEqual(30, result.Bx.GetLength(1));

        // Vector potential should peak near the wire
        Assert.IsTrue(Math.Abs(result.VectorPotential[15, 15]) > 0);
    }

    [TestMethod]
    public void MagneticField_FieldDecaysFromSource()
    {
        var result = SimulationType.Create(MultiphysicsType.MagneticField)
            .WithMaterial(EngineeringLibrary.Copper)
            .WithGeometry(width: 0.2, height: 0.2, nx: 40, ny: 40)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 0)
            .AddSource(20, 20, 1e5)
            .Solve(maxIterations: 5000, tolerance: 1e-6);

        // |B| at center should be larger than at edge
        double bCenter = result.Field[20, 20];
        double bEdge = result.Field[5, 5];
        Assert.IsTrue(bCenter > bEdge, "Magnetic field should be stronger near the wire.");
    }

    [TestMethod]
    public void MagneticField_UsesMaterialPermeability()
    {
        // Steel (μ_r=100) should give stronger field than Aluminum (μ_r≈1)
        var resultSteel = SimulationType.Create(MultiphysicsType.MagneticField)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 0)
            .AddSource(10, 10, 1e5)
            .Solve(maxIterations: 3000, tolerance: 1e-5);

        var resultAlum = SimulationType.Create(MultiphysicsType.MagneticField)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.1, height: 0.1, nx: 20, ny: 20)
            .WithBoundary(top: 0, bottom: 0, left: 0, right: 0)
            .AddSource(10, 10, 1e5)
            .Solve(maxIterations: 3000, tolerance: 1e-5);

        Assert.IsTrue(resultSteel.MaxValue > resultAlum.MaxValue,
            "Steel (high μ_r) should produce stronger field than Aluminum.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Issue 3: PlaneStress Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void PlaneStress_CantileverPointLoad()
    {
        var result = SimulationType.Create(MultiphysicsType.PlaneStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 1.0, height: 0.1, nx: 40, ny: 10)
            .AddSource(30, 5, 1e6) // Point force at interior node
            .Solve(maxIterations: 5000, tolerance: 1e-8);

        Assert.AreEqual(MultiphysicsType.PlaneStress, result.Type);
        Assert.IsNotNull(result.Ux);
        Assert.IsNotNull(result.Uy);
        Assert.IsNotNull(result.StressXX);
        Assert.IsNotNull(result.StressYY);
        Assert.IsNotNull(result.StressXY);
        Assert.AreEqual(40, result.Ux.GetLength(0));
        Assert.AreEqual(10, result.Ux.GetLength(1));

        // Left boundary should be fixed (ux=0)
        Assert.AreEqual(0.0, result.Ux[0, 5]);
        Assert.AreEqual(0.0, result.Uy[0, 5]);
    }

    [TestMethod]
    public void PlaneStress_DisplacementIncreasesWithLoad()
    {
        var resultSmall = SimulationType.Create(MultiphysicsType.PlaneStress)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.5, height: 0.1, nx: 20, ny: 8)
            .AddSource(15, 4, 1e4)
            .Solve(maxIterations: 10000, tolerance: 1e-12);

        var resultLarge = SimulationType.Create(MultiphysicsType.PlaneStress)
            .WithMaterial(EngineeringLibrary.Aluminum)
            .WithGeometry(width: 0.5, height: 0.1, nx: 20, ny: 8)
            .AddSource(15, 4, 1e6)
            .Solve(maxIterations: 10000, tolerance: 1e-12);

        // Larger load → larger displacement
        double maxDispSmall = MaxAbs(resultSmall.Ux);
        double maxDispLarge = MaxAbs(resultLarge.Ux);
        Assert.IsTrue(maxDispLarge > maxDispSmall,
            $"Larger load should produce larger displacement. Small={maxDispSmall}, Large={maxDispLarge}");
    }

    [TestMethod]
    public void PlaneStress_VonMisesFieldPopulated()
    {
        var result = SimulationType.Create(MultiphysicsType.PlaneStress)
            .WithMaterial(EngineeringLibrary.Steel)
            .WithGeometry(width: 0.5, height: 0.1, nx: 20, ny: 8)
            .AddSource(15, 4, 1e5)
            .Solve(maxIterations: 5000, tolerance: 1e-9);

        // Field should contain von Mises stress > 0 somewhere
        Assert.IsTrue(result.MaxValue > 0, "Von Mises stress field should have positive values.");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Helpers
    // ═══════════════════════════════════════════════════════════════

    private static double MaxAbs(double[,] arr)
    {
        double max = 0;
        for (int i = 0; i < arr.GetLength(0); i++)
            for (int j = 0; j < arr.GetLength(1); j++)
            {
                double v = Math.Abs(arr[i, j]);
                if (v > max) max = v;
            }
        return max;
    }
}
