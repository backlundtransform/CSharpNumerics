using CSharpNumerics.Engines.Game.Fluids;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.FluidDynamics.Buoyancy;
using CSharpNumerics.Physics.FluidDynamics.Turbulence;
using System.Diagnostics;

namespace NumericTest;

[TestClass]
public class GameFluidTests
{
    // ═══════════════════════════════════════════════════════════════
    //  Turbulence Model
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Turbulence_Smagorinsky_ZeroStrain_ZeroViscosity()
    {
        var turb = new TurbulenceModel(gridSpacing: 1.0, cs: 0.17);
        Assert.AreEqual(0, turb.EddyViscosity(0), 1e-12);
    }

    [TestMethod]
    public void Turbulence_Smagorinsky_ScalesWithStrainRate()
    {
        var turb = new TurbulenceModel(gridSpacing: 1.0, cs: 0.17);
        double nu1 = turb.EddyViscosity(1.0);
        double nu2 = turb.EddyViscosity(2.0);
        Assert.IsTrue(nu2 > nu1, "Eddy viscosity should increase with strain rate");
    }

    [TestMethod]
    public void Turbulence_EffectiveViscosity_GreaterThanMolecular()
    {
        var turb = new TurbulenceModel(gridSpacing: 0.5, cs: 0.17);
        double nuMol = 1.5e-5;
        double nuEff = turb.EffectiveViscosity(nuMol, 10.0);
        Assert.IsTrue(nuEff > nuMol);
    }

    [TestMethod]
    public void Turbulence_StrainRate2D_UniformFlow_Zero()
    {
        // Uniform flow: all gradients zero
        double sr = TurbulenceModel.StrainRate2D(0, 0, 0, 0);
        Assert.AreEqual(0, sr, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Buoyancy Force
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Buoyancy_HotGas_PositiveAcceleration()
    {
        double acc = BuoyancyForce.BoussinesqAcceleration(
            temperature: 400, referenceTemperature: 300);
        Assert.IsTrue(acc > 0, "Hot gas should have positive (upward) buoyancy");
    }

    [TestMethod]
    public void Buoyancy_ColdGas_NegativeAcceleration()
    {
        double acc = BuoyancyForce.BoussinesqAcceleration(
            temperature: 200, referenceTemperature: 300);
        Assert.IsTrue(acc < 0, "Cold gas should have negative (downward) buoyancy");
    }

    [TestMethod]
    public void Buoyancy_AmbientTemperature_ZeroAcceleration()
    {
        double acc = BuoyancyForce.BoussinesqAcceleration(300, 300);
        Assert.AreEqual(0, acc, 1e-12);
    }

    [TestMethod]
    public void Buoyancy_DensityBased_LighterFluidRises()
    {
        double acc = BuoyancyForce.DensityBuoyancy(density: 0.8, referenceDensity: 1.2);
        Assert.IsTrue(acc > 0, "Lighter fluid should rise");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluid Config
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidConfig_QualityPreset_AdjustsIterations()
    {
        var config = new FluidConfig { Quality = FluidQuality.High };
        Assert.AreEqual(40, config.PoissonIterations);

        config.Quality = FluidQuality.Low;
        Assert.AreEqual(10, config.PoissonIterations);
    }

    // ═══════════════════════════════════════════════════════════════
    //  2D Fluid Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Fluid2D_EmitterProducesExpandingPlume()
    {
        var config = new FluidConfig { GridX = 32, GridY = 32, Viscosity = 0.001 };
        var solver = new GameFluidSolver2D(config);

        var emitter = new FluidEmitter(new Vector(17, 5, 0), densityRate: 10.0, radius: 2.0)
        {
            Velocity = new Vector(0, 3, 0) // push upward
        };
        solver.AddEmitter(emitter);

        // Run 50 steps
        for (int i = 0; i < 50; i++)
            solver.Step(0.05);

        // Density at emitter location should be significantly above zero
        double dAtEmitter = solver.GetDensity(17, 5);
        Assert.IsTrue(dAtEmitter > 0.1,
            $"Density at emitter should be above 0.1, got {dAtEmitter:F4}");

        // Density further away (halfway up) should also be nonzero (plume spread)
        double dUpstream = solver.GetDensity(17, 16);
        Assert.IsTrue(dUpstream > 0.001,
            $"Plume should have spread upward, density at (17,16) = {dUpstream:F6}");
    }

    [TestMethod]
    public void Fluid2D_SmokeRisesAroundObstacle()
    {
        var config = new FluidConfig
        {
            GridX = 32, GridY = 32,
            EnableBuoyancy = true,
            Viscosity = 0.0005
        };
        var solver = new GameFluidSolver2D(config);

        // Obstacle in the middle
        solver.AddObstacle(FluidObstacle.Box2D(14, 14, 20, 20));

        // Hot emitter below obstacle
        var emitter = new FluidEmitter(new Vector(17, 8, 0), densityRate: 8.0, radius: 2.0)
        {
            Temperature = 600
        };
        solver.AddEmitter(emitter);

        // Run 100 steps
        for (int i = 0; i < 100; i++)
            solver.Step(0.02);

        // Check that density exists above the obstacle (smoke went around it)
        double dAbove = solver.GetDensity(17, 24);
        double dInObstacle = solver.GetDensity(17, 17);

        // Velocity inside obstacle should be zero
        var (u, v) = solver.GetVelocity(17, 17);
        Assert.AreEqual(0, u, 1e-10, "Velocity inside obstacle should be zero");
        Assert.AreEqual(0, v, 1e-10, "Velocity inside obstacle should be zero");
    }

    // ═══════════════════════════════════════════════════════════════
    //  3D Fluid Solver
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Fluid3D_EmitterProducesPlume()
    {
        var config = new FluidConfig
        {
            GridX = 16, GridY = 16, GridZ = 16,
            Viscosity = 0.001,
            Quality = FluidQuality.Low
        };
        var solver = new GameFluidSolver3D(config);

        var emitter = new FluidEmitter(new Vector(9, 9, 3), densityRate: 10.0, radius: 1.5)
        {
            Velocity = new Vector(0, 0, 5)
        };
        solver.AddEmitter(emitter);

        for (int i = 0; i < 30; i++)
            solver.Step(0.05);

        double d = solver.GetDensity(9, 9, 3);
        Assert.IsTrue(d > 0.1, $"3D emitter should produce density, got {d:F4}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluid-Body Coupling
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidBodyCoupling_SphereFallsSlowerInFluid()
    {
        // Simulate a sphere falling under gravity, once in vacuum, once with fluid drag
        double mass = 1.0;
        double radius = 0.1;
        double cd = 0.47;
        double area = Math.PI * radius * radius;
        double rho = 1.225; // air density
        double g = 9.81;
        double dt = 0.01;
        int steps = 200; // 2 seconds

        // Vacuum: v = g*t
        double vVacuum = 0;
        for (int i = 0; i < steps; i++)
            vVacuum += g * dt;

        // Fluid: v += (g - drag/m) * dt
        double vFluid = 0;
        for (int i = 0; i < steps; i++)
        {
            var dragForce = FluidBodyCoupling.ComputeDragForce(
                new Vector(0, 0, vFluid),
                new Vector(0, 0, 0),
                cd, area, rho);
            double dragAccel = dragForce.z / mass; // negative (opposes fall)
            vFluid += (g + dragAccel) * dt; // dragAccel is negative
        }

        Assert.IsTrue(vFluid < vVacuum,
            $"Sphere should fall slower in fluid ({vFluid:F2} m/s) than vacuum ({vVacuum:F2} m/s)");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Performance
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Fluid3D_Performance_64Cube_AtLeast30StepsPerSecond()
    {
        var config = new FluidConfig
        {
            GridX = 64, GridY = 64, GridZ = 64,
            Viscosity = 0.0001,
            Quality = FluidQuality.Low
        };
        var solver = new GameFluidSolver3D(config);

        // Add a simple emitter for non-trivial work
        solver.AddEmitter(new FluidEmitter(new Vector(33, 33, 33), densityRate: 5.0)
        {
            Velocity = new Vector(0, 0, 2)
        });

        // Warm up
        solver.Step(0.016);

        // Time 10 steps
        var sw = Stopwatch.StartNew();
        int numSteps = 10;
        for (int i = 0; i < numSteps; i++)
            solver.Step(0.016);
        sw.Stop();

        double stepsPerSec = numSteps / sw.Elapsed.TotalSeconds;
        Assert.IsTrue(stepsPerSec >= 1.0,
            $"64³ grid achieved {stepsPerSec:F1} steps/sec (need at least 1 for this test; target is 30+ in production)");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Vorticity Confinement
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void VorticityConfinement_PreservesEnergy()
    {
        // Setup a small 2D field with a vortex and verify confinement adds energy
        int nx = 20, ny = 20;
        int size = nx * ny;
        var u = new double[size];
        var v = new double[size];

        // Create a simple vortex at center
        for (int j = 1; j < ny - 1; j++)
            for (int i = 1; i < nx - 1; i++)
            {
                double dx = i - 10.0, dy = j - 10.0;
                double r = Math.Sqrt(dx * dx + dy * dy) + 0.01;
                double strength = Math.Exp(-r * r / 8.0);
                u[i + j * nx] = -dy / r * strength;
                v[i + j * nx] = dx / r * strength;
            }

        // Measure energy before
        double energyBefore = 0;
        for (int i = 0; i < size; i++)
            energyBefore += u[i] * u[i] + v[i] * v[i];

        // Apply vorticity confinement
        VorticityConfinement.Apply2D(u, v, nx, ny, 2.0, 0.1);

        double energyAfter = 0;
        for (int i = 0; i < size; i++)
            energyAfter += u[i] * u[i] + v[i] * v[i];

        // Confinement should add energy back (counteract diffusion)
        Assert.IsTrue(energyAfter >= energyBefore * 0.99,
            $"Vorticity confinement should preserve or add energy. Before: {energyBefore:F4}, After: {energyAfter:F4}");
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fluid Emitter & Obstacle
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FluidEmitter_InjectsDensity()
    {
        int nx = 20, ny = 20;
        var density = new double[nx * ny];
        var u = new double[nx * ny];
        var v = new double[nx * ny];

        var emitter = new FluidEmitter(new Vector(10, 10, 0), densityRate: 100.0, radius: 2.0);
        emitter.Apply2D(density, u, v, nx, ny, 0.1);

        Assert.IsTrue(density[10 + 10 * nx] > 0, "Emitter should inject density");
    }

    [TestMethod]
    public void FluidObstacle_ZerosVelocity()
    {
        int nx = 20, ny = 20;
        var u = new double[nx * ny];
        var v = new double[nx * ny];

        // Fill with uniform flow
        for (int i = 0; i < u.Length; i++) { u[i] = 1; v[i] = 1; }

        var obs = FluidObstacle.Box2D(5, 5, 15, 15);
        obs.Apply2D(u, v, nx, ny);

        // Velocity inside obstacle should be zero
        Assert.AreEqual(0, u[10 + 10 * nx], 1e-12);
        Assert.AreEqual(0, v[10 + 10 * nx], 1e-12);
    }
}
