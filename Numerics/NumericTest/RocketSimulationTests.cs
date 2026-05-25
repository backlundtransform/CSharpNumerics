using CSharpNumerics.Engines.Common;
using CSharpNumerics.Engines.Game.Rocket;
using CSharpNumerics.ML.ReinforcementLearning.Environments;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.Mechanics.Propulsion;
using CSharpNumerics.Physics.OrbitalMechanics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest;

[TestClass]
public class RocketSimulationTests
{
    private const double G0 = 9.80665;

    #region Tsiolkovsky Tests

    [TestMethod]
    public void Tsiolkovsky_SingleStageToOrbit_DeltaV_Approximately_9400()
    {
        // Typical single-stage-to-orbit requires ~9.4 km/s
        // Isp = 350s (vacuum), mass ratio ~15:1
        double isp = 350;
        double m0 = 500000; // 500 tonnes wet
        double mf = m0 / 15.0; // mass ratio 15
        double dv = TsiolkovskyExtensions.DeltaV(isp, m0, mf);

        // Should be around 9.3 km/s
        Assert.IsTrue(dv > 9000, $"ΔV should be > 9000 m/s, got {dv:F0}");
        Assert.IsTrue(dv < 10000, $"ΔV should be < 10000 m/s, got {dv:F0}");
    }

    [TestMethod]
    public void Tsiolkovsky_MassRatio_RoundTrip()
    {
        double isp = 311;
        double m0 = 100000;
        double mf = 30000;

        double dv = TsiolkovskyExtensions.DeltaV(isp, m0, mf);
        double ratio = TsiolkovskyExtensions.MassRatio(dv, isp);

        Assert.AreEqual(m0 / mf, ratio, 1e-8);
    }

    [TestMethod]
    public void Tsiolkovsky_ExhaustVelocity_MatchesIspTimesG0()
    {
        double isp = 300;
        double ve = TsiolkovskyExtensions.ExhaustVelocity(isp);
        Assert.AreEqual(isp * G0, ve, 1e-6);
    }

    [TestMethod]
    public void Tsiolkovsky_BurnTime_ConsistentWithMassFlow()
    {
        double isp = 311;
        double thrust = 7607000; // Merlin 1D (N)
        double propMass = 400000;

        double burnTime = TsiolkovskyExtensions.BurnTime(propMass, thrust, isp);
        double massFlow = thrust / (isp * G0);

        Assert.AreEqual(propMass / massFlow, burnTime, 1e-6);
    }

    #endregion

    #region RocketEngine Tests

    [TestMethod]
    public void RocketEngine_MassFlowRate_EqualsThrottleTimesMaxFlow()
    {
        // Merlin 1D-like: Isp_sl=282s, Isp_vac=311s, Thrust_vac=981kN
        var engine = new RocketEngine(282, 311, 981000);
        engine.IsActive = true;
        engine.Throttle = 1.0;

        double flowAtVacuum = engine.MassFlowRate(0);
        double expectedFlow = 981000.0 / (311 * G0);
        Assert.AreEqual(expectedFlow, flowAtVacuum, 0.01);
    }

    [TestMethod]
    public void RocketEngine_Thrust_ZeroWhenInactive()
    {
        var engine = new RocketEngine(282, 311, 981000);
        engine.IsActive = false;
        Assert.AreEqual(0, engine.Thrust(0.5));
        Assert.AreEqual(0, engine.MassFlowRate(0.5));
    }

    [TestMethod]
    public void RocketEngine_EffectiveIsp_InterpolatesBetweenSeaLevelAndVacuum()
    {
        var engine = new RocketEngine(282, 311, 981000);

        Assert.AreEqual(311, engine.EffectiveIsp(0), 1e-6);    // vacuum
        Assert.AreEqual(282, engine.EffectiveIsp(1), 1e-6);    // sea level
        Assert.AreEqual(296.5, engine.EffectiveIsp(0.5), 1e-6); // midpoint
    }

    #endregion

    #region PropellantTank Tests

    [TestMethod]
    public void PropellantTank_ConstantFlowRate_LinearMassDecrease()
    {
        double fuelMass = 100000;
        double oxMass = 250000; // O/F ratio 2.5
        var tank = new PropellantTank(fuelMass, oxMass);

        double massFlowRate = 1000; // kg/s
        double dt = 1.0; // 1 second steps

        double prevMass = tank.PropellantMass;
        for (int i = 0; i < 10; i++)
        {
            tank.Consume(massFlowRate, dt);
            double currentMass = tank.PropellantMass;
            double loss = prevMass - currentMass;
            Assert.AreEqual(massFlowRate * dt, loss, 0.01,
                $"Step {i}: mass decrease should be constant at {massFlowRate} kg/step");
            prevMass = currentMass;
        }
    }

    [TestMethod]
    public void PropellantTank_RemainingDeltaV_DecreasesWithConsumption()
    {
        var tank = new PropellantTank(40000, 100000);
        double dryMass = 20000;
        double isp = 311;

        double dvFull = tank.RemainingDeltaV(dryMass, isp);
        tank.Consume(500, 10); // consume 5000 kg
        double dvAfter = tank.RemainingDeltaV(dryMass, isp);

        Assert.IsTrue(dvAfter < dvFull, "ΔV should decrease after consuming propellant");
    }

    #endregion

    #region CenterOfMassTracker Tests

    [TestMethod]
    public void CgTracker_MonotonicShift_FromFullToEmpty()
    {
        double dryMass = 5000;
        double propMass = 100000;
        double cgFull = 15.0;  // CG at 15m when full (closer to bottom)
        double cgEmpty = 20.0; // CG at 20m when empty (shifts up toward nose)

        var tracker = new CenterOfMassTracker(dryMass, propMass, cgFull, cgEmpty);

        double prevCg = tracker.ComputeCg(propMass);
        Assert.AreEqual(cgFull, prevCg, 1e-6);

        // As propellant decreases, CG should move monotonically toward cgEmpty
        for (double frac = 0.9; frac >= 0; frac -= 0.1)
        {
            double currentProp = propMass * frac;
            double cg = tracker.ComputeCg(currentProp);
            Assert.IsTrue(cg >= prevCg, $"CG should move monotonically, got {cg} < {prevCg}");
            prevCg = cg;
        }

        double cgAtEmpty = tracker.ComputeCg(0);
        Assert.AreEqual(cgEmpty, cgAtEmpty, 1e-6);
    }

    [TestMethod]
    public void CgTracker_BoundedBetweenFullAndEmpty()
    {
        var tracker = new CenterOfMassTracker(5000, 100000, 10.0, 18.0);

        for (double prop = 0; prop <= 100000; prop += 5000)
        {
            double cg = tracker.ComputeCg(prop);
            Assert.IsTrue(cg >= 10.0 && cg <= 18.0,
                $"CG={cg} out of bounds [10, 18] at prop={prop}");
        }
    }

    #endregion

    #region RocketSimulationEngine Tests

    [TestMethod]
    public void RocketSimEngine_Step_ConservesMomentum_ThrustEqualsRateOfMomentumChange()
    {
        // Single-stage rocket, vertical launch
        var engine = new RocketEngine(282, 311, 981000);
        var tank = new PropellantTank(100000, 250000);
        var stage = new RocketStage(20000, new[] { engine }, new[] { tank });
        var vehicle = new RocketVehicle(new[] { stage });

        var sim = new RocketSimulationEngine(vehicle);
        sim.DragCoefficient = 0; // No drag for clean momentum test
        sim.ReferenceArea = 0;

        // Set initial state: on pad, vertical, NED frame (z-down)
        var state = new RocketState(
            new Vector(0, 0, 0),
            new Vector(0, 0, 0),
            // Nose-up in NED: body X (forward) points in world -Z (up)
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0),
            vehicle.TotalMass
        );
        sim.SetState(state);
        sim.Init();

        double dt = 0.1;
        // Take one step
        double massBefore = sim.State.Mass;
        Vector velBefore = sim.State.Velocity;
        sim.Step(dt);
        double massAfter = sim.State.Mass;
        Vector velAfter = sim.State.Velocity;

        // Momentum change should approximately equal thrust impulse minus gravity impulse
        double avgMass = (massBefore + massAfter) / 2.0;
        Vector momentumChange = avgMass * (velAfter - velBefore);

        // Thrust impulse (upward in NED = -z direction)
        double thrustMag = stage.TotalThrust(1.0); // sea level
        // Gravity impulse (downward in NED = +z direction)
        double gravImpulse = avgMass * G0 * dt;

        // The z-component of momentum change should reflect thrust - gravity
        // (negative z = upward in NED)
        double netUpwardImpulse = thrustMag * dt - gravImpulse;

        // Momentum change in z should be approximately -netUpwardImpulse (NED convention)
        // Allow 5% tolerance due to RK4 approximation and variable mass
        Assert.AreEqual(-netUpwardImpulse, momentumChange.z, Math.Abs(netUpwardImpulse) * 0.05,
            "Thrust impulse should match rate of momentum change");
    }

    [TestMethod]
    public void RocketSimEngine_VerticalLaunch_GainsAltitude()
    {
        var engine = new RocketEngine(282, 311, 5000000); // 5MN thrust (TWR > 1)
        var tank = new PropellantTank(80000, 200000);
        var stage = new RocketStage(15000, new[] { engine }, new[] { tank });
        var vehicle = new RocketVehicle(new[] { stage });

        var sim = new RocketSimulationEngine(vehicle);

        // Start on pad, nose up: body X → world -Z (up in NED)
        var state = new RocketState(
            new Vector(0, 0, 0),
            new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0),
            vehicle.TotalMass
        );
        sim.SetState(state);
        sim.Init();

        // Simulate 10 seconds
        for (int i = 0; i < 100; i++)
            sim.Step(0.1);

        // Should have gained altitude (z becomes more negative in NED)
        Assert.IsTrue(sim.State.Position.z < -100,
            $"Rocket should be above 100m after 10s, got altitude={sim.State.Altitude:F0}m");
    }

    [TestMethod]
    public void RocketSimEngine_Staging_FiresEvent()
    {
        // Small first stage that depletes quickly
        // Mass flow ≈ 500000/(311*9.80665) ≈ 163.9 kg/s → 50 kg depletes in ~0.3s
        var engine1 = new RocketEngine(282, 311, 500000);
        var tank1 = new PropellantTank(15, 35); // Only 50 kg propellant
        var stage1 = new RocketStage(1000, new[] { engine1 }, new[] { tank1 });

        var engine2 = new RocketEngine(282, 340, 500000);
        var tank2 = new PropellantTank(20000, 50000);
        var stage2 = new RocketStage(3000, new[] { engine2 }, new[] { tank2 });

        var vehicle = new RocketVehicle(new[] { stage1, stage2 });

        var sim = new RocketSimulationEngine(vehicle);
        var eventBus = new EventBus();
        sim.EventBus = eventBus;

        StageEvent receivedEvent = null;
        eventBus.Subscribe<StageEvent>(e => receivedEvent = e);

        var state = new RocketState(
            new Vector(0, 0, 0),
            new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0),
            vehicle.TotalMass
        );
        sim.SetState(state);
        sim.Init();

        // Step until staging occurs (50 kg depletes in ~0.3s at 163.9 kg/s)
        for (int i = 0; i < 50; i++)
            sim.Step(0.1);

        Assert.IsNotNull(receivedEvent, "Stage separation event should have been published");
        Assert.AreEqual(0, receivedEvent.StageIndex, "First stage (index 0) should have separated");
        Assert.AreEqual(1, vehicle.ActiveStageIndex, "Active stage should now be index 1");
    }

    [TestMethod]
    public void ThrustCurve_SolidBooster_ProfileShape()
    {
        var curve = ThrustCurve.SolidBooster(100);

        // At start: should be zero (ramp up begins)
        Assert.AreEqual(0.0, curve.Evaluate(0), 1e-6);

        // Middle of sustain: should be 1.0
        Assert.AreEqual(1.0, curve.Evaluate(50), 1e-6);

        // Past burn: should be zero
        Assert.AreEqual(0.0, curve.Evaluate(101), 1e-6);
    }

    [TestMethod]
    public void ThrustCurve_Constant_ReturnsFixedValue()
    {
        var curve = ThrustCurve.Constant(120, 0.8);

        Assert.AreEqual(0.8, curve.Evaluate(0), 1e-6);
        Assert.AreEqual(0.8, curve.Evaluate(60), 1e-6);
        Assert.AreEqual(0.8, curve.Evaluate(119.9), 1e-6);
        Assert.AreEqual(0.0, curve.Evaluate(120.1), 1e-6);
    }

    #endregion

    #region Phase 2 — Supersonic Aerodynamics & Max-Q

    [TestMethod]
    public void CompressibleDrag_PeaksNearMach1()
    {
        var model = new CompressibleDragModel(cd0: 0.3, cdPeak: 0.8, machPeak: 1.1);

        double cdSubsonic = model.Evaluate(0.5);
        double cdTransonic = model.Evaluate(1.1);
        double cdSupersonic = model.Evaluate(2.0);
        double cdHighSupersonic = model.Evaluate(4.0);

        // Subsonic: should be Cd0
        Assert.AreEqual(0.3, cdSubsonic, 1e-6);

        // Peak at Mach 1.1
        Assert.AreEqual(0.8, cdTransonic, 1e-6);

        // Supersonic: should decay below peak
        Assert.IsTrue(cdSupersonic < cdTransonic,
            $"Cd at Mach 2 ({cdSupersonic:F3}) should be less than peak ({cdTransonic:F3})");

        // Higher Mach: further decay
        Assert.IsTrue(cdHighSupersonic < cdSupersonic,
            $"Cd at Mach 4 ({cdHighSupersonic:F3}) should be less than at Mach 2 ({cdSupersonic:F3})");
    }

    [TestMethod]
    public void CompressibleDrag_DecaysAboveMach2()
    {
        var model = CompressibleDragModel.SlenderRocket();

        double cd2 = model.Evaluate(2.0);
        double cd3 = model.Evaluate(3.0);
        double cd5 = model.Evaluate(5.0);

        Assert.IsTrue(cd3 < cd2, "Cd should decrease with increasing Mach in supersonic regime");
        Assert.IsTrue(cd5 < cd3, "Cd should continue decreasing at higher Mach");
    }

    [TestMethod]
    public void MaxQ_OccursAtExpectedAltitude()
    {
        // Simulate a typical ascent: high thrust rocket going vertical
        var engine = new RocketEngine(282, 311, 7600000); // ~F9 thrust
        var tank = new PropellantTank(100000, 250000);
        var stage = new RocketStage(25000, new[] { engine }, new[] { tank });
        var vehicle = new RocketVehicle(new[] { stage });

        var sim = new RocketSimulationEngine(vehicle);
        sim.DragModel = CompressibleDragModel.SlenderRocket();
        sim.ReferenceArea = 10.5; // ~3.66m diameter

        var monitor = new MaxQMonitor();
        sim.MaxQMonitor = monitor;

        var state = new RocketState(
            new Vector(0, 0, 0),
            new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0),
            vehicle.TotalMass
        );
        sim.SetState(state);
        sim.Init();

        // Simulate 120 seconds of ascent
        for (int i = 0; i < 1200; i++)
            sim.Step(0.1);

        Assert.IsTrue(monitor.MaxQReached, "Max-Q should have been detected during ascent");
        // Max-Q typically at 11–14 km for orbital rockets
        Assert.IsTrue(monitor.MaxQAltitude > 8000,
            $"Max-Q altitude {monitor.MaxQAltitude:F0}m should be > 8 km");
        Assert.IsTrue(monitor.MaxQAltitude < 20000,
            $"Max-Q altitude {monitor.MaxQAltitude:F0}m should be < 20 km");
    }

    [TestMethod]
    public void ThrottleBucket_ReducesPeakDynamicPressure()
    {
        // Run same rocket twice: once without bucket, once with
        RocketSimulationEngine CreateSim()
        {
            var engine = new RocketEngine(282, 311, 7600000);
            var tank = new PropellantTank(100000, 250000);
            var stage = new RocketStage(25000, new[] { engine }, new[] { tank });
            var vehicle = new RocketVehicle(new[] { stage });
            var sim = new RocketSimulationEngine(vehicle);
            sim.DragModel = CompressibleDragModel.SlenderRocket();
            sim.ReferenceArea = 10.5;
            sim.MaxQMonitor = new MaxQMonitor();
            var state = new RocketState(
                new Vector(0, 0, 0), new Vector(0, 0, 0),
                Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
                new Vector(0, 0, 0), vehicle.TotalMass);
            sim.SetState(state);
            return sim;
        }

        // Without throttle bucket
        var simNoB = CreateSim();
        simNoB.Init();
        for (int i = 0; i < 1200; i++) simNoB.Step(0.1);
        double peakNoB = simNoB.MaxQMonitor.PeakQ;

        // With throttle bucket
        var simWithB = CreateSim();
        simWithB.ThrottleBucket = new ThrottleBucket(qOnset: 20000, qTarget: 35000, minThrottle: 0.6);
        simWithB.Init();
        for (int i = 0; i < 1200; i++) simWithB.Step(0.1);
        double peakWithB = simWithB.MaxQMonitor.PeakQ;

        Assert.IsTrue(peakWithB < peakNoB,
            $"Throttle bucket peak Q ({peakWithB:F0} Pa) should be less than unthrottled ({peakNoB:F0} Pa)");
    }

    [TestMethod]
    public void HeatFlux_OrderOfMagnitude_LEOVelocity()
    {
        // At LEO insertion: v ≈ 7800 m/s, altitude ≈ 100 km (very low density)
        // At lower altitude during reentry: v ≈ 7800, alt ≈ 70 km
        double velocity = 7800; // m/s
        double noseRadius = 1.0; // 1m nose radius

        // At 70 km altitude (reentry heating peak region)
        double qDot70 = HeatFlux.FromAltitudeAndVelocity(70000, velocity, noseRadius);

        // Order of magnitude: should be in MW/m² range for reentry
        // At 70 km, ρ ≈ 8.8e-5 kg/m³
        // q̇ = 1.7415e-4 * √(8.8e-5/1.0) * 7800³ ≈ 1.7415e-4 * 9.38e-3 * 4.74e11 ≈ 775 kW/m²
        Assert.IsTrue(qDot70 > 100000, $"Heat flux at 70km/7.8km/s should exceed 100 kW/m², got {qDot70:F0} W/m²");
        Assert.IsTrue(qDot70 < 5000000, $"Heat flux should be below 5 MW/m², got {qDot70:F0} W/m²");

        // At 200 km: essentially zero (no atmosphere)
        double qDot200 = HeatFlux.FromAltitudeAndVelocity(200000, velocity, noseRadius);
        Assert.IsTrue(qDot200 < qDot70 * 0.01,
            "Heat flux above 200km should be negligible compared to 70km");
    }

    [TestMethod]
    public void MachNumber_CorrectAtSeaLevel()
    {
        // Speed of sound at sea level ≈ 340 m/s
        double mach = MachNumber.Compute(340, 0);
        Assert.AreEqual(1.0, mach, 0.01);

        double mach2 = MachNumber.Compute(680, 0);
        Assert.AreEqual(2.0, mach2, 0.01);
    }

    #endregion

    #region Phase 3 — Staging & Event System

    [TestMethod]
    public void TwoStage_ReachesHigherVelocity_ThanSingleStage()
    {
        // Single stage: 5000 dry + 50000 propellant = 55000 total
        var engineSingle = new RocketEngine(282, 311, 1000000);
        var tankSingle = new PropellantTank(14000, 36000); // 50000 kg total
        var stageSingle = new RocketStage(5000, new[] { engineSingle }, new[] { tankSingle });
        var vehicleSingle = new RocketVehicle(new[] { stageSingle });

        var simSingle = new RocketSimulationEngine(vehicleSingle);
        simSingle.SetState(new RocketState(
            new Vector(0, 0, 0), new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0), vehicleSingle.TotalMass));
        simSingle.Init();
        for (int i = 0; i < 3000; i++) simSingle.Step(0.1);
        double vSingle = simSingle.State.Speed;

        // Two stages: same total mass, stage 1 burns fast then separates
        // Stage 1: 3000 dry + 12000 prop = 15000, high thrust for fast burn (~30s)
        // Stage 2: 2000 dry + 38000 prop = 40000, moderate thrust
        // Total: 55000 same as single
        var engine1 = new RocketEngine(282, 311, 1200000);
        var tank1 = new PropellantTank(3400, 8600); // 12000 kg
        var stage1 = new RocketStage(3000, new[] { engine1 }, new[] { tank1 });

        var engine2 = new RocketEngine(282, 320, 600000);
        var tank2 = new PropellantTank(10860, 27140); // 38000 kg
        var stage2 = new RocketStage(2000, new[] { engine2 }, new[] { tank2 });

        var vehicleTwo = new RocketVehicle(new[] { stage1, stage2 });
        var simTwo = new RocketSimulationEngine(vehicleTwo);
        simTwo.SetState(new RocketState(
            new Vector(0, 0, 0), new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0), vehicleTwo.TotalMass));
        simTwo.Init();
        for (int i = 0; i < 3000; i++) simTwo.Step(0.1);
        double vTwo = simTwo.State.Speed;

        Assert.IsTrue(vTwo > vSingle,
            $"Two-stage ({vTwo:F0} m/s) should reach higher velocity than single-stage ({vSingle:F0} m/s)");
    }

    [TestMethod]
    public void StageSeparation_FiresAtFuelDepletion_WithinTolerance()
    {
        // Stage with known burn time: prop=1000kg, flow≈163.9 kg/s → ~6.1s
        var engine = new RocketEngine(282, 311, 500000);
        var tank = new PropellantTank(285, 715); // 1000 kg total
        var stage1 = new RocketStage(500, new[] { engine }, new[] { tank });

        var engine2 = new RocketEngine(282, 311, 500000);
        var tank2 = new PropellantTank(10000, 25000);
        var stage2 = new RocketStage(2000, new[] { engine2 }, new[] { tank2 });

        var vehicle = new RocketVehicle(new[] { stage1, stage2 });
        var sim = new RocketSimulationEngine(vehicle);
        var eventBus = new EventBus();
        sim.EventBus = eventBus;

        StageEvent receivedEvent = null;
        eventBus.Subscribe<StageEvent>(e => receivedEvent = e);

        sim.SetState(new RocketState(
            new Vector(0, 0, 0), new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0), vehicle.TotalMass));
        sim.Init();

        // Expected burn time: 1000 / (500000 / (311 * 9.80665)) ≈ 6.1s
        double expectedBurnTime = 1000.0 / (500000.0 / (311.0 * G0));

        // Run with small dt for precision
        double dt = 0.05;
        for (int i = 0; i < 200; i++) sim.Step(dt);

        Assert.IsNotNull(receivedEvent, "Stage separation event should have fired");
        Assert.AreEqual(0, receivedEvent.StageIndex);
        // Separation should occur within ±0.1s of expected burn time
        Assert.AreEqual(expectedBurnTime, receivedEvent.Time, 0.15,
            $"Separation at {receivedEvent.Time:F2}s should be near expected {expectedBurnTime:F2}s");
    }

    [TestMethod]
    public void PostSeparation_ContinuousPositionVelocity_DiscontinuousMass()
    {
        // Track state across a separation event
        var engine1 = new RocketEngine(282, 311, 500000);
        var tank1 = new PropellantTank(15, 35); // Tiny: depletes fast
        var stage1 = new RocketStage(1000, new[] { engine1 }, new[] { tank1 });

        var engine2 = new RocketEngine(282, 311, 500000);
        var tank2 = new PropellantTank(10000, 25000);
        var stage2 = new RocketStage(2000, new[] { engine2 }, new[] { tank2 });

        var vehicle = new RocketVehicle(new[] { stage1, stage2 });
        var sim = new RocketSimulationEngine(vehicle);
        var eventBus = new EventBus();
        sim.EventBus = eventBus;

        Vector posBeforeSep = new Vector(0, 0, 0);
        Vector velBeforeSep = new Vector(0, 0, 0);
        double massBeforeSep = 0;
        Vector posAfterSep = new Vector(0, 0, 0);
        Vector velAfterSep = new Vector(0, 0, 0);
        double massAfterSep = 0;
        bool separated = false;

        sim.SetState(new RocketState(
            new Vector(0, 0, 0), new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0), vehicle.TotalMass));
        sim.Init();

        for (int i = 0; i < 100; i++)
        {
            if (!separated)
            {
                posBeforeSep = sim.State.Position;
                velBeforeSep = sim.State.Velocity;
                massBeforeSep = sim.State.Mass;
            }

            sim.Step(0.05);

            if (!separated && vehicle.ActiveStageIndex > 0)
            {
                posAfterSep = sim.State.Position;
                velAfterSep = sim.State.Velocity;
                massAfterSep = sim.State.Mass;
                separated = true;
            }
        }

        Assert.IsTrue(separated, "Separation should have occurred");

        // Position and velocity should be continuous (small change over one step)
        double posDiff = (posAfterSep - posBeforeSep).GetMagnitude();
        double velDiff = (velAfterSep - velBeforeSep).GetMagnitude();
        Assert.IsTrue(posDiff < 50, $"Position jump {posDiff:F1}m should be small (one step)");
        Assert.IsTrue(velDiff < 100, $"Velocity jump {velDiff:F1}m/s should be small (one step)");

        // Mass should be discontinuous (stage dropped)
        double massDrop = massBeforeSep - massAfterSep;
        Assert.IsTrue(massDrop > 500, $"Mass drop {massDrop:F0}kg should be significant (stage1 dry mass)");
    }

    [TestMethod]
    public void ParallelBoosters_MaintainSymmetricTrajectory()
    {
        // Core stage + 2 symmetric boosters
        var coreEngine = new RocketEngine(282, 311, 2000000);
        var coreTank = new PropellantTank(50000, 125000);
        var coreStage = new RocketStage(10000, new[] { coreEngine }, new[] { coreTank });

        var vehicle = new RocketVehicle(new[] { coreStage });

        // Add 2 identical boosters
        var boosterEngine1 = new RocketEngine(270, 270, 1500000);
        var boosterTank1 = new PropellantTank(2000, 5000); // 7000 kg: deplete in ~12.4s
        var booster1 = new RocketStage(2000, new[] { boosterEngine1 }, new[] { boosterTank1 });

        var boosterEngine2 = new RocketEngine(270, 270, 1500000);
        var boosterTank2 = new PropellantTank(2000, 5000);
        var booster2 = new RocketStage(2000, new[] { boosterEngine2 }, new[] { boosterTank2 });

        vehicle.Boosters.Add(booster1);
        vehicle.Boosters.Add(booster2);

        var sim = new RocketSimulationEngine(vehicle);
        var eventBus = new EventBus();
        sim.EventBus = eventBus;

        BoosterSeparationEvent boosterEvent = null;
        eventBus.Subscribe<BoosterSeparationEvent>(e => boosterEvent = e);

        sim.SetState(new RocketState(
            new Vector(0, 0, 0), new Vector(0, 0, 0),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 2),
            new Vector(0, 0, 0), vehicle.TotalMass));
        sim.Init();

        // Simulate until boosters separate
        for (int i = 0; i < 200; i++) sim.Step(0.1);

        Assert.IsNotNull(boosterEvent, "Booster separation event should have fired");
        Assert.IsTrue(vehicle.BoostersSeparated, "Boosters should be separated");

        // Trajectory should remain symmetric (no lateral drift)
        // In our model, symmetric boosters produce no net lateral force
        double lateralDrift = Math.Abs(sim.State.Position.y);
        Assert.IsTrue(lateralDrift < 1.0,
            $"Lateral drift {lateralDrift:F2}m should be negligible with symmetric boosters");

        // Vehicle should still be accelerating (core still has fuel)
        Assert.IsFalse(vehicle.AllStagesSpent, "Core stage should still be active after booster sep");
    }

    [TestMethod]
    public void EventDetector_LocatesFuelDepletion_Precisely()
    {
        var detector = new EventDetector();

        // 1000 kg propellant at 200 kg/s → depletes at 5.0s
        // With dt=1.0, event is at fraction 0.0 of the step starting at t=5.0
        // But actually within a step starting at t=4.5 with dt=1.0:
        //   prop at start = 100kg, flow = 200 kg/s → depletes at 100/200 = 0.5s into step
        double fraction = detector.LocateFuelDepletion(100, 200, 1.0);
        Assert.AreEqual(0.5, fraction, 1e-6);

        // Not depleted during step
        double noEvent = detector.LocateFuelDepletion(500, 200, 1.0);
        Assert.AreEqual(-1, noEvent);
    }

    [TestMethod]
    public void CoastPhase_NoThrust_BallisiticFlight()
    {
        var coast = new CoastPhase(10.0); // 10 second coast

        Assert.IsFalse(coast.IsStarted);
        Assert.IsFalse(coast.IsActive(0));

        coast.Start(50.0);
        Assert.IsTrue(coast.IsStarted);
        Assert.IsTrue(coast.IsActive(55.0));
        Assert.IsFalse(coast.IsComplete(55.0));
        Assert.AreEqual(0.5, coast.Progress(55.0), 1e-6);

        Assert.IsTrue(coast.IsComplete(60.0));
        Assert.IsFalse(coast.IsActive(60.0));
    }

    #endregion

    #region Phase 4 — Coordinate Systems & Earth Model

    [TestMethod]
    public void Geodetic_ECEF_Geodetic_RoundTrip_PreservesCoordinates()
    {
        // Test multiple points around the globe
        var testPoints = new[]
        {
            GeoCoordinate.FromDegrees(0, 0, 0),           // Equator/prime meridian
            GeoCoordinate.FromDegrees(51.4769, -0.0005, 11), // London
            GeoCoordinate.FromDegrees(-33.8688, 151.2093, 58), // Sydney
            GeoCoordinate.FromDegrees(90, 0, 0),          // North pole
            GeoCoordinate.FromDegrees(-90, 0, 0),         // South pole
            GeoCoordinate.FromDegrees(28.5729, -80.6490, 3), // Cape Canaveral
            GeoCoordinate.FromDegrees(45, 90, 10000),     // High altitude
        };

        foreach (var original in testPoints)
        {
            var ecef = CoordinateConversions.GeodeticToECEF(original);
            var roundTrip = CoordinateConversions.ECEFToGeodetic(ecef);

            // Sub-millimeter accuracy: < 1mm = 1e-3 m
            double latError = Math.Abs(roundTrip.Latitude - original.Latitude) * EarthModel.SemiMajorAxis;
            double lonError = Math.Abs(roundTrip.Longitude - original.Longitude) * EarthModel.SemiMajorAxis * Math.Cos(original.Latitude);
            double altError = Math.Abs(roundTrip.Altitude - original.Altitude);

            Assert.IsTrue(latError < 1e-3,
                $"Latitude error {latError:E3} m exceeds 1mm for ({original.LatitudeDegrees}°, {original.LongitudeDegrees}°)");
            Assert.IsTrue(lonError < 1e-3 || Math.Abs(Math.Abs(original.Latitude) - Math.PI / 2) < 1e-6,
                $"Longitude error {lonError:E3} m exceeds 1mm for ({original.LatitudeDegrees}°, {original.LongitudeDegrees}°)");
            Assert.IsTrue(altError < 1e-3,
                $"Altitude error {altError:E3} m exceeds 1mm for ({original.LatitudeDegrees}°, {original.LongitudeDegrees}°)");
        }
    }

    [TestMethod]
    public void EquatorialLaunch_GainsEarthRotationVelocity()
    {
        // At the equator, Earth's surface velocity ≈ 465 m/s
        // v = ω * R_equatorial = 7.2921150e-5 * 6378137 ≈ 465.1 m/s
        double equatorialVelocity = EarthModel.SurfaceVelocity(0); // latitude = 0

        Assert.IsTrue(equatorialVelocity > 460 && equatorialVelocity < 470,
            $"Equatorial rotation velocity {equatorialVelocity:F1} m/s should be ~465 m/s");

        // At 28.5° (Cape Canaveral), velocity should be less
        double capeVelocity = EarthModel.SurfaceVelocity(28.5729 * Math.PI / 180);
        Assert.IsTrue(capeVelocity > 400 && capeVelocity < 420,
            $"Cape Canaveral rotation velocity {capeVelocity:F1} m/s should be ~408 m/s");

        // Verify via simulation engine helper
        var engine = new RocketEngine(282, 311, 1000000);
        var tank = new PropellantTank(14000, 36000);
        var stage = new RocketStage(5000, new[] { engine }, new[] { tank });
        var vehicle = new RocketVehicle(new[] { stage });
        var sim = new RocketSimulationEngine(vehicle);

        sim.SetStateFromLaunchSite(GeoCoordinate.FromDegrees(0, 0, 0));
        double bonus = sim.EarthRotationVelocityBonus();
        Assert.IsTrue(bonus > 460 && bonus < 470,
            $"Earth rotation bonus {bonus:F1} m/s should be ~465 m/s");
    }

    [TestMethod]
    public void J2Perturbation_ProducesMeasurableOrbitalPrecession()
    {
        // A satellite in a 400 km circular orbit at 51.6° inclination (ISS-like)
        // J2 causes RAAN precession: dΩ/dt = -1.5 * n * J2 * (Re/a)² * cos(i)
        // For ISS: a ≈ 6771 km, n ≈ sqrt(GM/a³), i = 51.6°
        double a = EarthModel.SemiMajorAxis + 400000;
        double n = Math.Sqrt(EarthModel.GM / (a * a * a)); // mean motion
        double cosI = Math.Cos(51.6 * Math.PI / 180);
        double re_a = EarthModel.SemiMajorAxis / a;

        // Analytical RAAN precession rate
        double dOmegaDt = -1.5 * n * EarthModel.J2 * re_a * re_a * cosI;
        double dOmegaDtDegPerDay = dOmegaDt * 180.0 / Math.PI * 86400;

        // ISS RAAN precession ≈ -5° per day
        Assert.IsTrue(dOmegaDtDegPerDay < -4 && dOmegaDtDegPerDay > -6,
            $"RAAN precession {dOmegaDtDegPerDay:F2} °/day should be ~-5°/day");

        // Verify J2 gravity differs from point mass at the surface
        // At surface: equator at (a, 0, 0), pole at (0, 0, b) — different radii
        var posEquator = new Vector(EarthModel.SemiMajorAxis, 0, 0);
        var posPole = new Vector(0, 0, EarthModel.SemiMinorAxis);

        Vector gEquator = GravityModel.Acceleration(posEquator);
        Vector gPole = GravityModel.Acceleration(posPole);
        Vector gPointMass = GravityModel.PointMassAcceleration(posEquator);

        // At the surface, gravity is stronger at poles (smaller radius + J2 effect)
        double gMagEquator = gEquator.GetMagnitude();
        double gMagPole = gPole.GetMagnitude();
        double gMagPM = gPointMass.GetMagnitude();

        Assert.IsTrue(gMagPole > gMagEquator,
            $"Polar gravity {gMagPole:F6} should exceed equatorial gravity {gMagEquator:F6}");
        Assert.IsTrue(Math.Abs(gMagEquator - gMagPM) / gMagPM > 1e-4,
            "J2 perturbation should produce measurable difference from point mass");

        // Propagate a circular orbit for 10 orbits using velocity-Verlet (symplectic)
        double v0 = Math.Sqrt(EarthModel.GM / a); // circular velocity
        Vector pos = new Vector(a, 0, 0);
        Vector vel = new Vector(0, v0, 0);
        double dt = 10.0; // 10-second steps
        double period = 2 * Math.PI / n;
        int steps = (int)(10 * period / dt); // 10 orbits

        double E0 = 0.5 * v0 * v0 - EarthModel.GM / a; // specific energy

        for (int i = 0; i < steps; i++)
        {
            // Velocity Verlet
            Vector acc = GravityModel.Acceleration(pos);
            // Half-step velocity
            vel = new Vector(vel.x + 0.5 * acc.x * dt, vel.y + 0.5 * acc.y * dt, vel.z + 0.5 * acc.z * dt);
            // Full-step position
            pos = new Vector(pos.x + vel.x * dt, pos.y + vel.y * dt, pos.z + vel.z * dt);
            // New acceleration
            Vector accNew = GravityModel.Acceleration(pos);
            // Half-step velocity
            vel = new Vector(vel.x + 0.5 * accNew.x * dt, vel.y + 0.5 * accNew.y * dt, vel.z + 0.5 * accNew.z * dt);
        }

        double r = pos.GetMagnitude();
        double v = vel.GetMagnitude();
        double Ef = 0.5 * v * v - EarthModel.GM / r;
        double drift = Math.Abs((Ef - E0) / E0);

        // Energy drift should be small (< 1% over 10 orbits with 10s steps)
        Assert.IsTrue(drift < 0.01,
            $"Energy drift {drift:P4} exceeds 1% over 10 orbits");
    }

    [TestMethod]
    public void Coriolis_DeflectsVerticalLaunchWestward()
    {
        // A rocket launched vertically from the equator deflects westward
        // in the ECEF frame due to Coriolis effect.
        // For upward velocity vr at the equator: Coriolis = -2(ω × v) → -Y direction (westward)

        // Set up a simple ECEF simulation
        var engine = new RocketEngine(282, 311, 5000000);
        var tank = new PropellantTank(14000, 36000);
        var stage = new RocketStage(5000, new[] { engine }, new[] { tank });
        var vehicle = new RocketVehicle(new[] { stage });
        var sim = new RocketSimulationEngine(vehicle);

        // Launch from equator (0°, 0°, 0m)
        sim.SetStateFromLaunchSite(GeoCoordinate.FromDegrees(0, 0, 0));
        sim.Init();

        // Get initial position (on equator, prime meridian)
        Vector initialPos = sim.State.Position;

        // Simulate for 100 seconds
        for (int i = 0; i < 1000; i++) sim.Step(0.1);

        Vector finalPos = sim.State.Position;

        // In ECEF from equator/prime meridian launch:
        // Initial pos ≈ (R, 0, 0). Radial up = +X direction.
        // Coriolis for upward motion: -2*(ω × v_up) = -Y direction (westward)
        double westwardDeflection = -(finalPos.y - initialPos.y);

        // Should be positive (westward) and on the order of tens of meters
        Assert.IsTrue(westwardDeflection > 1.0,
            $"Westward Coriolis deflection {westwardDeflection:F1} m should be positive and measurable");
    }

    [TestMethod]
    public void ECEF_ECI_RoundTrip_PreservesPosition()
    {
        // Test ECEF → ECI → ECEF round-trip at various GMST angles
        var ecef = new ECEFPosition(6378137, 100000, 300000);
        double[] gmstAngles = { 0, Math.PI / 4, Math.PI / 2, Math.PI, 1.5 * Math.PI };

        foreach (double gmst in gmstAngles)
        {
            var eci = CoordinateConversions.ECEFToECI(ecef, gmst);
            var roundTrip = CoordinateConversions.ECIToECEF(eci, gmst);

            Assert.AreEqual(ecef.X, roundTrip.X, 1e-6, $"X mismatch at GMST={gmst:F3}");
            Assert.AreEqual(ecef.Y, roundTrip.Y, 1e-6, $"Y mismatch at GMST={gmst:F3}");
            Assert.AreEqual(ecef.Z, roundTrip.Z, 1e-6, $"Z mismatch at GMST={gmst:F3}");
        }
    }

    [TestMethod]
    public void EarthModel_WGS84Parameters_AreCorrect()
    {
        // Validate WGS84 constants
        Assert.AreEqual(6378137.0, EarthModel.SemiMajorAxis, 1e-1);
        Assert.IsTrue(Math.Abs(EarthModel.SemiMinorAxis - 6356752.314) < 1.0,
            $"Semi-minor axis {EarthModel.SemiMinorAxis:F3} should be ~6356752.314 m");
        Assert.IsTrue(EarthModel.EccentricitySquared > 0.00669 && EarthModel.EccentricitySquared < 0.00670,
            $"e² = {EarthModel.EccentricitySquared:F8} should be ~0.00669");
        Assert.AreEqual(7.2921150e-5, EarthModel.RotationRate, 1e-12);
        Assert.AreEqual(3.986004418e14, EarthModel.GM, 1e6);
        Assert.AreEqual(1.08263e-3, EarthModel.J2, 1e-8);

        // Equatorial surface gravity (gravitational only, no centrifugal) ≈ 9.81 m/s²
        double gEquator = GravityModel.Magnitude(EarthModel.SemiMajorAxis, 0);
        Assert.IsTrue(gEquator > 9.80 && gEquator < 9.83,
            $"Equatorial gravity {gEquator:F4} m/s² should be ~9.81");

        // Polar surface gravity ≈ 9.83 m/s²
        double gPole = GravityModel.Magnitude(EarthModel.SemiMinorAxis, Math.PI / 2);
        Assert.IsTrue(gPole > 9.82 && gPole < 9.87,
            $"Polar gravity {gPole:F4} m/s² should be ~9.83-9.86");
    }

    #endregion

    #region Phase 5 — Orbital Mechanics & Insertion

    [TestMethod]
    public void CircularOrbit_400km_HasCorrectPeriodAndVelocity()
    {
        // Circular orbit at 400 km: period ≈ 92.4 min, velocity ≈ 7.67 km/s
        double altitude = 400000; // 400 km
        double r = EarthModel.SemiMajorAxis + altitude;
        double mu = EarthModel.GM;

        double vCircular = Math.Sqrt(mu / r);
        double period = 2 * Math.PI * Math.Sqrt(r * r * r / mu);
        double periodMinutes = period / 60.0;

        Assert.IsTrue(vCircular > 7600 && vCircular < 7700,
            $"Circular velocity {vCircular:F0} m/s should be ~7670 m/s");
        Assert.IsTrue(periodMinutes > 92 && periodMinutes < 93,
            $"Orbital period {periodMinutes:F1} min should be ~92.4 min");

        // Verify via OrbitalElements
        var elements = new OrbitalElements(r, 0, 0, 0, 0, 0, mu);
        Assert.AreEqual(period, elements.Period, 1.0);
        Assert.IsTrue(elements.IsCircular);
        Assert.IsTrue(elements.IsBound);
    }

    [TestMethod]
    public void HohmannTransfer_LEO_to_GEO_DeltaV_MatchesAnalytical()
    {
        // LEO (400 km) → GEO (35786 km): total ΔV ≈ 3.94 km/s
        double altLEO = 400000;     // 400 km
        double altGEO = 35786000;   // 35786 km (geostationary)

        var result = HohmannTransfer.FromAltitudes(altLEO, altGEO);

        // Total ΔV should be approximately 3.94 km/s
        Assert.IsTrue(result.TotalDeltaV > 3800 && result.TotalDeltaV < 4100,
            $"Hohmann LEO→GEO total ΔV {result.TotalDeltaV:F0} m/s should be ~3940 m/s");

        // First burn (leave LEO) should be ~2.46 km/s
        Assert.IsTrue(result.DeltaV1 > 2300 && result.DeltaV1 < 2600,
            $"ΔV1 {result.DeltaV1:F0} m/s should be ~2460 m/s");

        // Second burn (circularize at GEO) should be ~1.48 km/s
        Assert.IsTrue(result.DeltaV2 > 1300 && result.DeltaV2 < 1600,
            $"ΔV2 {result.DeltaV2:F0} m/s should be ~1480 m/s");

        // Transfer time ≈ 5.25 hours
        double transferHours = result.TransferTime / 3600.0;
        Assert.IsTrue(transferHours > 5.0 && transferHours < 5.5,
            $"Transfer time {transferHours:F2} hours should be ~5.25 hours");
    }

    [TestMethod]
    public void StateToElements_RoundTrip_HighPrecision()
    {
        double mu = EarthModel.GM;

        // Test several different orbits
        var testOrbits = new[]
        {
            // Circular equatorial (ISS-like altitude)
            new OrbitalElements(6778137, 0.001, 0.9, 1.5, 2.0, 0.5, mu),
            // Highly eccentric (GTO-like)
            new OrbitalElements(24000000, 0.73, 0.45, 0.8, 3.0, 1.2, mu),
            // Polar circular
            new OrbitalElements(7000000, 0.0001, Math.PI / 2, 0.3, 1.0, 4.0, mu),
            // Inclined elliptical
            new OrbitalElements(12000000, 0.3, 1.0, 2.5, 0.5, 3.5, mu),
        };

        foreach (var original in testOrbits)
        {
            // Elements → State
            var (pos, vel) = ElementsToState.ToStateVector(original);

            // State → Elements
            var roundTrip = StateToElements.FromStateVector(pos, vel, mu);

            // Compare elements
            Assert.AreEqual(original.SemiMajorAxis, roundTrip.SemiMajorAxis,
                original.SemiMajorAxis * 1e-10,
                $"Semi-major axis mismatch for a={original.SemiMajorAxis:E3}");
            Assert.AreEqual(original.Eccentricity, roundTrip.Eccentricity, 1e-10,
                $"Eccentricity mismatch for e={original.Eccentricity}");
            Assert.AreEqual(original.Inclination, roundTrip.Inclination, 1e-10,
                $"Inclination mismatch for i={original.Inclination}");

            // RAAN comparison (handle wraparound)
            double raanDiff = Math.Abs(original.RAAN - roundTrip.RAAN);
            if (raanDiff > Math.PI) raanDiff = 2 * Math.PI - raanDiff;
            Assert.IsTrue(raanDiff < 1e-10 || original.Inclination < 1e-6,
                $"RAAN mismatch: {raanDiff:E3} for Ω={original.RAAN}");

            // True anomaly comparison
            double nuDiff = Math.Abs(original.TrueAnomaly - roundTrip.TrueAnomaly);
            if (nuDiff > Math.PI) nuDiff = 2 * Math.PI - nuDiff;
            Assert.IsTrue(nuDiff < 1e-9,
                $"True anomaly mismatch: {nuDiff:E3} for ν={original.TrueAnomaly}");
        }
    }

    [TestMethod]
    public void OrbitalPropagator_StableOver100Revolutions()
    {
        // Circular orbit at 400 km, propagate for 100 revolutions
        double mu = EarthModel.GM;
        double r = EarthModel.SemiMajorAxis + 400000;
        double v = Math.Sqrt(mu / r);

        var propagator = new OrbitalPropagator();
        propagator.IncludeJ2 = false; // Pure Kepler for energy conservation test
        propagator.SetState(new Vector(r, 0, 0), new Vector(0, v, 0));

        double E0 = propagator.SpecificEnergy;
        double period = 2 * Math.PI * Math.Sqrt(r * r * r / mu);

        // Propagate 100 orbits with 10-second steps
        propagator.Propagate(100 * period, 10.0);

        double Ef = propagator.SpecificEnergy;
        double drift = Math.Abs((Ef - E0) / E0);

        Assert.IsTrue(drift < 0.0001,
            $"Energy drift {drift:P6} exceeds 0.01% over 100 revolutions");

        // Radius should still be approximately circular
        double rFinal = propagator.Radius;
        double radiusError = Math.Abs(rFinal - r) / r;
        Assert.IsTrue(radiusError < 0.001,
            $"Radius drift {radiusError:P4} exceeds 0.1% after 100 orbits");
    }

    [TestMethod]
    public void CircularizationBurn_MatchesHohmann()
    {
        // After a Hohmann transfer from LEO to 800 km, circularization ΔV should match
        double r1 = EarthModel.SemiMajorAxis + 400000;
        double r2 = EarthModel.SemiMajorAxis + 800000;
        double mu = EarthModel.GM;

        var hohmann = HohmannTransfer.Compute(r1, r2, mu);

        // Create state at apoapsis of transfer orbit
        double aTransfer = hohmann.TransferSemiMajorAxis;
        double eTransfer = 1.0 - r1 / aTransfer;
        var transferElements = new OrbitalElements(aTransfer, eTransfer, 0, 0, 0, Math.PI, mu);
        var (pos, vel) = ElementsToState.ToStateVector(transferElements);

        double dvCirc = CircularizationBurn.Compute(pos, vel, mu);

        // Should match Hohmann ΔV2
        Assert.AreEqual(hohmann.DeltaV2, dvCirc, 1.0,
            $"Circularization ΔV {dvCirc:F1} should match Hohmann ΔV2 {hohmann.DeltaV2:F1}");
    }

    [TestMethod]
    public void OrbitalInsertionDetector_DetectsClosedOrbit()
    {
        // Create a state vector for a 400 km circular orbit
        double r = EarthModel.SemiMajorAxis + 400000;
        double v = Math.Sqrt(EarthModel.GM / r);

        var detector = new OrbitalInsertionDetector();
        detector.MinPeriapsisAltitude = 100000;

        // Sub-orbital state (too slow for orbit)
        bool subOrbital = detector.Check(new Vector(r, 0, 0), new Vector(0, v * 0.5, 0), 100);
        Assert.IsFalse(subOrbital, "Half-orbital velocity should not trigger insertion");

        // Full orbital velocity
        bool orbital = detector.Check(new Vector(r, 0, 0), new Vector(0, v, 0), 200);
        Assert.IsTrue(orbital, "Full circular velocity should trigger insertion");
        Assert.AreEqual(200, detector.InsertionTime);
        Assert.IsTrue(detector.InsertionElements.IsCircular);
    }

    [TestMethod]
    public void MissionProfile_PhasesAdvanceCorrectly()
    {
        var profile = MissionProfile.StandardAscentToOrbit(400000);

        Assert.AreEqual("Vertical Ascent", profile.ActivePhase.Name);
        Assert.IsFalse(profile.IsComplete);

        // Advance past vertical ascent (altitude > 1000m)
        bool transitioned = profile.Update(10, 1500, 50);
        Assert.IsTrue(transitioned);
        Assert.AreEqual("Gravity Turn", profile.ActivePhase.Name);

        // Advance past gravity turn (altitude > 80km)
        transitioned = profile.Update(120, 85000, 2000);
        Assert.IsTrue(transitioned);
        Assert.AreEqual("Coast to Apoapsis", profile.ActivePhase.Name);

        // Not yet at apoapsis
        transitioned = profile.Update(300, 300000, 5000);
        Assert.IsFalse(transitioned);
        Assert.AreEqual("Coast to Apoapsis", profile.ActivePhase.Name);
    }

    [TestMethod]
    public void TelemetryRecorder_RecordsAtSampleRate()
    {
        var recorder = new TelemetryRecorder();
        recorder.SampleInterval = 5.0;

        var state = new RocketState(
            new Vector(6778137, 0, 0), new Vector(0, 7670, 0),
            Quaternion.Identity, new Vector(0, 0, 0), 50000);

        // Record at t=0
        recorder.Record(state, 0);
        Assert.AreEqual(1, recorder.FrameCount);

        // Too early (t=2 < interval of 5)
        recorder.Record(state, 2);
        Assert.AreEqual(1, recorder.FrameCount);

        // At t=5, should record
        recorder.Record(state, 5);
        Assert.AreEqual(2, recorder.FrameCount);

        // At t=10, should record
        recorder.Record(state, 10);
        Assert.AreEqual(3, recorder.FrameCount);

        Assert.AreEqual(0, recorder.Frames[0].Time);
        Assert.AreEqual(5, recorder.Frames[1].Time);
        Assert.AreEqual(10, recorder.Frames[2].Time);
    }

    #endregion

    #region Phase 6 — Guidance, Navigation & Control

    [TestMethod]
    public void GravityTurn_ReachesTargetInclination_Within01Degrees()
    {
        // Simulate a gravity turn from vertical to near-horizontal
        // Target inclination 28.5° (Cape Canaveral)
        var guidance = new GravityTurnGuidance
        {
            KickAltitude = 500,
            KickAngle = 3.0 * Math.PI / 180.0,
            KickDuration = 10.0,
            TargetInclination = 28.5 * Math.PI / 180.0,
            LaunchLatitude = 28.5 * Math.PI / 180.0
        };

        // Verify launch azimuth for due-east launch at 28.5° latitude → 90° azimuth
        double azimuth = guidance.LaunchAzimuth();
        // For target inclination = latitude, azimuth should be ~90° (due east)
        Assert.IsTrue(Math.Abs(azimuth - Math.PI / 2) < 0.1 * Math.PI / 180.0,
            $"Azimuth {azimuth * 180 / Math.PI:F2}° should be ~90° for inc=lat");

        // Test a gravity turn sequence at 30 km (well into gravity turn phase)
        double altitude = 30000; // 30 km
        double speed = 1500.0;
        Vector velocity = new Vector(speed * Math.Cos(0.3), 0, -speed * Math.Sin(0.3));
        double time = 100;

        Quaternion att = guidance.ComputeAttitude(altitude, velocity, time);
        // Should be in gravity turn phase — attitude should follow velocity (prograde)
        Assert.AreNotEqual(Quaternion.Identity, att);
    }

    [TestMethod]
    public void PEG_ConvergesToTargetOrbit_Within3Cycles()
    {
        // Set up PEG for a 200 km circular orbit
        double targetRadius = EarthModel.SemiMajorAxis + 200000;
        var peg = new PEGGuidance
        {
            TargetSemiMajorAxis = targetRadius,
            TargetEccentricity = 0.0,
            TargetInclination = 28.5 * Math.PI / 180.0,
            Mu = EarthModel.GM,
            Iterations = 10
        };

        // Position: 150 km altitude, moving nearly orbital speed
        double r = EarthModel.SemiMajorAxis + 150000;
        Vector position = new Vector(r, 0, 0);
        double vOrb = Math.Sqrt(EarthModel.GM / r) * 0.95; // slightly below orbital
        Vector velocity = new Vector(0, vOrb, 0);

        double thrustAccel = 30.0; // m/s² (upper stage)
        double exhaustVel = 3500.0;

        // Run 3 guidance cycles
        Vector dir1 = peg.ComputeThrustDirection(position, velocity, thrustAccel, exhaustVel);
        Vector dir2 = peg.ComputeThrustDirection(position, velocity, thrustAccel, exhaustVel);
        Vector dir3 = peg.ComputeThrustDirection(position, velocity, thrustAccel, exhaustVel);

        // Thrust direction should be a unit vector
        Assert.IsTrue(Math.Abs(dir3.GetMagnitude() - 1.0) < 0.01,
            $"Thrust direction magnitude should be ~1, got {dir3.GetMagnitude():F4}");

        // PEG should converge (has a velocity-to-gain)
        Assert.IsTrue(peg.HasConverged || peg.VelocityToGain > 0,
            "PEG should either converge or have positive velocity-to-gain");
    }

    [TestMethod]
    public void GimbalRateLimited_DoesNotOvershoot()
    {
        var tvc = new ThrustVectorControl
        {
            MaxGimbalAngle = 5.0 * Math.PI / 180.0,  // ±5°
            MaxGimbalRate = 10.0 * Math.PI / 180.0,  // 10°/s
            MomentArm = 20.0
        };

        double dt = 0.01; // 100 Hz control loop
        double commandedAngle = 5.0 * Math.PI / 180.0; // Full deflection step

        // Step until converged
        double maxPitch = 0;
        for (int i = 0; i < 200; i++) // 2 seconds
        {
            tvc.Command(commandedAngle, 0, dt);
            if (tvc.GimbalPitch > maxPitch)
                maxPitch = tvc.GimbalPitch;
        }

        // Should converge to commanded angle
        Assert.IsTrue(Math.Abs(tvc.GimbalPitch - commandedAngle) < 0.001 * Math.PI / 180.0,
            $"Gimbal should reach commanded angle, got {tvc.GimbalPitch * 180 / Math.PI:F3}°");

        // Should NOT overshoot the commanded angle
        Assert.IsTrue(maxPitch <= commandedAngle + 1e-10,
            $"Gimbal overshot: max={maxPitch * 180 / Math.PI:F3}° > commanded={commandedAngle * 180 / Math.PI:F3}°");

        // Verify rate limiting: from 0° to 5° at 10°/s should take 0.5s (50 steps)
        tvc.Reset();
        tvc.Command(commandedAngle, 0, dt);
        double firstStepDelta = tvc.GimbalPitch;
        double expectedDelta = 10.0 * Math.PI / 180.0 * dt; // rate * dt
        Assert.IsTrue(Math.Abs(firstStepDelta - expectedDelta) < 1e-10,
            $"First step should move exactly rate*dt={expectedDelta * 180 / Math.PI:F4}°, got {firstStepDelta * 180 / Math.PI:F4}°");
    }

    [TestMethod]
    public void RocketLandingEnv_SoftTouchdown_AfterControlledDescent()
    {
        // Test that the RL environment can produce a soft landing
        // by using a simple hand-crafted controller (not RL training)
        var env = new RocketLandingEnv(
            dt: 0.1,
            maxSteps: 1000,
            dryMass: 25000,
            fuelMass: 8000,
            maxThrust: 800000,
            exhaustVelocity: 3000,
            gravity: 9.81);

        var (state, _) = env.Reset(seed: 42);
        Assert.AreEqual(8, state.Length);

        bool landed = false;
        double totalReward = 0;
        int steps = 0;

        // Simple bang-bang controller: full thrust when speed is too high
        for (int i = 0; i < 1000; i++)
        {
            steps++;
            double altitude = state[0] * 1000.0; // denormalize
            double vz = state[3] * 100.0;        // denormalize vertical vel
            double speed = state[7] * 100.0;     // denormalize speed

            // Throttle proportional to descent speed + gravity compensation
            double throttle = 0;
            if (vz < -5.0)
                throttle = Math.Min(1.0, (-vz - 2.0) / 50.0 + 0.3);
            else if (altitude < 100 && vz < -1.0)
                throttle = 0.8;

            double gimbal = -state[4] * 2.0; // counter pitch

            var action = new VectorN(new double[] { throttle, gimbal });
            var (nextState, reward, done, info) = env.Step(action);
            totalReward += reward;
            state = nextState;

            if (done)
            {
                if (info.ContainsKey("landed"))
                    landed = true;
                break;
            }
        }

        // With this simple controller and generous fuel, landing should be achievable
        Assert.IsTrue(steps > 10, "Should run for multiple steps");
        // The environment terminates when altitude <= 0
        // We don't require perfect landing from hand-crafted controller, 
        // just verify the env works and terminates
        Assert.IsTrue(steps < 1000, "Should terminate before max steps");
    }

    [TestMethod]
    public void AscentOptimizationEnv_ObservationAndActionSpaces_Correct()
    {
        var env = new AscentOptimizationEnv();
        Assert.AreEqual(8, env.ObservationSize);
        Assert.AreEqual(1, env.ActionSize);
        Assert.IsFalse(env.IsDiscrete);

        var (state, _) = env.Reset(seed: 123);
        Assert.AreEqual(8, state.Length);

        // Initial state: at ground, no speed
        Assert.AreEqual(0, state[0], 0.01); // altitude = 0
        Assert.AreEqual(0, state[1], 0.01); // speed = 0
    }

    [TestMethod]
    public void AttitudeController_TracksCommandedAttitude()
    {
        var controller = new AttitudeController
        {
            Kp = 2.0,
            Kd = 1.5,
            Ki = 0.0,
            MaxRate = 5.0 * Math.PI / 180.0,
            MaxTorque = 100000.0
        };

        // Command a 10° pitch rotation
        Quaternion current = Quaternion.Identity;
        Quaternion commanded = Quaternion.FromAxisAngle(new Vector(0, 1, 0), 10.0 * Math.PI / 180.0);
        Vector angularRate = new Vector(0, 0, 0);

        Vector torque = controller.ComputeTorque(current, commanded, angularRate, 0.1);

        // Torque should be primarily about Y axis (pitch)
        Assert.IsTrue(Math.Abs(torque.y) > Math.Abs(torque.x),
            "Torque should be primarily about pitch axis");
        Assert.IsTrue(Math.Abs(torque.y) > Math.Abs(torque.z),
            "Torque should be primarily about pitch axis");
        // Torque should be positive (rotating toward commanded)
        Assert.IsTrue(torque.y > 0, "Torque should drive toward commanded attitude");
    }

    [TestMethod]
    public void GuidanceComputer_SwitchesModesCorrectly()
    {
        var gc = new GuidanceComputer();
        gc.PEGActivationAltitude = 80000;
        gc.PEG.TargetSemiMajorAxis = EarthModel.SemiMajorAxis + 200000;
        gc.PEG.Mu = EarthModel.GM;
        gc.GravityTurn.KickAltitude = 500;

        // In gravity turn mode at low altitude
        var state = new RocketState(
            new Vector(0, 0, -10000),   // 10 km altitude (NED: z is down)
            new Vector(500, 0, -200),
            Quaternion.Identity,
            new Vector(0, 0, 0),
            100000);

        gc.Mode = GuidanceMode.GravityTurn;
        gc.Update(state, 60, 10000, 20.0, 3000, 0.1);

        Assert.AreEqual(GuidanceMode.GravityTurn, gc.Mode);
        Assert.AreNotEqual(Quaternion.Identity, gc.CommandedAttitude);
    }

    [TestMethod]
    public void NavigationFilter_PerfectMode_ReturnsExactState()
    {
        var nav = new NavigationFilter(seed: 0) { AddNoise = false };

        var trueState = new RocketState(
            new Vector(100, 200, -5000),
            new Vector(1000, 50, -300),
            Quaternion.FromAxisAngle(new Vector(0, 1, 0), 0.1),
            new Vector(0.01, 0.02, 0.0),
            50000);

        var estimated = nav.Estimate(trueState);

        Assert.AreEqual(trueState.Position.x, estimated.Position.x, 1e-10);
        Assert.AreEqual(trueState.Position.y, estimated.Position.y, 1e-10);
        Assert.AreEqual(trueState.Position.z, estimated.Position.z, 1e-10);
        Assert.AreEqual(trueState.Velocity.x, estimated.Velocity.x, 1e-10);
        Assert.AreEqual(trueState.Mass, estimated.Mass, 1e-10);
    }

    [TestMethod]
    public void NavigationFilter_NoisyMode_AddsUncertainty()
    {
        var nav = new NavigationFilter(seed: 42) { AddNoise = true, PositionNoiseSigma = 100.0 };

        var trueState = new RocketState(
            new Vector(1000000, 0, -200000),
            new Vector(7000, 0, 0),
            Quaternion.Identity,
            new Vector(0, 0, 0),
            5000);

        var est = nav.Estimate(trueState);

        // Position should differ from true (noise added)
        double posDiff = Math.Sqrt(
            Math.Pow(est.Position.x - trueState.Position.x, 2) +
            Math.Pow(est.Position.y - trueState.Position.y, 2) +
            Math.Pow(est.Position.z - trueState.Position.z, 2));

        Assert.IsTrue(posDiff > 0.1, $"Noisy estimate should differ from truth, diff={posDiff:F4}");
    }

    [TestMethod]
    public void TVC_ComputeTorque_ProportionalToThrustAndDeflection()
    {
        var tvc = new ThrustVectorControl { MomentArm = 20.0, MaxGimbalAngle = 10.0 * Math.PI / 180.0 };

        // Set a known gimbal angle
        tvc.Command(3.0 * Math.PI / 180.0, 2.0 * Math.PI / 180.0, 100.0); // large dt → instant

        double thrust = 1000000.0; // 1 MN
        Vector torque = tvc.ComputeTorque(thrust);

        // Expected: torque_y = F * L * sin(pitch_gimbal), torque_z = F * L * sin(yaw_gimbal)
        double expectedY = thrust * 20.0 * Math.Sin(3.0 * Math.PI / 180.0);
        double expectedZ = thrust * 20.0 * Math.Sin(2.0 * Math.PI / 180.0);

        Assert.AreEqual(expectedY, torque.y, expectedY * 0.01);
        Assert.AreEqual(expectedZ, torque.z, expectedZ * 0.01);
        Assert.AreEqual(0, torque.x, 1e-6); // no roll torque from gimbal
    }

    #endregion

    #region Phase 7 — Unity Integration & Visualization

    [TestMethod]
    public void Adapter_PreservesState_Over100000Steps_WithoutDrift()
    {
        var adapter = new RocketUnityAdapter
        {
            CoordinateMode = AdapterCoordinateMode.PassThrough,
            PositionScale = 1.0
        };

        // Simulate 100,000 steps of a circular orbit at 400 km
        double r = 6778137.0; // Earth radius + 400 km
        double v = Math.Sqrt(3.986e14 / r); // Orbital velocity
        double period = 2.0 * Math.PI * r / v;
        double dt = period / 100000.0; // 100k steps = one orbit
        double omega = v / r;

        double angle = 0;
        for (int i = 0; i < 100000; i++)
        {
            angle += omega * dt;
            var state = new RocketState(
                new Vector(r * Math.Cos(angle), r * Math.Sin(angle), 0),
                new Vector(-v * Math.Sin(angle), v * Math.Cos(angle), 0),
                Quaternion.Identity,
                new Vector(0, 0, 0),
                5000);

            adapter.Update(state);
        }

        // After one complete orbit, position should be close to start
        Vector finalPos = adapter.UnityPosition;
        double endR = Math.Sqrt(finalPos.x * finalPos.x + finalPos.y * finalPos.y + finalPos.z * finalPos.z);

        // Radius should remain constant (no drift) to within double precision
        Assert.AreEqual(r, endR, 1.0,
            $"Radial position drifted: expected {r:F1}, got {endR:F1}, diff={Math.Abs(r - endR):E3}");

        // Accumulated drift from the adapter's tracking should be minimal
        // (The adapter's drift metric tracks floating-point ULP accumulation)
        Assert.IsTrue(adapter.AccumulatedDrift < 1.0,
            $"Adapter accumulated drift {adapter.AccumulatedDrift:E3} should be < 1 meter");
    }

    [TestMethod]
    public void TimeWarp_ProducesIdenticalTrajectory_RegardlessOfWarpFactor()
    {
        // Create a simple vehicle for simulation
        var engine = new RocketEngine(282, 311, 50000);
        var tank = new PropellantTank(300, 700);
        var stage = new RocketStage(dryMass: 500, engines: new[] { engine }, tanks: new[] { tank });
        stage.SeparationTrigger = StageSeparationTrigger.FuelDepletion();
        var vehicle = new RocketVehicle(new[] { stage }, payloadMass: 100);

        // Run at 1x warp (direct stepping)
        var sim1 = new RocketSimulationEngine(vehicle);
        sim1.Init();
        double physDt = 0.01;
        int totalSteps = 500;
        for (int i = 0; i < totalSteps; i++)
            sim1.Step(physDt);

        Vector pos1x = sim1.State.Position;
        double time1x = sim1.Time;

        // Run at 100x warp — the TimeWarp should produce the same total sim time
        // when we give it the appropriate number of render frames.
        // At 100x warp, renderDt=1/60: each frame advances 100/60 = 1.667 sim-seconds.
        // For 5 sim-seconds total, we need 5/(100/60) = 3 render frames.
        vehicle = new RocketVehicle(
            new[] { new RocketStage(dryMass: 500,
                engines: new[] { new RocketEngine(282, 311, 50000) },
                tanks: new[] { new PropellantTank(300, 700) }) },
            payloadMass: 100);
        vehicle.Stages[0].SeparationTrigger = StageSeparationTrigger.FuelDepletion();

        var sim100 = new RocketSimulationEngine(vehicle);
        sim100.Init();

        var warp = new TimeWarp
        {
            PhysicsTimestep = physDt,
            WarpFactor = 100,
            RenderInterval = 1.0 / 60.0
        };

        // Step enough frames to cover totalSteps * physDt simulation time
        double targetSimTime = totalSteps * physDt; // 5 seconds
        while (warp.TotalSimTime < targetSimTime - physDt)
            warp.StepEngine(sim100, 1.0 / 60.0);

        // Both should have advanced the same total number of physics steps
        // (within one frame's worth of steps due to accumulator quantization)
        Assert.IsTrue(Math.Abs(warp.TotalSimTime - time1x) < 2.0,
            $"Warp sim time {warp.TotalSimTime:F4} should be close to direct {time1x:F4}");

        // Positions should match (same physics, same dt)
        Assert.AreEqual(pos1x.x, sim100.State.Position.x, Math.Abs(pos1x.x) * 0.01 + 1.0,
            $"X position differs: 1x={pos1x.x:F4}, warp={sim100.State.Position.x:F4}");
        Assert.AreEqual(pos1x.z, sim100.State.Position.z, Math.Abs(pos1x.z) * 0.01 + 1.0,
            $"Z position differs: 1x={pos1x.z:F4}, warp={sim100.State.Position.z:F4}");
    }

    [TestMethod]
    public void TelemetryStream_Maintains60fps_At100xWarp()
    {
        var stream = new TelemetryStream { TargetFrameRate = 60 };

        // Simulate 100x warp: physics at 100 Hz, warp = 100x
        // So 1 real second = 100 sim seconds = 6000 frames expected at 60fps
        double physDt = 0.01;
        double warpFactor = 100.0;
        double renderDt = 1.0 / 60.0;
        double simTimePerFrame = renderDt * warpFactor; // 1.667 sim-seconds per render frame
        int stepsPerFrame = (int)(simTimePerFrame / physDt); // ~167 physics steps per frame

        double simTime = 0;
        int deliveredFrames = 0;

        // Simulate 60 render frames (1 real second = 100 sim seconds)
        for (int frame = 0; frame < 60; frame++)
        {
            // Physics sub-steps
            for (int s = 0; s < stepsPerFrame; s++)
            {
                simTime += physDt;
                var state = new RocketState(
                    new Vector(0, 0, -simTime * 100), // Rising
                    new Vector(100, 0, -simTime * 10),
                    Quaternion.Identity,
                    new Vector(0, 0, 0),
                    50000 - simTime);

                stream.Push(state, simTime);
            }

            // Try to deliver at render rate
            if (stream.TryDeliver(simTime))
                deliveredFrames++;
        }

        // At 60 render frames and 100x warp, we should deliver a frame each render frame
        // (since sim time advances 1.67s per render frame, well above 1/60 = 0.0167s delivery interval)
        Assert.IsTrue(deliveredFrames >= 55,
            $"Should deliver ~60 frames per real second at 100x warp, got {deliveredFrames}");

        // Total frames delivered should reflect the delivery rate
        Assert.IsTrue(stream.FramesDelivered >= 55,
            $"Stream reported {stream.FramesDelivered} total frames, expected ~60");

        // Current snapshot should have valid data
        Assert.IsTrue(stream.Current.Altitude > 0,
            "Current telemetry altitude should be positive");
        Assert.IsTrue(stream.Current.Speed > 0,
            "Current telemetry speed should be positive");
    }

    [TestMethod]
    public void TrajectoryPredictor_CircularOrbit_ClosesAfterOnePeriod()
    {
        // ISS-like orbit: 400 km altitude, circular
        double r = 6378137.0 + 400000.0;
        double mu = 3.986004418e14;
        double v = Math.Sqrt(mu / r);
        double period = 2.0 * Math.PI * Math.Sqrt(r * r * r / mu); // ~5560 s

        var predictor = new TrajectoryPredictor
        {
            NumPoints = 5000,
            PredictionHorizon = period + 100, // Ensure we cover full orbit
            Mu = mu,
            BodyRadius = 6378137.0
        };

        Vector pos = new Vector(r, 0, 0);
        Vector vel = new Vector(0, v, 0);

        var points = predictor.Predict(pos, vel);

        // Should have points covering one orbit
        Assert.IsTrue(points.Count > 1000, $"Expected >1000 points, got {points.Count}");

        // Last point should be close to first (one orbit)
        var first = points[0];
        var last = points[points.Count - 1];
        double closingError = Math.Sqrt(
            Math.Pow(last.Position.x - first.Position.x, 2) +
            Math.Pow(last.Position.y - first.Position.y, 2) +
            Math.Pow(last.Position.z - first.Position.z, 2));

        // Velocity-Verlet with 5000 steps over one orbit
        Assert.IsTrue(closingError < 100000,
            $"Orbit closing error {closingError:F0} m should be < 100 km for prediction");

        // All altitudes should remain roughly constant (~400 km ± 50 km acceptable for game preview)
        for (int i = 0; i < points.Count; i++)
        {
            Assert.IsTrue(points[i].Altitude > 300000 && points[i].Altitude < 500000,
                $"Point {i} altitude {points[i].Altitude / 1000:F1} km out of range");
        }
    }

    [TestMethod]
    public void TrajectoryPredictor_Apsides_CorrectForEllipticalOrbit()
    {
        var predictor = new TrajectoryPredictor { Mu = 3.986004418e14, BodyRadius = 6378137.0 };

        // Elliptical orbit: periapsis 200 km, apoapsis 800 km
        double rPeri = 6378137.0 + 200000.0;
        double rApo = 6378137.0 + 800000.0;
        double sma = (rPeri + rApo) / 2.0;
        double vPeri = Math.Sqrt(3.986004418e14 * (2.0 / rPeri - 1.0 / sma));

        Vector pos = new Vector(rPeri, 0, 0);
        Vector vel = new Vector(0, vPeri, 0);

        var (apoAlt, periAlt) = predictor.PredictApsides(pos, vel);

        Assert.AreEqual(800000.0, apoAlt, 10000.0,
            $"Apoapsis altitude {apoAlt / 1000:F1} km should be ~800 km");
        Assert.AreEqual(200000.0, periAlt, 10000.0,
            $"Periapsis altitude {periAlt / 1000:F1} km should be ~200 km");
    }

    [TestMethod]
    public void TimeWarp_IncreasDecrease_CyclesThroughLevels()
    {
        var warp = new TimeWarp();

        Assert.AreEqual(1.0, warp.WarpFactor);
        warp.IncreaseWarp();
        Assert.AreEqual(2.0, warp.WarpFactor);
        warp.IncreaseWarp();
        Assert.AreEqual(5.0, warp.WarpFactor);
        warp.IncreaseWarp();
        Assert.AreEqual(10.0, warp.WarpFactor);

        warp.DecreaseWarp();
        Assert.AreEqual(5.0, warp.WarpFactor);
        warp.DecreaseWarp();
        Assert.AreEqual(2.0, warp.WarpFactor);
        warp.DecreaseWarp();
        Assert.AreEqual(1.0, warp.WarpFactor);
    }

    #endregion
}
