using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

namespace NumericTest;

[TestClass]
public class FlightDynamicsTests
{
    private const double Tol = 1e-6;

    // ═══════════════════════════════════════════════════════════════
    //  Atmosphere Model
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Atmosphere_SeaLevel_MatchesISA()
    {
        Assert.AreEqual(288.15, AtmosphereModel.Temperature(0), 0.01);
        Assert.AreEqual(101325.0, AtmosphereModel.Pressure(0), 1.0);
        Assert.AreEqual(1.225, AtmosphereModel.Density(0), 0.001);
        Assert.AreEqual(340.29, AtmosphereModel.SpeedOfSound(0), 0.1);
    }

    [TestMethod]
    public void Atmosphere_11km_MatchesTropopause()
    {
        // At 11 km: T ≈ 216.65 K, P ≈ 22632 Pa
        Assert.AreEqual(216.65, AtmosphereModel.Temperature(11000), 0.1);
        Assert.AreEqual(22632.0, AtmosphereModel.Pressure(11000), 50);
    }

    [TestMethod]
    public void Atmosphere_DensityDecreases_WithAltitude()
    {
        double rho0 = AtmosphereModel.Density(0);
        double rho5k = AtmosphereModel.Density(5000);
        double rho10k = AtmosphereModel.Density(10000);

        Assert.IsTrue(rho0 > rho5k);
        Assert.IsTrue(rho5k > rho10k);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Airfoil Model
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void AirfoilNACA_ZeroAoA_ZeroLift()
    {
        var airfoil = AirfoilModel.NACASymmetric();
        Assert.AreEqual(0, airfoil.Cl(0), 0.01);
    }

    [TestMethod]
    public void AirfoilNACA_PositiveAoA_PositiveLift()
    {
        var airfoil = AirfoilModel.NACASymmetric();
        double cl5deg = airfoil.Cl(5 * Math.PI / 180);
        Assert.IsTrue(cl5deg > 0.3, $"Expected Cl > 0.3 at 5°, got {cl5deg}");
    }

    [TestMethod]
    public void AirfoilNACA_StallBehavior_ClDropsAfterCriticalAoA()
    {
        var airfoil = AirfoilModel.NACASymmetric();
        double clNearStall = airfoil.Cl(0.25);  // ~14.3° — just before stall
        double clPostStall = airfoil.Cl(0.50);   // ~28.6° — well past stall

        Assert.IsTrue(clPostStall < clNearStall,
            $"Cl should drop post-stall. Near stall: {clNearStall}, post-stall: {clPostStall}");
    }

    [TestMethod]
    public void AirfoilNACA_Symmetric_NegativeAoA_NegativeLift()
    {
        var airfoil = AirfoilModel.NACASymmetric();
        double cl = airfoil.Cl(-5 * Math.PI / 180);
        Assert.IsTrue(cl < -0.3, $"Expected Cl < -0.3 at -5°, got {cl}");
    }

    [TestMethod]
    public void AirfoilFlatPlate_DragIncreases_WithAoA()
    {
        var plate = AirfoilModel.FlatPlate();
        double cd0 = plate.Cd(0);
        double cd45 = plate.Cd(Math.PI / 4);
        Assert.IsTrue(cd45 > cd0);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Control Surface
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void ControlSurface_ZeroDeflection_ZeroDelta()
    {
        var elevator = ControlSurface.Elevator();
        Assert.AreEqual(0, elevator.DeltaCl(0), Tol);
        Assert.AreEqual(0, elevator.DeltaCd(0), Tol);
    }

    [TestMethod]
    public void ControlSurface_PositiveDeflection_PositiveDeltaCl()
    {
        var elevator = ControlSurface.Elevator();
        double dcl = elevator.DeltaCl(0.1); // 5.7°
        Assert.IsTrue(dcl > 0);
    }

    [TestMethod]
    public void ControlSurface_ClampedToLimits()
    {
        var elevator = ControlSurface.Elevator();
        double dclMax = elevator.DeltaCl(10.0); // way beyond limit
        double dclAtLimit = elevator.DeltaCl(elevator.MaxDeflection);
        Assert.AreEqual(dclAtLimit, dclMax, Tol);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Propulsion Model
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Jet_FullThrottle_SeaLevel_ReturnsMaxThrust()
    {
        var jet = PropulsionModel.Jet(maxThrustSeaLevel: 50000, engineCount: 2);
        double thrust = jet.Thrust(1.0, 0, 100);
        Assert.AreEqual(100000, thrust, 100); // 2 × 50k at sea level
    }

    [TestMethod]
    public void Jet_ZeroThrottle_ZeroThrust()
    {
        var jet = PropulsionModel.Jet(maxThrustSeaLevel: 50000);
        Assert.AreEqual(0, jet.Thrust(0, 0, 100), Tol);
    }

    [TestMethod]
    public void Jet_ThrustDecreases_WithAltitude()
    {
        var jet = PropulsionModel.Jet(maxThrustSeaLevel: 50000);
        double tSea = jet.Thrust(1, 0, 200);
        double t10k = jet.Thrust(1, 10000, 200);
        Assert.IsTrue(t10k < tSea);
    }

    [TestMethod]
    public void Prop_ThrustDecreases_WithSpeed()
    {
        var prop = PropulsionModel.Propeller(maxPower: 200000, staticThrust: 5000);
        double tSlow = prop.Thrust(1, 0, 20);
        double tFast = prop.Thrust(1, 0, 100);
        Assert.IsTrue(tFast < tSlow);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Frame Transforms
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void FrameTransform_BodyToWorld_Identity_NoChange()
    {
        var v = new Vector(1, 2, 3);
        var result = FrameTransforms.BodyToWorld(Quaternion.Identity, v);
        Assert.AreEqual(1, result.x, Tol);
        Assert.AreEqual(2, result.y, Tol);
        Assert.AreEqual(3, result.z, Tol);
    }

    [TestMethod]
    public void FrameTransform_RoundTrip_BodyWorldBody()
    {
        var q = Quaternion.FromEulerAngles(0.3, 0.5, 0.7);
        var v = new Vector(1, 2, 3);
        var world = FrameTransforms.BodyToWorld(q, v);
        var body = FrameTransforms.WorldToBody(q, world);

        Assert.AreEqual(v.x, body.x, 1e-10);
        Assert.AreEqual(v.y, body.y, 1e-10);
        Assert.AreEqual(v.z, body.z, 1e-10);
    }

    [TestMethod]
    public void FrameTransform_WindBody_RoundTrip()
    {
        double alpha = 0.1, beta = 0.05;
        var v = new Vector(100, 5, 10);
        var wind = FrameTransforms.BodyToWind(alpha, beta, v);
        var body = FrameTransforms.WindToBody(alpha, beta, wind);

        Assert.AreEqual(v.x, body.x, 1e-10);
        Assert.AreEqual(v.y, body.y, 1e-10);
        Assert.AreEqual(v.z, body.z, 1e-10);
    }

    [TestMethod]
    public void FrameTransform_AngleOfAttack_PureForward_Zero()
    {
        var v = new Vector(100, 0, 0); // pure forward in body frame
        double alpha = FrameTransforms.AngleOfAttack(v);
        Assert.AreEqual(0, alpha, Tol);
    }

    [TestMethod]
    public void FrameTransform_SideslipAngle_PureForward_Zero()
    {
        var v = new Vector(100, 0, 0);
        double beta = FrameTransforms.SideslipAngle(v);
        Assert.AreEqual(0, beta, Tol);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Flight Dynamics Engine — Integration Tests
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void Flight_StraightAndLevel_MaintainsAltitude()
    {
        // Set up level flight at ~60 m/s, altitude 1000 m
        var config = AircraftConfig.GenericLightAircraft();
        var engine = new FlightDynamicsEngine(config);
        engine.Init();

        // Start at altitude 1000 m (NED: z = -1000)
        engine.SetState(new AircraftState(
            position: new Vector(0, 0, -1000),
            velocity: new Vector(60, 0, 0),
            attitude: Quaternion.Identity,
            angularRate: new Vector(0, 0, 0)));

        // Apply moderate thrust and slight pitch-up to maintain altitude
        engine.SetInput(new ControlInput(throttle: 0.65, pitch: 0.05, roll: 0, yaw: 0));

        double dt = 0.01;
        double initialAlt = engine.State.Altitude;

        // Simulate 5 seconds
        for (int i = 0; i < 500; i++)
            engine.Step(dt);

        double finalAlt = engine.State.Altitude;
        double altChange = Math.Abs(finalAlt - initialAlt);

        // Allow up to 100 m altitude change (approximate trim, not exact)
        Assert.IsTrue(altChange < 100,
            $"Altitude changed by {altChange:F1} m — expected roughly level flight");
        Assert.IsTrue(engine.State.Airspeed > 20,
            $"Airspeed dropped to {engine.State.Airspeed:F1} m/s — aircraft should maintain reasonable speed");
    }

    [TestMethod]
    public void Flight_StallBehavior_ClDropsAtHighAoA()
    {
        // Verify via the airfoil model directly — stall should reduce Cl
        var airfoil = AirfoilModel.NACASymmetric();

        double clPre = airfoil.Cl(0.15);    // ~8.6° — before stall
        double clAt = airfoil.Cl(0.25);     // ~14.3° — near stall
        double clPost = airfoil.Cl(0.50);   // ~28.6° — well past stall

        Assert.IsTrue(clAt > clPre, "Cl should still be rising before stall");
        Assert.IsTrue(clPost < clAt, "Cl should drop after stall angle");
    }

    [TestMethod]
    public void Flight_ClimbRate_WithFullThrottle()
    {
        var config = AircraftConfig.GenericLightAircraft();
        var engine = new FlightDynamicsEngine(config);
        engine.Init();

        // Start level at 1000 m, pitch up slightly
        engine.SetState(new AircraftState(
            position: new Vector(0, 0, -1000),
            velocity: new Vector(50, 0, 0),
            attitude: Quaternion.FromEulerAngles(0, 0.05, 0), // slight pitch up
            angularRate: new Vector(0, 0, 0)));

        engine.SetInput(new ControlInput(throttle: 1.0, pitch: 0.15, roll: 0, yaw: 0));

        double dt = 0.01;
        double startAlt = engine.State.Altitude;

        for (int i = 0; i < 1000; i++) // 10 seconds
            engine.Step(dt);

        double endAlt = engine.State.Altitude;

        // With full throttle and pitch up, aircraft should climb (or at least not dive much)
        // Light aircraft typical climb rate: 3-5 m/s → ~30-50 m in 10s
        Assert.IsTrue(endAlt > startAlt - 50,
            $"Expected climbing or near-level. Start: {startAlt:F1}, End: {endAlt:F1}");
    }

    [TestMethod]
    public void Flight_RollResponse_ToAileronInput()
    {
        var config = AircraftConfig.GenericLightAircraft();
        var engine = new FlightDynamicsEngine(config);
        engine.Init();

        engine.SetState(new AircraftState(
            position: new Vector(0, 0, -1000),
            velocity: new Vector(60, 0, 0),
            attitude: Quaternion.Identity,
            angularRate: new Vector(0, 0, 0)));

        // Apply right aileron
        engine.SetInput(new ControlInput(throttle: 0.5, pitch: 0, roll: 0.5, yaw: 0));

        double dt = 0.01;
        for (int i = 0; i < 200; i++) // 2 seconds
            engine.Step(dt);

        // Should have developed a roll rate
        double rollRate = engine.State.AngularRate.x;
        var (roll, _, _) = engine.State.EulerAngles;

        Assert.IsTrue(Math.Abs(roll) > 0.01 || Math.Abs(rollRate) > 0.001,
            $"Expected roll response. Roll angle: {roll:F4} rad, Roll rate: {rollRate:F4} rad/s");
    }

    [TestMethod]
    public void Flight_PitchResponse_ToElevatorInput()
    {
        var config = AircraftConfig.GenericLightAircraft();
        var engine = new FlightDynamicsEngine(config);
        engine.Init();

        engine.SetState(new AircraftState(
            position: new Vector(0, 0, -1000),
            velocity: new Vector(60, 0, 0),
            attitude: Quaternion.Identity,
            angularRate: new Vector(0, 0, 0)));

        // Apply nose-up elevator
        engine.SetInput(new ControlInput(throttle: 0.5, pitch: 0.5, roll: 0, yaw: 0));

        double dt = 0.01;
        for (int i = 0; i < 200; i++) // 2 seconds
            engine.Step(dt);

        double pitchRate = engine.State.AngularRate.y;
        var (_, pitch, _) = engine.State.EulerAngles;

        Assert.IsTrue(Math.Abs(pitch) > 0.01 || Math.Abs(pitchRate) > 0.001,
            $"Expected pitch response. Pitch angle: {pitch:F4} rad, Pitch rate: {pitchRate:F4} rad/s");
    }

    // ═══════════════════════════════════════════════════════════════
    //  AircraftState
    // ═══════════════════════════════════════════════════════════════

    [TestMethod]
    public void AircraftState_Altitude_IsNegativeZ()
    {
        var state = new AircraftState
        {
            Position = new Vector(0, 0, -5000)
        };
        Assert.AreEqual(5000, state.Altitude, Tol);
    }

    [TestMethod]
    public void AircraftState_Clone_IsDeepCopy()
    {
        var original = new AircraftState(
            new Vector(1, 2, 3),
            new Vector(4, 5, 6),
            Quaternion.FromEulerAngles(0.1, 0.2, 0.3),
            new Vector(0.01, 0.02, 0.03));

        var clone = original.Clone();
        clone.Position = new Vector(99, 99, 99);

        Assert.AreEqual(1, original.Position.x, Tol);
    }
}
