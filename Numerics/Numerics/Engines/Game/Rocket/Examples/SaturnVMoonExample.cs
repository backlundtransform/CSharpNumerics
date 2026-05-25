using System;
using CSharpNumerics.Engines.Game.Rocket;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.Mechanics.Propulsion;
using CSharpNumerics.Physics.OrbitalMechanics;

namespace CSharpNumerics.Engines.Game.Rocket.Examples;

/// <summary>
/// Example: Saturn V Moon trajectory (trans-lunar injection).
/// 
/// Demonstrates:
/// - Three-stage vehicle (S-IC, S-II, S-IVB)
/// - Earth parking orbit followed by TLI burn
/// - TrajectoryPredictor showing Moon-intercept trajectory
/// - TimeWarp at 100x+ for coast phases
/// - Telemetry showing altitude, speed, distance to Moon
///
/// <code>
/// var example = new SaturnVMoonExample();
/// example.Configure();
/// 
/// // Fast-forward through parking orbit
/// example.Warp.WarpFactor = 100;
/// 
/// // Each frame:
/// example.UpdateFrame(Time.deltaTime);
/// 
/// // Check mission status:
/// if (example.TLIComplete) { /* switch to Moon-relative view */ }
/// </code>
/// </summary>
public class SaturnVMoonExample
{
    /// <summary>The simulation engine.</summary>
    public RocketSimulationEngine Engine { get; private set; }

    /// <summary>Unity adapter for rendering.</summary>
    public RocketUnityAdapter Adapter { get; } = new RocketUnityAdapter();

    /// <summary>Telemetry stream for HUD.</summary>
    public TelemetryStream Telemetry { get; } = new TelemetryStream();

    /// <summary>Trajectory predictor.</summary>
    public TrajectoryPredictor Predictor { get; } = new TrajectoryPredictor();

    /// <summary>Time warp controller.</summary>
    public TimeWarp Warp { get; } = new TimeWarp();

    /// <summary>Whether Trans-Lunar Injection burn is complete.</summary>
    public bool TLIComplete { get; private set; }

    /// <summary>Distance to Moon (meters). Updated each frame.</summary>
    public double DistanceToMoon { get; private set; }

    /// <summary>Moon position vector (simplified circular orbit).</summary>
    public Vector MoonPosition { get; private set; }

    private const double MoonDistance = 384400000.0; // 384,400 km
    private const double MoonOrbitalPeriod = 27.3 * 86400.0; // ~27.3 days in seconds

    /// <summary>
    /// Configures the Saturn V vehicle and initializes the simulation.
    /// </summary>
    public void Configure()
    {
        // S-IC (First stage): 5x F-1 engines
        // RocketEngine(ispSL, ispVac, maxThrustVac)
        var sIcEngines = new RocketEngine[5];
        for (int i = 0; i < 5; i++)
            sIcEngines[i] = new RocketEngine(263, 304, 7770000);

        // PropellantTank(fuelMass, oxidizerMass) — RP-1/LOX ratio ~2.27
        var sIcTank = new PropellantTank(660000, 1500000); // ~2160 tonnes total
        var sIc = new RocketStage(dryMass: 130000, engines: sIcEngines, tanks: new[] { sIcTank });
        sIc.SeparationTrigger = StageSeparationTrigger.FuelDepletion();

        // S-II (Second stage): 5x J-2 engines (LH2/LOX)
        var j2Engines = new RocketEngine[5];
        for (int i = 0; i < 5; i++)
            j2Engines[i] = new RocketEngine(200, 421, 1033000);

        var sIiTank = new PropellantTank(74000, 370000); // ~444 tonnes LOX/LH2 (ratio ~5:1)
        var sIi = new RocketStage(dryMass: 36000, engines: j2Engines, tanks: new[] { sIiTank });
        sIi.SeparationTrigger = StageSeparationTrigger.FuelDepletion();

        // S-IVB (Third stage): 1x J-2 (restartable for TLI)
        var sIvbEngine = new RocketEngine(200, 421, 1033000);
        var sIvbTank = new PropellantTank(18000, 89000); // ~107 tonnes
        var sIvb = new RocketStage(dryMass: 11500, engines: new[] { sIvbEngine }, tanks: new[] { sIvbTank });
        sIvb.SeparationTrigger = StageSeparationTrigger.FuelDepletion();

        // Vehicle (payload = Apollo CSM + LM ≈ 45 tonnes)
        var vehicle = new RocketVehicle(
            new[] { sIc, sIi, sIvb },
            payloadMass: 45000);

        Engine = new RocketSimulationEngine(vehicle);
        Engine.UseEarthModel = true;
        Engine.LaunchSite = GeoCoordinate.FromDegrees(28.6, -80.6, 0); // KSC

        // Initialize
        Engine.SetStateFromLaunchSite(Engine.LaunchSite.Value);
        Engine.Init();

        // Configure subsystems
        Warp.PhysicsTimestep = 0.1;
        Warp.MaxWarpFactor = 10000; // Need high warp for Moon coast

        Adapter.CoordinateMode = AdapterCoordinateMode.ECIToUnity;
        Adapter.PositionScale = 1e-6; // Mm scale for Moon trajectory

        Predictor.PredictionHorizon = 3 * 86400; // 3 days prediction
        Predictor.NumPoints = 500;

        Telemetry.TargetFrameRate = 60;
    }

    /// <summary>
    /// Called each render frame.
    /// </summary>
    public void UpdateFrame(double renderDt)
    {
        if (Engine == null || !Engine.IsInitialized) return;

        Warp.StepEngine(Engine, renderDt);

        RocketState state = Engine.State;
        Adapter.Update(state);

        double fuelFraction = Engine.Vehicle.AllStagesSpent ? 0 :
            Engine.Vehicle.ActiveStage.PropellantMass /
            Math.Max(1, Engine.Vehicle.ActiveStage.InitialPropellantMass);
        Telemetry.Push(state, Engine.Time, fuelFraction: fuelFraction);
        Telemetry.TryDeliver(Engine.Time);

        // Update Moon position (simplified circular orbit in ECI XY plane)
        double moonAngle = 2.0 * Math.PI * Engine.Time / MoonOrbitalPeriod;
        MoonPosition = new Vector(
            MoonDistance * Math.Cos(moonAngle),
            MoonDistance * Math.Sin(moonAngle),
            0);

        Vector toMoon = new Vector(
            MoonPosition.x - state.Position.x,
            MoonPosition.y - state.Position.y,
            MoonPosition.z - state.Position.z);
        DistanceToMoon = toMoon.GetMagnitude();

        // TLI detection: vehicle beyond 50000 km from Earth with positive radial velocity
        double r = state.Position.GetMagnitude();
        if (r > 50000000 && !TLIComplete)
        {
            double radialVel = (state.Position.x * state.Velocity.x +
                                state.Position.y * state.Velocity.y +
                                state.Position.z * state.Velocity.z) / r;
            if (radialVel > 0)
                TLIComplete = true;
        }
    }
}
