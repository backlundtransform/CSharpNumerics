using System;
using CSharpNumerics.Engines.Game.Rocket;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics.Propulsion;

namespace CSharpNumerics.Engines.Game.Rocket.Examples;

/// <summary>
/// Example: Falcon 9-style two-stage to LEO with booster landing.
/// 
/// Demonstrates:
/// - Two-stage vehicle configuration (9 Merlin first stage + 1 Merlin Vac second stage)
/// - Full ascent to 200 km LEO
/// - Boostback and landing burn on first stage (conceptual — uses TimeWarp for fast-forward)
/// - TelemetryStream for HUD data
/// - TrajectoryPredictor for orbital line display
/// - RocketUnityAdapter for rendering
///
/// <code>
/// // Setup (in Unity MonoBehaviour.Start):
/// var example = new Falcon9ToLEOExample();
/// example.Configure();
/// 
/// // Each frame (in Unity MonoBehaviour.Update):
/// example.UpdateFrame(Time.deltaTime);
/// Vector pos = example.Adapter.UnityPosition;
/// Quaternion rot = example.Adapter.UnityRotation;
/// transform.position = new UnityEngine.Vector3((float)pos.x, (float)pos.y, (float)pos.z);
/// 
/// // HUD data:
/// double alt = example.Telemetry.Current.Altitude;
/// double spd = example.Telemetry.Current.Speed;
/// double fuel = example.Telemetry.Current.FuelPercent;
/// double apo = example.Telemetry.Current.ApoapsisAltitude;
/// </code>
/// </summary>
public class Falcon9ToLEOExample
{
    /// <summary>The simulation engine.</summary>
    public RocketSimulationEngine Engine { get; private set; }

    /// <summary>Unity adapter for rendering.</summary>
    public RocketUnityAdapter Adapter { get; } = new RocketUnityAdapter();

    /// <summary>Telemetry stream for HUD.</summary>
    public TelemetryStream Telemetry { get; } = new TelemetryStream();

    /// <summary>Trajectory predictor for orbit line.</summary>
    public TrajectoryPredictor Predictor { get; } = new TrajectoryPredictor();

    /// <summary>Time warp controller.</summary>
    public TimeWarp Warp { get; } = new TimeWarp();

    /// <summary>Mission profile with phase sequencing.</summary>
    public MissionProfile Profile { get; private set; }

    /// <summary>Whether MECO (Main Engine Cut-Off) has occurred.</summary>
    public bool MecoOccurred { get; private set; }

    /// <summary>Whether the second stage has achieved orbit.</summary>
    public bool OrbitAchieved { get; private set; }

    /// <summary>
    /// Configures the Falcon 9-like vehicle and initializes the simulation.
    /// </summary>
    public void Configure()
    {
        // First stage: 9x Merlin 1D (sea level)
        // RocketEngine(ispSL, ispVac, maxThrustVac)
        var firstStageEngines = new RocketEngine[9];
        for (int i = 0; i < 9; i++)
            firstStageEngines[i] = new RocketEngine(282, 311, 914000);

        // PropellantTank(fuelMass, oxidizerMass) — RP-1/LOX ratio ~2.56
        var firstStageTank = new PropellantTank(111000, 284700); // ~396 tonnes total
        var firstStage = new RocketStage(
            dryMass: 25600,
            engines: firstStageEngines,
            tanks: new[] { firstStageTank });
        firstStage.SeparationTrigger = StageSeparationTrigger.FuelDepletion();

        // Second stage: 1x Merlin Vacuum
        var secondStageTank = new PropellantTank(26000, 66670); // ~93 tonnes total
        var secondStage = new RocketStage(
            dryMass: 4000,
            engines: new[] { new RocketEngine(200, 348, 934000) },
            tanks: new[] { secondStageTank });
        secondStage.SeparationTrigger = StageSeparationTrigger.FuelDepletion();

        // Vehicle
        var vehicle = new RocketVehicle(
            new[] { firstStage, secondStage },
            payloadMass: 22800); // 22.8 tonnes to LEO

        Engine = new RocketSimulationEngine(vehicle);
        Engine.UseEarthModel = true;
        Engine.LaunchSite = GeoCoordinate.FromDegrees(28.5, -80.5, 0); // Cape Canaveral

        // Mission profile — standard ascent to 200 km LEO
        Profile = MissionProfile.StandardAscentToOrbit(200000);

        // Initialize
        Engine.SetStateFromLaunchSite(Engine.LaunchSite.Value);
        Engine.Init();

        // Configure warp
        Warp.PhysicsTimestep = 0.05; // 20 Hz physics
        Warp.WarpFactor = 1;

        // Adapter in ECI mode for orbital flight
        Adapter.CoordinateMode = AdapterCoordinateMode.ECIToUnity;
        Adapter.PositionScale = 0.001; // km scale for rendering

        Telemetry.TargetFrameRate = 60;
    }

    /// <summary>
    /// Called each render frame. Advances simulation and updates all subsystems.
    /// </summary>
    /// <param name="renderDt">Render frame time (seconds). Use Time.deltaTime in Unity.</param>
    public void UpdateFrame(double renderDt)
    {
        if (Engine == null || !Engine.IsInitialized) return;

        // Advance physics
        Warp.StepEngine(Engine, renderDt);

        // Update subsystems
        RocketState state = Engine.State;
        Adapter.Update(state);

        double fuelFraction = Engine.Vehicle.AllStagesSpent ? 0 :
            Engine.Vehicle.ActiveStage.PropellantMass /
            Math.Max(1, Engine.Vehicle.ActiveStage.InitialPropellantMass);
        Telemetry.Push(state, Engine.Time, fuelFraction: fuelFraction);
        Telemetry.TryDeliver(Engine.Time);

        // Update trajectory prediction every ~1 second sim time
        if (((int)(Engine.Time * 10)) % 10 == 0)
        {
            Predictor.Predict(state.Position, state.Velocity);
        }
    }
}
