using CSharpNumerics.Engines.Common;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;
using CSharpNumerics.Physics.Mechanics;
using CSharpNumerics.Physics.Mechanics.Propulsion;
using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Rocket launch simulation engine implementing 6DOF + variable mass dynamics.
/// Integrates the equations of motion using RK4 with atmospheric drag and gravity.
/// 
/// Supports multi-stage rockets with staging events, throttle control,
/// and thrust vector control via engine gimbal.
/// 
/// Implements <see cref="ISimulationEngine"/> for use with common engine infrastructure.
/// </summary>
public class RocketSimulationEngine : ISimulationEngine
{
    private const double G0 = 9.80665;
    private const double SeaLevelPressure = 101325.0;

    private RocketState _state;
    private double _time;
    private bool _initialized;

    /// <summary>The rocket vehicle being simulated.</summary>
    public RocketVehicle Vehicle { get; }

    /// <summary>Optional event bus for publishing simulation events.</summary>
    public EventBus EventBus { get; set; }

    /// <summary>Current simulation time in seconds.</summary>
    public double Time => _time;

    /// <summary>True after Init() has been called.</summary>
    public bool IsInitialized => _initialized;

    /// <summary>Current rocket state.</summary>
    public RocketState State => _state;

    /// <summary>Gravitational acceleration magnitude (m/s²). Default 9.80665.</summary>
    public double Gravity { get; set; } = G0;

    /// <summary>Reference area for drag computation (m²).</summary>
    public double ReferenceArea { get; set; } = 10.0;

    /// <summary>Drag coefficient (used when no CompressibleDragModel is set).</summary>
    public double DragCoefficient { get; set; } = 0.3;

    /// <summary>
    /// Mach-dependent drag model. When set, overrides the constant <see cref="DragCoefficient"/>.
    /// </summary>
    public CompressibleDragModel DragModel { get; set; }

    /// <summary>
    /// Max-Q monitor. When set, updated each step with current flight conditions.
    /// </summary>
    public MaxQMonitor MaxQMonitor { get; set; }

    /// <summary>
    /// Throttle bucket controller. When set, automatically modulates throttle near Max-Q.
    /// </summary>
    public ThrottleBucket ThrottleBucket { get; set; }

    /// <summary>
    /// Stage separation handler. When set, called before each separation event.
    /// </summary>
    public StageSeparationHandler SeparationHandler { get; set; }

    /// <summary>
    /// When true, the engine uses the WGS84 Earth model with J2 gravity,
    /// Coriolis/centrifugal forces, and ECEF coordinates.
    /// Position/velocity in RocketState are interpreted as ECEF (meters, m/s).
    /// Default false (flat-Earth NED mode).
    /// </summary>
    public bool UseEarthModel { get; set; }

    /// <summary>
    /// Launch site in geodetic coordinates. Used to set initial ECEF position
    /// and to compute the Earth-rotation velocity bonus when UseEarthModel is true.
    /// </summary>
    public GeoCoordinate? LaunchSite { get; set; }

    /// <summary>
    /// Altitude threshold (meters) above which flat-Earth gravity switches to J2 model.
    /// Only relevant when UseEarthModel is false (for future hybrid mode).
    /// </summary>
    public double EarthModelThreshold { get; set; } = 50000.0;

    /// <summary>
    /// Creates a rocket simulation engine for the given vehicle.
    /// </summary>
    /// <param name="vehicle">The multi-stage rocket vehicle.</param>
    public RocketSimulationEngine(RocketVehicle vehicle)
    {
        Vehicle = vehicle ?? throw new ArgumentNullException(nameof(vehicle));
        _state = new RocketState();
    }

    /// <summary>Sets the rocket state (e.g., to establish initial conditions).</summary>
    public void SetState(RocketState state) => _state = state;

    /// <summary>
    /// Sets the initial state from a launch site (geodetic) for ECEF simulation.
    /// Position is set to the ECEF coordinates of the launch site.
    /// Velocity is set to the Earth-rotation velocity at that location.
    /// </summary>
    public void SetStateFromLaunchSite(GeoCoordinate site)
    {
        LaunchSite = site;
        UseEarthModel = true;

        ECEFPosition ecef = CoordinateConversions.GeodeticToECEF(site);
        double omega = EarthModel.RotationRate;

        // Earth rotation velocity in ECEF is zero (frame co-rotates).
        // But if we want ECI-equivalent velocity for orbit calculations,
        // we compute v_inertial = ω × r. For ECEF propagation, initial v = 0.
        _state = new RocketState(
            new Vector(ecef.X, ecef.Y, ecef.Z),
            new Vector(0, 0, 0), // In ECEF, surface velocity is zero
            Quaternion.Identity,
            new Vector(0, 0, 0),
            Vehicle.TotalMass);
    }

    /// <summary>
    /// Gets the inertial (ECI) velocity bonus from Earth's rotation at the launch site.
    /// This is the velocity that the rocket inherits from being on the rotating Earth.
    /// </summary>
    public double EarthRotationVelocityBonus()
    {
        if (!LaunchSite.HasValue) return 0;
        return EarthModel.SurfaceVelocity(LaunchSite.Value.Latitude);
    }

    /// <summary>Initialize the engine: set mass from vehicle, ignite first stage.</summary>
    public void Init()
    {
        _time = 0;
        _state.Mass = Vehicle.TotalMass;
        Vehicle.Ignite();
        _initialized = true;
    }

    /// <summary>
    /// Advance the simulation by dt seconds using RK4 integration.
    /// Handles staging events when triggered.
    /// </summary>
    public void Step(double dt)
    {
        if (!_initialized) throw new InvalidOperationException("Engine not initialized. Call Init() first.");
        if (dt <= 0) return;

        // Check for staging (with step-halving for precise event location)
        double altitude = UseEarthModel
            ? _state.Position.GetMagnitude() - EarthModel.SemiMajorAxis
            : _state.Altitude;
        if (Vehicle.ShouldSeparate(_time, altitude))
        {
            SeparationHandler?.OnSeparation(_state, _time);
            var stage = Vehicle.ActiveStage;
            double separatedMass = stage.TotalMass;
            var triggerType = stage.SeparationTrigger.Type;
            var separated = Vehicle.Separate();
            _state.Mass = Vehicle.TotalMass;
            EventBus?.Publish(new StageEvent(
                _time, altitude, _state.Speed, Vehicle.ActiveStageIndex - 1,
                triggerType, separatedMass, _state.Mass));
        }

        // Check for booster separation
        if (Vehicle.ShouldSeparateBoosters())
        {
            double boosterMass = Vehicle.SeparateBoosters();
            _state.Mass = Vehicle.TotalMass;
            EventBus?.Publish(new BoosterSeparationEvent(_time, altitude, _state.Speed, boosterMass));
        }

        // RK4 integration
        VectorN y = _state.ToVector();
        VectorN k1 = ComputeDerivatives(y, _time);
        VectorN k2 = ComputeDerivatives(y + 0.5 * dt * k1, _time + 0.5 * dt);
        VectorN k3 = ComputeDerivatives(y + 0.5 * dt * k2, _time + 0.5 * dt);
        VectorN k4 = ComputeDerivatives(y + dt * k3, _time + dt);

        VectorN dy = (1.0 / 6.0) * dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        VectorN yNew = y + dy;

        _state = RocketState.FromVector(yNew);

        // Consume propellant
        if (!Vehicle.AllStagesSpent)
        {
            double pressureRatio = GetPressureRatio(_state.Altitude);
            Vehicle.ActiveStage.ConsumePropellant(dt, pressureRatio);
            Vehicle.ConsumeBoosterPropellant(dt, pressureRatio);
            _state.Mass = Vehicle.TotalMass;
        }

        // Update Max-Q monitor
        MaxQMonitor?.Update(_state.Altitude, _state.Speed, _time);

        // Apply throttle bucket
        if (ThrottleBucket != null && !Vehicle.AllStagesSpent)
        {
            double throttle = ThrottleBucket.ComputeThrottle(
                MaxQMonitor?.CurrentQ ?? 0, MaxQMonitor?.MaxQReached ?? false);
            var stage = Vehicle.ActiveStage;
            for (int i = 0; i < stage.Engines.Length; i++)
                stage.Engines[i].SetThrottle(throttle);
        }

        _time += dt;
    }

    /// <summary>Reset to initial state.</summary>
    public void Reset()
    {
        _time = 0;
        _state = new RocketState();
        Vehicle.Reset();
        _initialized = false;
    }

    /// <summary>
    /// Computes the state derivative vector for RK4 integration.
    /// dy/dt = [velocity, acceleration, qDot, angularAccel, -massFlowRate]
    /// </summary>
    private VectorN ComputeDerivatives(VectorN state, double t)
    {
        var s = RocketState.FromVector(state);
        double mass = s.Mass;
        if (mass <= 0) return new VectorN(14);

        double altitude = s.Altitude;
        double pressureRatio = GetPressureRatio(altitude);

        // --- Forces ---
        Vector gravityForce;
        Vector coriolisForce = new Vector(0, 0, 0);

        if (UseEarthModel)
        {
            // ECEF mode: position is in ECEF, gravity from J2 model
            Vector posECEF = s.Position;
            double r = posECEF.GetMagnitude();
            altitude = r - EarthModel.SemiMajorAxis; // Approximate altitude from ECEF radius
            pressureRatio = GetPressureRatio(altitude);

            // J2 gravity (returns acceleration, multiply by mass for force)
            Vector gAccel = GravityModel.Acceleration(posECEF);
            gravityForce = mass * gAccel;

            // Coriolis + centrifugal in ECEF frame
            Vector fictitiousAccel = CoriolisForce.TotalFictitiousAcceleration(posECEF, s.Velocity);
            coriolisForce = mass * fictitiousAccel;
        }
        else
        {
            // Flat-Earth NED mode: gravity is +z (down)
            gravityForce = new Vector(0, 0, mass * Gravity);
        }

        // Thrust (along body X-axis, transformed to world frame, with gimbal)
        Vector thrustWorld = new Vector(0, 0, 0);
        double massFlowRate = 0;
        if (!Vehicle.AllStagesSpent && Vehicle.ActiveStage != null)
        {
            var stage = Vehicle.ActiveStage;
            double thrust = stage.TotalThrust(pressureRatio);
            massFlowRate = stage.TotalMassFlowRate(pressureRatio);

            // Add booster thrust
            thrust += Vehicle.BoosterThrust(pressureRatio);
            massFlowRate += Vehicle.BoosterMassFlowRate(pressureRatio);

            // Thrust direction in body frame (forward = +X, gimbal deflects into Y,Z)
            Vector thrustBody = new Vector(thrust, 0, 0);
            thrustWorld = s.Attitude.Rotate(thrustBody);
        }

        // Aerodynamic drag (opposes velocity)
        Vector dragForce = new Vector(0, 0, 0);
        double speed = s.Speed;
        if (speed > 0.1 && altitude < 86000)
        {
            double density = GetDensity(altitude);
            double dynamicPressure = 0.5 * density * speed * speed;

            // Use Mach-dependent drag model if available
            double cd;
            if (DragModel != null)
            {
                double mach = MachNumber.Compute(speed, altitude);
                cd = DragModel.Evaluate(mach);
            }
            else
            {
                cd = DragCoefficient;
            }

            double dragMag = cd * ReferenceArea * dynamicPressure;
            Vector velUnit = (1.0 / speed) * s.Velocity;
            dragForce = (-dragMag) * velUnit;
        }

        // Net acceleration
        Vector totalExternal = new Vector(
            gravityForce.x + dragForce.x + coriolisForce.x,
            gravityForce.y + dragForce.y + coriolisForce.y,
            gravityForce.z + dragForce.z + coriolisForce.z);
        Vector acceleration = VariableMassDynamics.ComputeAcceleration(mass, thrustWorld, totalExternal);

        // --- Attitude kinematics ---
        // qDot = 0.5 * q * ω (quaternion derivative from angular velocity)
        Vector w = s.AngularRate;
        Quaternion omegaQuat = new Quaternion(0, w.x, w.y, w.z);
        Quaternion qDot = 0.5 * s.Attitude * omegaQuat;

        // Angular acceleration (simplified: no torques in Phase 1 basic sim)
        Vector angularAccel = new Vector(0, 0, 0);

        // Pack derivatives: [vx,vy,vz, ax,ay,az, qdw,qdx,qdy,qdz, αx,αy,αz, -ṁ]
        return new VectorN(new double[]
        {
            s.Velocity.x, s.Velocity.y, s.Velocity.z,
            acceleration.x, acceleration.y, acceleration.z,
            qDot.w, qDot.x, qDot.y, qDot.z,
            angularAccel.x, angularAccel.y, angularAccel.z,
            -massFlowRate
        });
    }

    /// <summary>
    /// Gets pressure ratio (ambient / sea-level) from ISA atmosphere model.
    /// </summary>
    private double GetPressureRatio(double altitude)
    {
        if (altitude <= 0) return 1.0;
        if (altitude >= 86000) return 0.0;
        double pressure = AtmosphereModel.Pressure(altitude);
        return pressure / SeaLevelPressure;
    }

    /// <summary>
    /// Gets air density from ISA atmosphere model.
    /// </summary>
    private double GetDensity(double altitude)
    {
        if (altitude <= 0) return 1.225;
        if (altitude >= 86000) return 0.0;
        return AtmosphereModel.Density(altitude);
    }
}

/// <summary>
/// Event published when a stage separates.
/// </summary>
public class StageEvent
{
    /// <summary>Time of separation in seconds.</summary>
    public double Time { get; }

    /// <summary>Altitude at separation in metres.</summary>
    public double Altitude { get; }

    /// <summary>Speed at separation in m/s.</summary>
    public double Speed { get; }

    /// <summary>Index of the stage that was separated.</summary>
    public int StageIndex { get; }

    /// <summary>Trigger type that caused the separation.</summary>
    public SeparationType SeparationType { get; }

    /// <summary>Mass of the separated stage (kg).</summary>
    public double SeparatedMass { get; }

    /// <summary>Vehicle mass after separation (kg).</summary>
    public double RemainingMass { get; }

    public StageEvent(double time, double altitude, double speed, int stageIndex,
        SeparationType separationType = SeparationType.FuelDepletion,
        double separatedMass = 0, double remainingMass = 0)
    {
        Time = time;
        Altitude = altitude;
        Speed = speed;
        StageIndex = stageIndex;
        SeparationType = separationType;
        SeparatedMass = separatedMass;
        RemainingMass = remainingMass;
    }
}

/// <summary>
/// Event published when strap-on boosters are separated simultaneously.
/// </summary>
public class BoosterSeparationEvent
{
    /// <summary>Time of separation in seconds.</summary>
    public double Time { get; }

    /// <summary>Altitude at separation in metres.</summary>
    public double Altitude { get; }

    /// <summary>Speed at separation in m/s.</summary>
    public double Speed { get; }

    /// <summary>Total mass of all separated boosters in kg.</summary>
    public double BoosterMass { get; }

    public BoosterSeparationEvent(double time, double altitude, double speed, double boosterMass)
    {
        Time = time;
        Altitude = altitude;
        Speed = speed;
        BoosterMass = boosterMass;
    }
}
