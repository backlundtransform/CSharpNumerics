using CSharpNumerics.Engines.Common;
using CSharpNumerics.Engines.Game.Fluids;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;
using System;

namespace CSharpNumerics.Engines.Game.Flight;

/// <summary>
/// 6-degree-of-freedom flight dynamics engine.
/// Integrates the nonlinear equations of motion for a rigid aircraft using RK4.
/// 
/// State: position (NED), velocity (NED), attitude quaternion, body angular rates (p,q,r).
/// Forces: gravity, aerodynamic lift/drag/side, thrust, control surface moments.
/// 
/// Implements <see cref="ISimulationEngine"/> for use with the common engine infrastructure.
/// </summary>
public class FlightDynamicsEngine : ISimulationEngine
{
    private AircraftConfig _config;
    private AircraftState _state;
    private ControlInput _input;
    private double _time;
    private bool _initialized;

    /// <summary>Current simulation time in seconds.</summary>
    public double Time => _time;

    /// <summary>True after Init() has been called.</summary>
    public bool IsInitialized => _initialized;

    /// <summary>Current aircraft state (read-only reference; use SetState to replace).</summary>
    public AircraftState State => _state;

    /// <summary>Current control input.</summary>
    public ControlInput Input => _input;

    /// <summary>Aircraft configuration.</summary>
    public AircraftConfig Config => _config;

    /// <summary>
    /// Optional fluid solver for wind-field coupling.
    /// When set, the aircraft samples wind from the fluid field and
    /// exhaust is fed back as a fluid emitter.
    /// </summary>
    public GameFluidSolver3D FluidField { get; set; }

    /// <summary>
    /// Wind velocity at the aircraft's current position, sampled from the fluid field.
    /// Zero if no fluid field is attached.
    /// </summary>
    public Vector WindVelocity { get; private set; } = new Vector(0, 0, 0);

    /// <summary>
    /// Optional exhaust emitter injected into the fluid field behind the aircraft.
    /// Created automatically when a fluid field is attached.
    /// </summary>
    public FluidEmitter ExhaustEmitter { get; private set; }

    /// <summary>
    /// Creates a flight dynamics engine with the given aircraft configuration.
    /// </summary>
    public FlightDynamicsEngine(AircraftConfig config)
    {
        _config = config;
        _state = new AircraftState();
        _input = new ControlInput();
    }

    /// <summary>Sets the aircraft state (e.g., to establish initial conditions).</summary>
    public void SetState(AircraftState state) => _state = state;

    /// <summary>Sets the control input for the next step.</summary>
    public void SetInput(ControlInput input) => _input = input;

    /// <summary>Replaces the aircraft configuration.</summary>
    public void SetConfig(AircraftConfig config) => _config = config;

    /// <summary>Initialize the engine for simulation.</summary>
    public void Init()
    {
        _time = 0;
        _initialized = true;
    }

    /// <summary>Reset to initial state.</summary>
    public void Reset()
    {
        _state = new AircraftState();
        _input = new ControlInput();
        _time = 0;
    }

    /// <summary>
    /// Advance the simulation by dt seconds using classical RK4 integration.
    /// </summary>
    public void Step(double dt)
    {
        if (!_initialized)
            throw new InvalidOperationException("Call Init() before stepping.");

        // Sample wind from fluid field if attached
        if (FluidField != null)
        {
            var pos = _state.Position;
            var (wu, wv, ww) = FluidField.SampleVelocity(pos.x, pos.y, pos.z);
            WindVelocity = new Vector(wu, wv, ww);
        }
        else
        {
            WindVelocity = new Vector(0, 0, 0);
        }

        // RK4 integration on the combined state
        var s0 = _state;
        var k1 = ComputeDerivatives(s0);
        var s1 = Advance(s0, k1, dt * 0.5);
        var k2 = ComputeDerivatives(s1);
        var s2 = Advance(s0, k2, dt * 0.5);
        var k3 = ComputeDerivatives(s2);
        var s3 = Advance(s0, k3, dt);
        var k4 = ComputeDerivatives(s3);

        // Combine: state += (dt/6)(k1 + 2k2 + 2k3 + k4)
        _state = RK4Combine(s0, k1, k2, k3, k4, dt);
        _state.Attitude = _state.Attitude.Normalize();
        _time += dt;

        // Feed exhaust back into fluid field
        UpdateExhaustEmitter();
    }

    // ═══════════════════════════════════════════════════════════════
    //  Equations of motion
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Derivatives of the 12-state vector: (dPos, dVel, dQuat, dOmega).
    /// </summary>
    private StateDerivatives ComputeDerivatives(AircraftState s)
    {
        var cfg = _config;
        double mass = cfg.Mass;

        // --- Body-frame velocity (airspeed = velocity relative to wind) ---
        var aeroVelocityWorld = s.Velocity - WindVelocity;
        var vBody = FrameTransforms.WorldToBody(s.Attitude, aeroVelocityWorld);
        double u = vBody.x, v = vBody.y, w = vBody.z;
        double airspeed = aeroVelocityWorld.GetMagnitude();

        // --- Aerodynamic angles ---
        double alpha = airspeed > 0.1 ? Math.Atan2(w, u) : 0.0;
        double beta = airspeed > 0.1 ? Math.Asin(Math.Clamp(v / airspeed, -1, 1)) : 0.0;

        // --- Atmosphere ---
        double alt = s.Altitude;
        double rho = AtmosphereModel.Density(alt);
        double qBar = 0.5 * rho * airspeed * airspeed; // dynamic pressure

        // --- Lift and drag coefficients ---
        double cl = cfg.Airfoil.Cl(alpha);
        double cd = cfg.Airfoil.Cd(alpha);

        // Control surface increments
        double elevatorDeflection = _input.Pitch * cfg.Elevator.MaxDeflection;
        double aileronDeflection = _input.Roll * cfg.Aileron.MaxDeflection;
        double rudderDeflection = _input.Yaw * cfg.Rudder.MaxDeflection;

        cl += cfg.Elevator.DeltaCl(elevatorDeflection);
        cd += cfg.Elevator.DeltaCd(elevatorDeflection);
        cd += cfg.Aileron.DeltaCd(aileronDeflection);
        cd += cfg.Rudder.DeltaCd(rudderDeflection);

        // Flaps
        cl += cfg.FlapClIncrement * _input.Flaps;
        cd += cfg.FlapCdIncrement * _input.Flaps;

        // Gear drag
        if (_input.GearDown)
            cd += cfg.GearDragArea / cfg.WingArea;

        // --- Aerodynamic forces in wind frame ---
        double S = cfg.WingArea;
        double lift = qBar * S * cl;
        double drag = qBar * S * cd;
        double sideForce = 0; // simplified — no lateral aero model yet

        // Wind frame: x_w = along velocity (drag opposes), z_w = perpendicular in plane of symmetry (lift)
        // Convert to body frame via AoA rotation
        double fxAero = -drag * Math.Cos(alpha) + lift * Math.Sin(alpha);
        double fzAero = -drag * Math.Sin(alpha) - lift * Math.Cos(alpha);
        double fyAero = sideForce;

        var aeroForceBody = new Vector(fxAero, fyAero, fzAero);

        // --- Thrust (along body X) ---
        double thrust = cfg.Propulsion.Thrust(_input.Throttle, alt, airspeed);
        var thrustForceBody = new Vector(thrust, 0, 0);

        // --- Gravity in body frame ---
        var gravityWorld = new Vector(0, 0, 9.80665); // NED: +z is down
        var gravityBody = FrameTransforms.WorldToBody(s.Attitude, gravityWorld);
        var weightForceBody = mass * gravityBody;

        // --- Total force in body frame → acceleration ---
        var totalForceBody = aeroForceBody + thrustForceBody + weightForceBody;
        var accelBody = totalForceBody / mass;

        // --- Acceleration in world frame ---
        var accelWorld = FrameTransforms.BodyToWorld(s.Attitude, accelBody);

        // --- Aerodynamic moments in body frame ---
        double p = s.AngularRate.x; // roll rate
        double q = s.AngularRate.y; // pitch rate
        double r = s.AngularRate.z; // yaw rate
        double b = cfg.WingSpan;
        double c = cfg.MeanChord;

        // Non-dimensional angular rates for damping
        double pHat = airspeed > 0.1 ? (p * b) / (2.0 * airspeed) : 0;
        double qHat = airspeed > 0.1 ? (q * c) / (2.0 * airspeed) : 0;
        double rHat = airspeed > 0.1 ? (r * b) / (2.0 * airspeed) : 0;

        // Roll moment: aileron + damping
        double Cl_moment = cfg.Aileron.DeltaCl(aileronDeflection) * 0.1 + cfg.Clp * pHat;
        double rollMoment = qBar * S * b * Cl_moment;

        // Pitch moment: elevator + damping
        double Cm = cfg.Elevator.DeltaCl(elevatorDeflection) * 0.25 + cfg.Cmq * qHat;
        double pitchMoment = qBar * S * c * Cm;

        // Yaw moment: rudder + damping
        double Cn = cfg.Rudder.DeltaCl(rudderDeflection) * 0.1 + cfg.Cnr * rHat;
        double yawMoment = qBar * S * b * Cn;

        var moment = new Vector(rollMoment, pitchMoment, yawMoment);

        // --- Angular acceleration (Euler's equation for rigid body) ---
        // I·ω̇ = M - ω × (I·ω)
        // Simplified for diagonal inertia tensor (Ixz ≈ 0 for now)
        double Ixx = cfg.Ixx, Iyy = cfg.Iyy, Izz = cfg.Izz;

        double pdot = (moment.x - (Izz - Iyy) * q * r) / Ixx;
        double qdot = (moment.y - (Ixx - Izz) * p * r) / Iyy;
        double rdot = (moment.z - (Iyy - Ixx) * p * q) / Izz;

        var angularAccel = new Vector(pdot, qdot, rdot);

        // --- Quaternion derivative ---
        // q̇ = ½ · ω_quat · q  where ω_quat = (0, p, q, r)
        var omegaQ = new Quaternion(0, p, q, r);
        var qDot = 0.5 * (omegaQ * s.Attitude);

        return new StateDerivatives
        {
            DVelocity = accelWorld,
            DAngularRate = angularAccel,
            QDot = qDot
        };
    }

    /// <summary>
    /// Advances a state by adding dt * derivatives.
    /// </summary>
    private static AircraftState Advance(AircraftState s, StateDerivatives d, double dt)
    {
        var q = new Quaternion(
            s.Attitude.w + d.QDot.w * dt,
            s.Attitude.x + d.QDot.x * dt,
            s.Attitude.y + d.QDot.y * dt,
            s.Attitude.z + d.QDot.z * dt).Normalize();

        return new AircraftState(
            position: s.Position + dt * s.Velocity,
            velocity: s.Velocity + dt * d.DVelocity,
            attitude: q,
            angularRate: s.AngularRate + dt * d.DAngularRate);
    }

    /// <summary>
    /// RK4 combination: s + (dt/6)(k1 + 2k2 + 2k3 + k4).
    /// </summary>
    private static AircraftState RK4Combine(
        AircraftState s0,
        StateDerivatives k1, StateDerivatives k2,
        StateDerivatives k3, StateDerivatives k4,
        double dt)
    {
        double h6 = dt / 6.0;

        var vel = s0.Velocity + h6 * (k1.DVelocity + 2 * k2.DVelocity + 2 * k3.DVelocity + k4.DVelocity);
        var omega = s0.AngularRate + h6 * (k1.DAngularRate + 2 * k2.DAngularRate + 2 * k3.DAngularRate + k4.DAngularRate);

        // For position, use midpoint velocity
        var pos = s0.Position + dt * s0.Velocity + 0.5 * dt * h6 *
            (k1.DVelocity + 2 * k2.DVelocity + 2 * k3.DVelocity + k4.DVelocity);

        // Quaternion: integrate using combined derivative
        var qDotAvg = (1.0 / 6.0) * (k1.QDot + 2.0 * k2.QDot + 2.0 * k3.QDot + k4.QDot);
        var q = new Quaternion(
            s0.Attitude.w + qDotAvg.w * dt,
            s0.Attitude.x + qDotAvg.x * dt,
            s0.Attitude.y + qDotAvg.y * dt,
            s0.Attitude.z + qDotAvg.z * dt).Normalize();

        return new AircraftState(pos, vel, q, omega);
    }

    /// <summary>
    /// Internal struct holding time derivatives of the aircraft state.
    /// </summary>
    private struct StateDerivatives
    {
        public Vector DVelocity;       // world-frame acceleration
        public Vector DAngularRate;    // body-frame angular acceleration
        public Quaternion QDot;        // quaternion time derivative
    }

    // ═══════════════════════════════════════════════════════════════
    //  Flight ↔ Fluid integration
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Updates the exhaust emitter position and velocity in the fluid field.
    /// Called at the end of each Step() when a fluid field is attached.
    /// </summary>
    private void UpdateExhaustEmitter()
    {
        if (FluidField == null) return;

        // Position exhaust emitter behind the aircraft (negative body X direction)
        var exhaustOffsetBody = new Vector(-3, 0, 0); // 3 m behind CG
        var exhaustOffsetWorld = FrameTransforms.BodyToWorld(_state.Attitude, exhaustOffsetBody);
        var exhaustPos = _state.Position + exhaustOffsetWorld;

        // Exhaust velocity: rearward relative to aircraft, scaled by throttle
        var exhaustDirBody = new Vector(-1, 0, 0);
        var exhaustDirWorld = FrameTransforms.BodyToWorld(_state.Attitude, exhaustDirBody);
        double exhaustSpeed = 20.0 * _input.Throttle; // simplified exhaust velocity

        if (ExhaustEmitter == null)
        {
            ExhaustEmitter = new FluidEmitter(exhaustPos, densityRate: 3.0, radius: 2.0)
            {
                Velocity = exhaustSpeed * exhaustDirWorld,
                Temperature = 600
            };
            FluidField.AddEmitter(ExhaustEmitter);
        }
        else
        {
            ExhaustEmitter.Position = exhaustPos;
            ExhaustEmitter.Velocity = exhaustSpeed * exhaustDirWorld;
            ExhaustEmitter.Active = _input.Throttle > 0.01;
        }
    }
}
