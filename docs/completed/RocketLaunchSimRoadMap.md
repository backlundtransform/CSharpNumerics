# Rocket Launch Simulation Roadmap

## Goal

Build a high-fidelity rocket launch simulation engine capable of modeling all flight phases from pad ignition through orbital insertion. Target use case: Unity asset for developers building rocket/space simulation games with physically accurate trajectory, staging, and orbital mechanics.

## Architecture

The rocket simulation orchestrates existing infrastructure without duplicating physics or math:

```
Engines.Game.Rocket       → Physics.FluidDynamics.Aerodynamics (atmosphere, drag)
                          → Physics.Mechanics (variable-mass dynamics, gravity)
                          → Physics.OrbitalMechanics (Kepler, propagation)
                          → Numerics.FiniteDifference.TimeStepping (RK4, adaptive)
                          → Engines.Common (EventBus, SimulationClock)
```

### Flight Phase Architecture

```
┌─────────────────────────────────────────────────────────┐
│  RocketSimulationEngine (ISimulationEngine)              │
├─────────────────────────────────────────────────────────┤
│  Phase 1: Atmospheric Ascent (0–86 km)                  │
│    • ISA atmosphere density/pressure                     │
│    • Mach-dependent aerodynamic drag                     │
│    • Gravity turn guidance                               │
│    • Max-Q monitoring                                    │
├─────────────────────────────────────────────────────────┤
│  Phase 2: Upper Atmosphere / Vacuum (86 km+)            │
│    • Vacuum thrust (Isp_vac)                             │
│    • Inverse-square gravity                              │
│    • Stage separation events                             │
├─────────────────────────────────────────────────────────┤
│  Phase 3: Orbital Insertion                              │
│    • Orbital element computation                         │
│    • Circularization burn                                │
│    • Propagation in ECI frame                            │
└─────────────────────────────────────────────────────────┘
```

---

## Phase 1 — Rocket Propulsion & Variable Mass

Core rocket physics: engines, fuel, and the Tsiolkovsky equation.

### Physics additions (`Physics.Mechanics.Propulsion`)
- [x] `RocketEngine` — thrust (sea-level + vacuum Isp), mass flow rate ṁ = T / (Isp · g₀), throttle range, gimbal angle limits
- [x] `ThrustCurve` — thrust vs. time profile (solid boosters: grain geometry; liquid: throttle schedule)
- [x] `PropellantTank` — fuel mass, oxidizer mass, mixture ratio, remaining ΔV calculation
- [x] `VariableMassDynamics` — extends rigid-body EOM with dm/dt term: F = ma + ṁ·v_exhaust
- [x] `CenterOfMassTracker` — CG shift as fuel burns (affects moment arms and stability)
- [x] `TsiolkovskyExtensions` — ΔV = Isp·g₀·ln(m₀/mf), burn time, mass ratio utilities

### Engine additions (`Engines.Game.Rocket`)
- [x] `RocketStage` — dry mass, engines[], propellant tanks[], separation trigger (time/altitude/fuel depletion)
- [x] `RocketVehicle` — stack of RocketStages, active stage index, total mass, composite CG
- [x] `RocketState` — 13-state vector: position(3), velocity(3), quaternion(4), angular rates(3) + mass scalar
- [x] `RocketSimulationEngine` — implements `ISimulationEngine`; integrates 6DOF + variable mass using RK4

### Tests
- [x] Tsiolkovsky ΔV matches analytical: single-stage to orbit requires ~9.4 km/s
- [x] Constant thrust + constant mass-flow → linear mass decrease over burn
- [x] CG shift is monotonic and bounded between full/empty positions
- [x] RocketSimulationEngine step conserves momentum (thrust = rate of momentum change)

---

## Phase 2 — Supersonic Aerodynamics & Max-Q

Extend atmosphere model for rocket-specific high-speed aerodynamics.

### Physics additions (`Physics.FluidDynamics.Aerodynamics`)
- [x] `MachNumber` — M = v / a(h), where a from AtmosphereModel.SpeedOfSound(altitude)
- [x] `CompressibleDragModel` — Cd(Mach): subsonic → transonic rise → supersonic decay curve (Prandtl-Glauert, wave drag)
- [x] `DynamicPressure` — q = ½ρv², Max-Q detection (dq/dt = 0)
- [x] `HeatFlux` — stagnation point heating: q̇ ∝ ρ^0.5 · v^3 (Sutton-Graves approximation)

### Engine additions (`Engines.Game.Rocket`)
- [x] `MaxQMonitor` — tracks dynamic pressure, fires event at peak, provides structural load factor
- [x] `ThrottleBucket` — automatic throttle reduction near Max-Q (throttle profile optimization)
- [x] Integrate `CompressibleDragModel` into `RocketSimulationEngine` force accumulation

### Tests
- [x] Cd peaks near Mach 1.0–1.2, decays above Mach 2
- [x] Max-Q occurs at expected altitude (~11–14 km for typical ascent profiles)
- [x] Throttle bucket reduces peak dynamic pressure by configured margin
- [x] Heat flux matches order-of-magnitude for LEO insertion velocity

---

## Phase 3 — Staging & Event System

Multi-stage separation logic with discontinuity handling.

### Engine additions (`Engines.Game.Rocket`)
- [x] `StageEvent` — data class: time, altitude, velocity, stage index, separation type
- [x] `StageSeparationHandler` — on trigger: deactivate current stage, drop mass, activate next stage, publish event
- [x] `EventDetector` — monitors conditions (fuel depletion, altitude threshold, timer) and triggers step-halving for precise event location
- [x] `UllageMotor` — small impulse after separation to settle propellant before main engine ignition
- [x] `CoastPhase` — ballistic coast between stages (no thrust, only gravity + drag)

### Integration
- [x] EventBus publishes `StageEvent` — listeners can trigger audio, VFX, telemetry
- [x] AdaptiveRK45 resets after discontinuity (new initial conditions from post-separation state)
- [x] Support parallel staging (strap-on boosters separate simultaneously)

### Tests
- [x] Two-stage rocket reaches higher velocity than single-stage with same total mass
- [x] Stage separation fires at correct fuel-depletion moment (±0.1s)
- [x] Post-separation state is continuous in position/velocity, discontinuous in mass
- [x] Parallel booster separation maintains symmetric trajectory

---

## Phase 4 — Coordinate Systems & Earth Model

Required for realistic trajectory beyond flat-Earth approximation.

### Numerics additions (`Numerics.Objects`)
- [x] `GeoCoordinate` — latitude, longitude, altitude (geodetic)
- [x] `ECEFPosition` — Earth-Centered Earth-Fixed Cartesian
- [x] `ECIPosition` — Earth-Centered Inertial Cartesian
- [x] `CoordinateConversions` — geodetic ↔ ECEF ↔ ECI transforms (WGS84 ellipsoid)

### Physics additions (`Physics.Mechanics`)
- [x] `EarthModel` — WGS84 parameters (semi-major axis, flattening, rotation rate)
- [x] `GravityModel` — J2 oblateness perturbation: g(r,φ) with latitude correction
- [x] `CoriolisForce` — 2m(ω × v) for Earth-rotating frame effects

### Engine integration
- [x] `RocketSimulationEngine` switches from flat-NED to ECEF/ECI when altitude > configurable threshold (~50 km)
- [x] Launch site position (lat/lon) determines initial ECEF state and Earth-rotation velocity bonus
- [x] Gravity vector computed from J2 model, not constant downward

### Tests
- [x] Geodetic → ECEF → geodetic round-trip preserves coordinates to < 1mm
- [x] Equatorial launch gains ~460 m/s from Earth rotation
- [x] J2 perturbation produces measurable orbital precession over multiple orbits
- [x] Coriolis deflects a vertical launch (matches analytical prediction)

---

## Phase 5 — Orbital Mechanics & Insertion

Full orbital propagation for post-insertion simulation.

### Physics additions (`Physics.OrbitalMechanics`)
- [x] `OrbitalElements` — a, e, i, Ω, ω, ν (classical Keplerian elements)
- [x] `StateToElements` — convert position + velocity → orbital elements
- [x] `ElementsToState` — convert orbital elements → position + velocity
- [x] `OrbitalPropagator` — numerical integration in ECI with J2, atmospheric drag (if perigee < 400 km)
- [x] `HohmannTransfer` — ΔV₁, ΔV₂, transfer time for circular orbit changes
- [x] `OrbitalManeuver` — generalized impulsive burn (direction, magnitude, timing)
- [x] `CircularizationBurn` — compute ΔV needed at apoapsis to circularize

### Engine additions (`Engines.Game.Rocket`)
- [x] `OrbitalInsertionDetector` — detects when trajectory becomes closed orbit (eccentricity < 1, perigee > atmosphere)
- [x] `MissionProfile` — sequence of phases: ascent → gravity turn → coast → circularize → orbit
- [x] `TelemetryRecorder` — logs state vector, orbital elements, events at configurable rate

### Tests
- [x] Circular orbit at 400 km: period ≈ 92.4 min, velocity ≈ 7.67 km/s
- [x] Hohmann transfer ΔV matches analytical (LEO → GEO ≈ 3.94 km/s total)
- [x] State ↔ Elements conversion round-trip error < 1e-10
- [x] Propagated orbit remains stable over 100 revolutions (energy drift < 0.01%)

---

## Phase 6 — Guidance, Navigation & Control (GNC)

Closed-loop guidance for realistic trajectory shaping.

### Physics additions (`Physics.Mechanics.Propulsion`)
- [x] `GravityTurnGuidance` — pitch-over program: initial kick angle → gravity-assisted turn
- [x] `PEGGuidance` — Powered Explicit Guidance (iterative, used by Shuttle/Falcon): targets orbital elements
- [x] `AttitudeController` — PID or quaternion-feedback for maintaining commanded attitude via gimbal

### Engine additions (`Engines.Game.Rocket`)
- [x] `GuidanceComputer` — selects guidance law per phase, computes commanded attitude each step
- [x] `ThrustVectorControl` — gimbal deflection → moment about CG, limited by actuator rate/angle
- [x] `NavigationFilter` — simple state estimator (for game: perfect knowledge or configurable noise)

### ML integration (`ML.ReinforcementLearning.Environments`)
- [x] `RocketLandingEnv` — train RL agent for propulsive landing (SpaceX-style): observation = state, action = gimbal + throttle
- [x] `AscentOptimizationEnv` — optimize ascent trajectory for max payload to orbit

### Tests
- [x] Gravity turn reaches target orbit inclination within 0.1°
- [x] PEG converges to target orbit within 3 guidance cycles
- [x] Gimbal rate-limited response doesn't overshoot
- [x] RL landing agent achieves soft touchdown (v < 2 m/s) after training

---

## Phase 7 — Unity Integration & Visualization

Package as Unity-ready asset with real-time visualization hooks.

### Engine additions (`Engines.Game.Rocket`)
- [x] `RocketUnityAdapter` — converts RocketState to Unity Transform (position, rotation)
- [x] `TrajectoryPredictor` — fast forward propagation for trajectory line rendering
- [x] `TelemetryStream` — provides velocity, altitude, acceleration, fuel%, orbital elements per frame for HUD
- [x] `TimeWarp` — configurable simulation speed (1x–1000x) with interpolated rendering

### Documentation & Examples
- [x] Unity sample: Falcon 9-style two-stage to LEO with booster landing
- [x] Unity sample: Saturn V Moon trajectory (trans-lunar injection)
- [x] Telemetry HUD prefab (altitude, velocity, acceleration, apoapsis/periapsis)
- [x] Update Engines/README.md with Rocket Simulation section

### Tests
- [x] Adapter preserves state over 100,000 steps without drift
- [x] TimeWarp produces identical trajectory regardless of warp factor
- [x] TelemetryStream maintains 60fps data delivery at 100x warp

---

## Summary — Namespace Placement

| Feature | Namespace | Rationale |
|---------|-----------|-----------|
| RocketEngine, ThrustCurve, PropellantTank | `Physics.Mechanics.Propulsion` | Physical models |
| VariableMassDynamics, CG tracker | `Physics.Mechanics.Propulsion` | Physical law (F=ma+ṁv) |
| CompressibleDragModel, MachNumber | `Physics.FluidDynamics.Aerodynamics` | Aerodynamic physics |
| HeatFlux | `Physics.FluidDynamics.Aerodynamics` | Thermodynamic model |
| EarthModel, GravityModel, Coriolis | `Physics.Mechanics` | Physical models |
| OrbitalElements, Propagator, Hohmann | `Physics.OrbitalMechanics` | Physical/mathematical |
| GeoCoordinate, ECEF, ECI conversions | `Numerics.Objects` | Pure coordinate math |
| RocketSimulationEngine, stages, events | `Engines.Game.Rocket` | Orchestration |
| GuidanceComputer, TVC | `Engines.Game.Rocket` | Engine-level control |
| RocketLandingEnv | `ML.ReinforcementLearning.Environments` | Training environment |
| Unity adapter, telemetry | `Engines.Game.Rocket` | Platform bridge |

---

## Dependencies on Advanced Game Engine Roadmap

This roadmap builds on Phase 1 of the Advanced Game Engine Roadmap (atmosphere model, flight dynamics infrastructure). The following must be complete before starting:

- [x] `AtmosphereModel` (ISA standard atmosphere) — already implemented
- [x] `FlightDynamicsEngine` with RK4 integration — already implemented
- [x] `FrameTransforms` (body ↔ world) — already implemented
- [x] Drag/lift force extensions — already implemented
- [ ] `AirfoilModel` (Phase 1 of AdvancedGameEngine) — useful but not blocking for rockets (rockets use Cd tables, not airfoil polars)
