# Advanced Game Engine Roadmap

## Goal

Extend the existing `Engines/Game` engine into a high-fidelity simulation platform targeting:
- **Flight simulators** (6DOF aircraft, aerodynamics, atmosphere)
- **Real-time fluid simulations** (wind, smoke, water — game-quality CFD)
- **ML-driven game AI** (NPC behavior, adaptive difficulty, physics-informed control)

The library will be distributed as Unity-compatible assets for developers building advanced simulation games.

## Architecture Principle

The engine **does not** replace or duplicate code in Numerics, Physics, or ML. It orchestrates them:

```
Engines.Game.Flight      → Physics.FluidDynamics + Physics.Mechanics + Numerics.ODE
Engines.Game.Fluids      → Physics.FluidDynamics + Numerics.FiniteDifference
Engines.Game.AI          → ML.ReinforcementLearning + ML.NeuralNetwork
Engines.Game (existing)  → Physics.Mechanics (RigidBody, collisions, constraints)
```

Namespace: `CSharpNumerics.Engines.Game.*`

---

## Phase 1 — Flight Dynamics Core

Flight simulation foundation: 6DOF rigid-body aircraft model with aerodynamic forces.

### Numerics additions (Numerics namespace)
- [ ] Quaternion-based rotation utilities (if not already sufficient — verify QuaternionAlgebra coverage for body↔world frame transforms)
- [ ] Frame transformation helpers (body frame ↔ world frame, wind frame ↔ body frame)

### Physics additions (Physics namespace)
- [ ] `Physics.FluidDynamics.Aerodynamics.AtmosphereModel` — ISA standard atmosphere (density, pressure, temperature vs altitude)
- [ ] `Physics.FluidDynamics.Aerodynamics.AirfoilModel` — lift/drag coefficient tables (AoA-based lookup with interpolation, flat plate, NACA symmetric)
- [ ] `Physics.FluidDynamics.Aerodynamics.ControlSurface` — deflection → ΔCl/ΔCd model (elevator, aileron, rudder)
- [ ] `Physics.FluidDynamics.Aerodynamics.PropulsionModel` — thrust as function of throttle, altitude, speed (jet + propeller)

### Engine additions (Engines.Game namespace)
- [ ] `Engines.Game.Flight.AircraftState` — 12-state vector (position, velocity, quaternion attitude, angular rates)
- [ ] `Engines.Game.Flight.FlightDynamicsEngine` — implements `ISimulationEngine`; integrates 6DOF equations using RK4Stepper from Numerics
- [ ] `Engines.Game.Flight.ControlInput` — throttle, pitch, roll, yaw, flaps, gear
- [ ] `Engines.Game.Flight.AircraftConfig` — wing area, span, mass, CG position, engine count, control surface geometry

### Tests
- [ ] Straight-and-level flight (trim condition)
- [ ] Stall behavior (AoA > critical → Cl drops)
- [ ] Climb/descent rates match expected performance
- [ ] Roll/pitch response to control inputs

---

## Phase 2 — Real-Time Fluid Simulation for Games

Adapt existing Navier-Stokes solvers for game-quality real-time performance. Focus on visual fidelity at interactive framerates rather than engineering accuracy.

### Physics additions (Physics namespace)
- [ ] `Physics.FluidDynamics.Turbulence.TurbulenceModel` — simple k-epsilon or Smagorinsky SGS for visual turbulence
- [ ] `Physics.FluidDynamics.Buoyancy.BuoyancyForce` — hot gas rises, cold sinks (for smoke/fire)

### Engine additions (Engines.Game namespace)
- [ ] `Engines.Game.Fluids.GameFluidSolver2D` — stripped-down 2D N-S optimized for speed (fewer Poisson iterations, coarser grid, LOD)
- [ ] `Engines.Game.Fluids.GameFluidSolver3D` — 3D variant with configurable resolution (32³–128³)
- [ ] `Engines.Game.Fluids.FluidConfig` — grid size, viscosity, timestep, boundary mode, quality preset (Low/Medium/High)
- [ ] `Engines.Game.Fluids.FluidBodyCoupling` — two-way: rigid bodies experience drag/lift from fluid; bodies displace fluid
- [ ] `Engines.Game.Fluids.VorticityConfinement` — add back vorticity lost to numerical diffusion (sharper smoke curls)
- [ ] `Engines.Game.Fluids.FluidEmitter` — inject velocity/density at a point or region (exhaust, wind sources, explosions)
- [ ] `Engines.Game.Fluids.FluidObstacle` — static/dynamic geometry that blocks flow (buildings, terrain, aircraft)

### Integration with Flight
- [ ] Aircraft samples wind field at its position → adds to aerodynamic velocity
- [ ] Engine exhaust feeds back into fluid as emitter
- [ ] Wake turbulence behind aircraft (simplified vortex model)

### Tests
- [ ] Smoke rising around obstacle (visual validation: vortex shedding)
- [ ] Fluid-body coupling: sphere falls slower in fluid vs vacuum
- [ ] Emitter produces expanding plume
- [ ] Performance: 64³ grid runs at ≥30 steps/sec on single thread

---

## Phase 3 — ML-Driven Game AI

Wire existing RL infrastructure into the game engine for intelligent NPC behavior and adaptive systems.

### ML additions (ML namespace)
- [ ] `ML.ReinforcementLearning.Environments.FlightEnv` — observation: aircraft state (12D), action: control surfaces (4D continuous), reward: waypoint tracking + fuel efficiency
- [ ] `ML.ReinforcementLearning.Environments.DogfightEnv` — multi-agent pursuit-evasion with flight dynamics
- [ ] `ML.ReinforcementLearning.Environments.FluidNavigationEnv` — agent navigates through wind/current field

### Engine additions (Engines.Game namespace)
- [ ] `Engines.Game.AI.GameAIAgent` — wraps a trained RL policy; takes game state observations, returns actions per tick
- [ ] `Engines.Game.AI.AITrainer` — offline training loop: runs FlightDynamicsEngine headless, trains PPO/DDPG policy
- [ ] `Engines.Game.AI.BehaviorTree` — lightweight behavior tree executor for hybrid AI (ML decisions + scripted fallbacks)
- [ ] `Engines.Game.AI.FormationController` — multi-agent coordination (wingmen following leader using RL + formation rules)
- [ ] `Engines.Game.AI.AdaptiveDifficulty` — monitors player performance, adjusts AI aggressiveness via policy parameters

### Tests
- [ ] Trained agent can maintain level flight (reward converges)
- [ ] Dogfight agent pursues target with reasonable intercept geometry
- [ ] Formation controller maintains spacing under wind perturbation
- [ ] Adaptive difficulty reduces AI skill when player struggles

---

## Phase 4 — Advanced Physics Integration

Deeper physics features needed for high-fidelity simulation games.

### Physics additions (Physics namespace)
- [ ] `Physics.Mechanics.SoftBody.DeformableMesh` — mass-spring network on triangulated mesh (uses CoupledOscillators)
- [ ] `Physics.Mechanics.SoftBody.ClothSimulation` — constrained particle system with self-collision
- [ ] `Physics.FluidDynamics.SPH.SPHSolver` — Smoothed Particle Hydrodynamics for water splashes, liquid in containers
- [ ] `Physics.FluidDynamics.FreeSurface.VOFTracker` — Volume-of-Fluid for water surface tracking

### Engine additions (Engines.Game namespace)
- [ ] `Engines.Game.Particles.ParticleSystem` — emission, lifetime, forces (gravity, drag, wind from fluid field), collision with world
- [ ] `Engines.Game.Particles.ParticleEmitter` — rate, cone angle, initial velocity, randomness
- [ ] `Engines.Game.Terrain.TerrainCollider` — heightmap-based collision for ground interaction
- [ ] `Engines.Game.Terrain.WindOverTerrain` — terrain deflects wind field (simple boundary layer model)
- [ ] Continuous collision detection (CCD) — swept sphere for fast projectiles
- [ ] Spatial partitioning upgrade — BVH or octree for large worlds

### Tests
- [ ] Cloth draped over sphere (visual: no penetration, natural drape)
- [ ] SPH water splashes when object dropped
- [ ] Particle system affected by fluid wind field
- [ ] Terrain collision for aircraft landing gear

---

## Phase 5 — Unity Integration Layer

Package the engine for consumption as a Unity asset. This layer is a thin adapter — all logic stays in the core library.

### Engine additions (Engines.Game namespace)
- [ ] `Engines.Game.Unity.UnityAdapter` — converts between Unity Vector3/Quaternion and CSharpNumerics VectorN/Matrix
- [ ] `Engines.Game.Unity.PhysicsSync` — synchronizes PhysicsWorld state with Unity Transform components (position, rotation)
- [ ] `Engines.Game.Unity.FluidRenderer` — provides density/velocity textures from GameFluidSolver3D for Unity VFX Graph or shader consumption
- [ ] `Engines.Game.Unity.FlightController` — MonoBehaviour-style interface: exposes ControlInput, reads AircraftState, drives Unity Transform
- [ ] `Engines.Game.Unity.AIBridge` — feeds Unity game state to GameAIAgent, applies returned actions

### Documentation & Examples
- [ ] Unity sample project: flight simulator scene with HUD
- [ ] Unity sample: real-time smoke simulation with fluid-body interaction
- [ ] Unity sample: ML-trained dogfight AI opponents
- [ ] API reference documentation for all public types
- [ ] Update Engines/README.md with Game Engine section

### Tests
- [ ] Adapter round-trips preserve values (Vector3 → VectorN → Vector3)
- [ ] PhysicsSync maintains frame coherence over 10,000 steps
- [ ] FluidRenderer produces valid texture data at 60fps

---

## Phase 6 — Performance & Polish

Optimize for production game use.

- [ ] SIMD acceleration for vector math hot paths (Vector128/256 where available)
- [ ] Span<T> / stackalloc for zero-allocation fluid solver inner loops
- [ ] Multi-threaded broad phase collision (partition spatial grid across threads)
- [ ] Fluid solver thread: run on background thread, game samples latest snapshot
- [ ] Memory pooling for particles, rigid bodies, constraint arrays
- [ ] Profiling benchmarks: measure ms/frame for key scenarios
- [ ] LOD system for fluid: coarse grid far from camera, fine near
- [ ] Deterministic replay: serialize full state for network sync / replays

---

## Summary — What Goes Where

| Feature | Namespace | Rationale |
|---------|-----------|-----------|
| Atmosphere model | `Physics.FluidDynamics.Aerodynamics` | Physical model, no game logic |
| Airfoil lift/drag | `Physics.FluidDynamics.Aerodynamics` | Physical model |
| Turbulence model | `Physics.FluidDynamics.Turbulence` | Physical model |
| SPH solver | `Physics.FluidDynamics.SPH` | Physical algorithm |
| Soft body / cloth | `Physics.Mechanics.SoftBody` | Physical model |
| Quaternion math | `Numerics.Objects` | Pure math |
| Flight dynamics engine | `Engines.Game.Flight` | Orchestration |
| Game fluid solver | `Engines.Game.Fluids` | Game-optimized orchestration |
| Game AI agent | `Engines.Game.AI` | Game-specific ML application |
| Particle system | `Engines.Game.Particles` | Game feature |
| Terrain interaction | `Engines.Game.Terrain` | Game feature |
| Unity adapter | `Engines.Game.Unity` | Platform integration |
| RL flight environment | `ML.ReinforcementLearning.Environments` | Training environment (ML domain) |

---

## Decision: New Engine or Extend Existing?

**Extend the existing `Engines/Game` engine.** Rationale:

1. The current PhysicsWorld + constraint solver + collision pipeline is solid and tested
2. Flight, fluids, and AI are **additive sub-systems** — they don't replace rigid-body physics
3. The `ISimulationEngine` + `SimulationClock` + `EventBus` infrastructure already supports multiple coordinated engines
4. Creating a separate engine would duplicate the rigid-body foundation and violate DRY
5. Unity integration is just an adapter layer — it doesn't require engine-level changes

The extended structure:
```
Engines/Game/
├── PhysicsWorld.cs          (existing — rigid body orchestrator)
├── Flight/                  (NEW — 6DOF flight dynamics)
├── Fluids/                  (NEW — real-time game fluid sim)
├── AI/                      (NEW — ML-driven game intelligence)
├── Particles/               (NEW — particle effects)
├── Terrain/                 (NEW — terrain interaction)
├── Unity/                   (NEW — Unity asset bridge)
├── Constraints/             (existing)
├── BroadPhase/              (existing)
└── Objects/                 (existing)
```
