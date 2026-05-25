## 🎮 Applied Physics

Simulation-specific code for game engines and real-time physics — separated from the fundamental physics in `Physics/`.

**Namespace:** `CSharpNumerics.Physics.Applied`

---

### Bounding Volumes

Fast broad-phase overlap tests for spatial queries.

**AABB (Axis-Aligned Bounding Box)**

```csharp
var box = new AABB(new Vector(0, 0, 0), new Vector(2, 2, 2));
var box2 = AABB.FromCenterExtents(new Vector(5, 5, 5), new Vector(1, 1, 1));

Vector center = box.Center;           // (1, 1, 1)
Vector half = box.HalfExtents;        // (1, 1, 1)
bool inside = box.Contains(new Vector(1, 1, 1));  // true
Vector closest = box.ClosestPoint(new Vector(5, 1, 1)); // (2, 1, 1)
AABB merged = box.Merge(box2);        // smallest box containing both
AABB padded = box.Expand(0.1);        // uniform margin
```

**BoundingSphere**

```csharp
var sphere = new BoundingSphere(new Vector(0, 0, 0), 5);

bool inside = sphere.Contains(new Vector(3, 0, 0));  // true
Vector closest = sphere.ClosestPoint(new Vector(10, 0, 0)); // (5, 0, 0)
BoundingSphere merged = sphere.Merge(otherSphere);
```

---

### Overlap Tests (Broad Phase)

```csharp
bool hit = boxA.Intersects(boxB);           // AABB vs AABB
bool hit2 = sphereA.Intersects(sphereB);    // Sphere vs Sphere
bool hit3 = box.Intersects(sphere);          // AABB vs Sphere
bool hit4 = sphere.Intersects(box);          // Sphere vs AABB
```

---

### Contact Generation (Narrow Phase)

Returns `ContactPoint?` — null if no overlap. Normal points from A toward B.

```csharp
// Sphere vs Sphere
ContactPoint? contact = sphereA.SphereSphereContact(sphereB);

// Sphere vs AABB
ContactPoint? contact2 = sphere.SphereAABBContact(box);

// AABB vs Sphere (reversed normal)
ContactPoint? contact3 = box.AABBSphereContact(sphere);

if (contact is ContactPoint c)
{
    Vector pos = c.Position;           // world-space contact point
    Vector normal = c.Normal;          // unit normal from A → B
    double depth = c.PenetrationDepth; // overlap distance
}
```

---

### Collision Response

Impulse-based resolution with angular effects, Coulomb friction, and Baumgarte positional correction.

```csharp
var a = RigidBody.CreateSolidSphere(mass: 2, radius: 1);
a.Position = new Vector(0, 0, 0);
a.Velocity = new Vector(5, 0, 0);

var b = RigidBody.CreateSolidSphere(mass: 3, radius: 1);
b.Position = new Vector(3, 0, 0);
b.Velocity = new Vector(-2, 0, 0);

// Detect
var sA = new BoundingSphere(a.Position, 1);
var sB = new BoundingSphere(b.Position, 1);
var contact = sA.SphereSphereContact(sB);

if (contact is ContactPoint c)
{
    // Resolve velocity (normal + friction impulse)
    CollisionResponse.ResolveCollision(ref a, ref b, c,
        restitution: 0.8,   // 0 = sticky, 1 = fully elastic
        friction: 0.3);     // Coulomb friction coefficient

    // Fix overlap (Baumgarte stabilization)
    CollisionResponse.CorrectPositions(ref a, ref b, c,
        correctionFraction: 0.4,  // how aggressive (0.2–0.8)
        slop: 0.01);              // ignore tiny penetrations
}
```

**Properties of the impulse solver:**
- Conserves linear momentum exactly
- Conserves kinetic energy when `restitution = 1.0`
- Handles static (immovable) bodies — only the dynamic body moves
- Includes angular velocity effects via `InverseInertiaTensorWorld`
- Friction impulse is clamped to the Coulomb cone (`|jt| ≤ μ·jn`)

---

### Constraints & Joints

Sequential impulse solver for connected bodies. Each constraint corrects velocities iteratively until convergence.

**Simulation loop with constraints:**

```csharp
// Setup
var bodies = new RigidBody[]
{
    RigidBody.CreateStatic(new Vector(0, 0, 10)),   // ceiling
    RigidBody.CreateSolidSphere(mass: 1, radius: 0.2),
    RigidBody.CreateSolidSphere(mass: 1, radius: 0.2),
};
bodies[1].Position = new Vector(0, 0, 7);
bodies[2].Position = new Vector(0, 0, 4);

var constraints = new IConstraint[]
{
    new DistanceConstraint(0, 1, new Vector(0,0,0), new Vector(0,0,0), distance: 3),
    new DistanceConstraint(1, 2, new Vector(0,0,0), new Vector(0,0,0), distance: 3),
};

// Each frame:
double dt = 0.001;
for (int i = 0; i < bodies.Length; i++)
    if (!bodies[i].IsStatic)
        bodies[i].Velocity += dt * new Vector(0, 0, -9.8);

ConstraintSolver.Solve(bodies, constraints, dt, iterations: 10);

for (int i = 0; i < bodies.Length; i++)
    if (!bodies[i].IsStatic)
        bodies[i].Position += dt * bodies[i].Velocity;
```

**Constraint types:**

```csharp
// Distance — rigid rod between two anchors
var rod = new DistanceConstraint(bodyA: 0, bodyB: 1,
    localAnchorA: new Vector(0, 0, 0),
    localAnchorB: new Vector(0, 0, 0),
    distance: 3.0);

// Ball-socket — shared pivot, free rotation (shoulder, ragdoll)
var socket = BallSocketJoint.FromWorldPivot(0, 1,
    worldPivot: new Vector(0, 0, 4), bodies);

// Hinge — single-axis rotation (door, elbow, wheel)
var hinge = HingeJoint.FromWorldPivot(0, 1,
    worldPivot: new Vector(0, 0, 5),
    hingeAxis: new Vector(0, 0, 1), bodies);

// Spring-damper — soft elastic connection
var spring = new SpringJoint(0, 1,
    new Vector(0, 0, 0), new Vector(0, 0, 0),
    stiffness: 20, damping: 2, restLength: 1);
```

**Solver properties:**
- Sequential impulse method — iterates for convergence (4–20 iterations typical)
- Baumgarte positional stabilization prevents drift
- Handles mixed static + dynamic bodies
- Hard constraints (Distance, BallSocket, Hinge) iterated; soft constraints (SpringJoint) applied once

---

### PhysicsWorld

The main orchestrator — manages bodies, constraints, collisions, and the full simulation pipeline in one call.

```csharp
var world = new PhysicsWorld
{
    Gravity = new Vector(0, 0, -9.8),
    DefaultRestitution = 0.7,
    DefaultFriction = 0.3,
    SolverIterations = 10,
    FixedTimeStep = 1.0 / 60.0,
};

// Add bodies (returns index for constraints/queries)
var floor = RigidBody.CreateStatic(new Vector(0, 0, 0));
int iFloor = world.AddBody(floor, boundingRadius: 100);

var ball = RigidBody.CreateSolidSphere(mass: 1, radius: 0.5);
ball.Position = new Vector(0, 0, 10);
int iBall = world.AddBody(ball, boundingRadius: 0.5);

// Add constraints
world.AddConstraint(new DistanceConstraint(0, 1,
    new Vector(0, 0, 0), new Vector(0, 0, 0), distance: 5));

// Collision callback
world.OnCollision = (a, b, contact) =>
    Console.WriteLine($"Collision: body {a} ↔ body {b}, depth={contact.PenetrationDepth:F3}");

// Modify bodies directly via ref
ref var b = ref world.Body(iBall);
b.Velocity = new Vector(3, 0, 0);
```

**Simulation step — full pipeline in one call:**

```csharp
// Option 1: Single fixed step
world.Step(dt: 0.01);
// Pipeline: gravity → constraints → positions → broadphase → collisions

// Option 2: Fixed-timestep accumulator (deterministic, framerate-independent)
int steps = world.Update(elapsed: deltaTime); // consumes time in fixed chunks

// Smooth rendering between physics steps
double alpha = world.Alpha; // interpolation factor [0, 1)
```

**Broad phase algorithms:**

```csharp
// Default: sweep-and-prune on X axis — O(n log n)
var world = new PhysicsWorld(); // uses SweepAndPruneBroadPhase

// Brute force — O(n²), simpler, for small scenes
var world2 = new PhysicsWorld(new BruteForceBroadPhase());

// BVH — top-down bounding volume hierarchy, good for large worlds
var world3 = new PhysicsWorld(new BVHBroadPhase());
```

---

### Continuous Collision Detection (CCD)

Swept-sphere tests prevent fast-moving objects from tunnelling through geometry.

```csharp
using CSharpNumerics.Engines.Game;

// Sphere moving from A to B — does it hit a static sphere?
var result = ContinuousCollisionDetection.SweptSphereVsSphere(
    startPos: new Vector(0, 0, 0),
    endPos: new Vector(10, 0, 0),
    radiusA: 0.5,
    center: new Vector(5, 0, 0),
    radiusB: 1.0);

if (result.Hit)
{
    double toi = result.TimeOfImpact;   // [0, 1]
    Vector hitPos = result.HitPosition;
    Vector hitNormal = result.HitNormal;
}

// Swept sphere vs AABB (Minkowski expansion + slab test)
var box = new AABB(new Vector(4, -1, -1), new Vector(6, 1, 1));
var r2 = ContinuousCollisionDetection.SweptSphereVsAABB(
    start, end, radius: 0.5, box);

// Swept sphere vs ground plane
var r3 = ContinuousCollisionDetection.SweptSphereVsPlane(
    start, end, radius: 0.5, planeHeight: 0);
```

---

## ✈️ Flight Dynamics

6-degree-of-freedom aircraft simulation with aerodynamic forces, atmosphere model, and RK4 integration.

**Namespace:** `CSharpNumerics.Engines.Game.Flight`

### Architecture

```
FlightDynamicsEngine
  ├── AircraftConfig (geometry, mass, control surfaces, propulsion)
  ├── AircraftState  (12-state vector: pos, vel, quaternion, angular rates)
  ├── ControlInput   (throttle, pitch, roll, yaw, flaps, gear)
  └── Physics layer
       ├── AtmosphereModel   (ISA standard atmosphere)
       ├── AirfoilModel      (Cl/Cd vs AoA: NACA symmetric, flat plate)
       ├── ControlSurface    (elevator, aileron, rudder ΔCl/ΔCd)
       └── PropulsionModel   (jet + propeller thrust models)
```

### Quick Start

```csharp
using CSharpNumerics.Engines.Game.Flight;

// Use a preset aircraft configuration
var config = AircraftConfig.GenericLightAircraft();
var engine = new FlightDynamicsEngine(config);
engine.Init();

// Set initial conditions: 1000 m altitude, 50 m/s airspeed heading north
engine.SetState(new AircraftState(
    position: new Vector(0, 0, -1000),  // NED: z = -altitude
    velocity: new Vector(50, 0, 0),
    attitude: Quaternion.Identity,
    angularRate: new Vector(0, 0, 0)));

// Fly with 60% throttle, slight nose-up
var input = new ControlInput(throttle: 0.6, pitch: 0.05, roll: 0, yaw: 0);
engine.SetInput(input);

// Simulate 10 seconds
for (int i = 0; i < 1000; i++)
    engine.Step(0.01);

double altitude = engine.State.Altitude;    // metres above ground
double airspeed = engine.State.Airspeed;    // m/s
var (roll, pitch, yaw) = engine.State.EulerAngles;
```

### Aircraft Presets

```csharp
var cessna = AircraftConfig.GenericLightAircraft();  // Cessna 172-like
var jet = AircraftConfig.GenericJet();               // Business jet
```

### Wind Coupling

Attach a fluid solver to inject wind into the flight model and feed exhaust back:

```csharp
var fluid = new GameFluidSolver3D(new FluidConfig { GridX = 64, GridY = 64, GridZ = 64 });
engine.FluidField = fluid;

// Aircraft now samples wind at its position
engine.Step(0.01);
Vector wind = engine.WindVelocity;  // sampled from fluid field
```

---

## 💨 Real-Time Fluid Simulation

Game-quality Navier-Stokes solver based on Jos Stam's Stable Fluids — unconditionally stable at any timestep.

**Namespace:** `CSharpNumerics.Engines.Game.Fluids`

### 2D Solver

```csharp
using CSharpNumerics.Engines.Game.Fluids;

var config = new FluidConfig
{
    GridX = 128, GridY = 128,
    Viscosity = 0.0001,
    VorticityConfinementStrength = 2.0,  // sharper smoke curls
    EnableBuoyancy = true,               // hot gas rises
};
config.Quality = FluidQuality.High;      // 40 Poisson iterations

var solver = new GameFluidSolver2D(config);

// Smoke emitter at centre
var emitter = new FluidEmitter(
    position: new Vector(64, 10, 0),
    densityRate: 10.0,
    radius: 3);
emitter.Velocity = new Vector(0, 5, 0);  // upward
emitter.Temperature = 600;               // hot → buoyant
solver.AddEmitter(emitter);

// Obstacle (building)
solver.AddObstacle(FluidObstacle.Box2D(50, 40, 70, 60));

// Step the simulation
solver.Step(dt: 0.016);

// Read density field for rendering
ReadOnlySpan<double> density = solver.Density;
(double u, double v) = solver.SampleVelocity(worldX: 5.0, worldY: 3.0);
```

### 3D Solver

```csharp
var config3D = new FluidConfig { GridX = 64, GridY = 64, GridZ = 64 };
var solver3D = new GameFluidSolver3D(config3D);

solver3D.AddEmitter(new FluidEmitter(new Vector(32, 32, 5), densityRate: 8));
solver3D.Step(0.016);

// Sample velocity at a world point
var (vx, vy, vz) = solver3D.SampleVelocity(worldX: 3.0, worldY: 3.0, worldZ: 1.0);
```

### Fluid-Body Coupling

Two-way interaction: bodies feel drag from the fluid; bodies displace fluid.

```csharp
using CSharpNumerics.Engines.Game.Fluids;

// Compute drag on a sphere in a fluid stream
Vector drag = FluidBodyCoupling.ComputeDragForce(
    bodyVelocity: new Vector(5, 0, 0),
    fluidVelocity: new Vector(-2, 0, 0),
    dragCoefficient: 0.47,     // sphere
    crossSectionArea: 0.785,   // π·r²
    fluidDensity: 1.225);      // air at sea level

// Body pushes fluid away
FluidBodyCoupling.DisplaceFluid2D(u, v, nx, ny,
    bodyGridX: 32, bodyGridY: 32, bodyRadius: 3,
    bodyVelX: 2, bodyVelY: 0);
```

---

## 🤖 ML-Driven Game AI

Wire existing RL infrastructure into the game engine for intelligent NPC behaviour.

**Namespace:** `CSharpNumerics.Engines.Game.AI`

### Training a Flight AI

```csharp
using CSharpNumerics.Engines.Game.AI;
using CSharpNumerics.ML.ReinforcementLearning.Environments;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;

// Create RL environment wrapping the flight engine
var env = new FlightEnv();

// Create a PPO agent
var agent = new PPOAgent(
    observationSize: env.ObservationSize,
    actionSize: env.ActionSize);

// Train offline
var trainer = new AITrainer(env, agent)
    .WithEpisodes(500, maxStepsPerEpisode: 1000)
    .WithSeed(42);

GameAIAgent trainedAgent = trainer.Train("FlightAI");
double avgReturn = trainer.Evaluate(numEpisodes: 10);
```

### Behavior Trees

Combine ML decisions with scripted fallback logic:

```csharp
var tree = new BehaviorTree("CombatAI",
    new SelectorNode("Root",
        new SequenceNode("Engage",
            new ConditionNode("InRange", ctx => ctx.Get<double>("distance") < 500),
            new ActionNode("MLAttack", ctx =>
            {
                // Use RL policy for attack maneuver
                var action = aiAgent.Act(ctx.Get<VectorN>("observation"));
                ctx.Set("action", action);
                return NodeStatus.Success;
            })),
        new SequenceNode("Evade",
            new ConditionNode("LowHealth", ctx => ctx.Get<double>("health") < 0.3),
            new ActionNode("RunAway", ctx => { /* flee logic */ return NodeStatus.Success; })),
        new ActionNode("Patrol", ctx => { /* patrol waypoints */ return NodeStatus.Running; })
    ));

var ctx = new BehaviorContext();
ctx.Set("distance", 300.0);
ctx.Set("health", 0.8);
tree.Tick(ctx);
```

### Formation Flying

```csharp
var leader = new FlightDynamicsEngine(AircraftConfig.GenericJet());
leader.Init();

var formation = new FormationController(leader);
var wingman1 = new FlightDynamicsEngine(AircraftConfig.GenericJet());
wingman1.Init();
formation.AddWingman("Wing 2", offset: new Vector(-30, 30, 0), engine: wingman1);

// Step all aircraft
formation.Step(dt: 0.01);
double[] errors = formation.GetPositionErrors();  // deviation from desired offset
```

### Adaptive Difficulty

```csharp
var difficulty = new AdaptiveDifficulty(initialDifficulty: 0.5, targetPerformance: 0.5);

// After each round, record player performance
difficulty.RecordPerformance(0.8);  // player doing well → difficulty increases
difficulty.RecordPerformance(0.2);  // player struggling → difficulty decreases

// Apply to AI parameters
var aiParams = difficulty.GetAIParameters();
double reactionDelay = aiParams["reactionDelay"];   // 0–0.5s
double accuracy = aiParams["accuracyScale"];         // 0.3–1.0
```

---

## 🧱 Soft Body & Cloth

Mass-spring deformable mesh and cloth simulation using Verlet integration.

**Namespace:** `CSharpNumerics.Physics.Mechanics.SoftBody`

### Deformable Mesh

```csharp
using CSharpNumerics.Physics.Mechanics.SoftBody;

// Create a 10×10 grid mesh, 1m spacing, each vertex = 0.1 kg
var mesh = DeformableMesh.CreateGrid(
    width: 10, height: 10, resX: 10, resY: 10,
    mass: 0.1, stiffness: 0.8);

// Pin two corners
mesh.Pin(0);
mesh.Pin(9);

// Add sphere obstacle at (5, 5, -3)
mesh.Step(dt: 0.01);
mesh.CollideWithSphere(center: new Vector(5, 5, -3), radius: 2);
mesh.CollideWithGround(groundZ: -10);
```

### Cloth Simulation

```csharp
var cloth = new ClothSimulation(
    width: 5, height: 5, resX: 20, resY: 20,
    mass: 0.05, stiffness: 0.9);
cloth.PinTopEdge();                          // hang from top
cloth.Wind = new Vector(2, 0, 0);           // apply wind
cloth.EnableSelfCollision = true;
cloth.SelfCollisionRadius = 0.02;

for (int i = 0; i < 500; i++)
    cloth.Step(0.01);

Vector pos = cloth.GetPosition(10, 15);     // query vertex position
```

---

## 💧 SPH Fluid & Free Surface

### Smoothed Particle Hydrodynamics

Lagrangian particle-based liquid simulation. Poly6/Spiky/Viscosity kernels.

**Namespace:** `CSharpNumerics.Physics.FluidDynamics.SPH`

```csharp
using CSharpNumerics.Physics.FluidDynamics.SPH;

// Create a block of water particles
var positions = new List<Vector>();
for (int x = 0; x < 10; x++)
    for (int y = 0; y < 10; y++)
        for (int z = 0; z < 10; z++)
            positions.Add(new Vector(x * 0.05, y * 0.05, z * 0.05 + 0.5));

var sph = new SPHSolver(positions)
{
    SmoothingRadius = 0.1,
    ParticleMass = 0.02,
    RestDensity = 1000,
    GasConstant = 2000,
    Viscosity = 1.0,
    Gravity = new Vector(0, 0, -9.81),
    BoundsMin = new Vector(-1, -1, 0),
    BoundsMax = new Vector(1, 1, 2),
    Dt = 0.001,
};

for (int i = 0; i < 100; i++)
    sph.Step();

ReadOnlySpan<SPHSolver.Particle> particles = sph.Particles;
```

### Volume-of-Fluid Tracker

Eulerian free-surface tracking on a grid. Each cell stores a liquid fraction F ∈ [0,1].

**Namespace:** `CSharpNumerics.Physics.FluidDynamics.FreeSurface`

```csharp
using CSharpNumerics.Physics.FluidDynamics.FreeSurface;

var vof = new VOFTracker(nx: 64, ny: 64, cellSize: 0.1);

// Fill a circular region with water
vof.FillCircle(cx: 3.2, cy: 3.2, radius: 1.0);

// Advect with a velocity field
vof.Advect(u, v, dt: 0.01);

// Surface normal at an interface cell
var (nx, ny) = vof.EstimateSurfaceNormal(i: 32, j: 32);
```

---

## 🌋 Particle System

Game particle system with emission, lifetime, physics forces, and ground collision.

**Namespace:** `CSharpNumerics.Engines.Game.Particles`

```csharp
using CSharpNumerics.Engines.Game.Particles;

var particles = new ParticleSystem(maxParticles: 10000, seed: 42)
{
    Gravity = new Vector(0, 0, -9.81),
    DragCoefficient = 0.1,
    GroundHeight = 0,
    GroundRestitution = 0.3,
};

// Cone emitter: sparks shooting upward
var emitter = new ParticleEmitter
{
    Position = new Vector(0, 0, 0),
    Direction = new Vector(0, 0, 1),
    EmissionRate = 200,
    ConeAngle = 0.4,
    InitialSpeed = 10,
    SpeedVariance = 0.3,
    ParticleLifetime = 2.0,
};
particles.AddEmitter(emitter);

// Simulate
particles.Update(dt: 0.016);

int alive = particles.AliveCount;
ReadOnlySpan<ParticleSystem.ParticleData> live = particles.Particles;

// Couple with fluid wind field
particles.ApplyWindField(
    solver.VelocityX, solver.VelocityY,
    gridNx: 66, gridNy: 66,
    gridOriginX: 0, gridOriginY: 0,
    cellSize: 1.0, strength: 5.0, dt: 0.016);
```

---

## ⛰️ Terrain Interaction

Heightmap-based terrain for collision and wind deflection.

**Namespace:** `CSharpNumerics.Engines.Game.Terrain`

### Terrain Collider

```csharp
using CSharpNumerics.Engines.Game.Terrain;

// Create terrain from a height function (hill)
var terrain = TerrainCollider.FromFunction(
    resX: 100, resY: 100, cellSize: 1.0,
    heightFn: (x, y) => 5.0 * Math.Exp(-(x * x + y * y) / 200));

// Sample height at a world point
double h = terrain.SampleHeight(wx: 10, wy: 15);
Vector normal = terrain.SurfaceNormal(10, 15);

// Sphere collision test
bool hit = terrain.TestSphere(position, radius: 0.5, out double pen, out Vector n);

// Resolve: push sphere above terrain, reflect velocity
terrain.ResolveSphere(ref position, ref velocity, radius: 0.5, restitution: 0.3);
```

### Wind Over Terrain

```csharp
// Modify a wind field based on terrain slope (Jackson-Hunt model)
WindOverTerrain.ApplyTerrainEffect(
    windX, windY, terrain,
    gridNx: 64, gridNy: 64,
    gridOriginX: 0, gridOriginY: 0,
    gridCellSize: 1.0,
    speedUpFactor: 1.5);   // ridge acceleration

// Zero out wind inside solid terrain
WindOverTerrain.ApplyTerrainBlocking(
    windX, windY, terrain,
    gridNx: 64, gridNy: 64,
    gridOriginX: 0, gridOriginY: 0,
    gridCellSize: 1.0,
    windHeight: 10);
```

---

## 🔌 Unity Integration

Thin adapter layer for using CSharpNumerics in Unity projects. No Unity dependency — defines surrogate types that mirror Unity's types.

**Namespace:** `CSharpNumerics.Engines.Game.Unity`

### Coordinate Conversion

CSharpNumerics uses right-handed Z-up; Unity uses left-handed Y-up.

```csharp
using CSharpNumerics.Engines.Game.Unity;
using static CSharpNumerics.Engines.Game.Unity.UnityAdapter;

// Position: CSN(x, y, z) → Unity(x, z, y)
UnityVector3 uPos = UnityAdapter.ToUnityVector3(new Vector(10, 20, 30));
Vector csnPos = UnityAdapter.FromUnityVector3(uPos);  // round-trip

// Rotation: 3×3 matrix ↔ Unity quaternion with axis swap
UnityQuaternion uQuat = UnityAdapter.ToUnityQuaternion(rotationMatrix);
Matrix mat = UnityAdapter.FromUnityQuaternion(uQuat);

// VectorN (abstract data, no coordinate swap)
UnityVector3 v3 = UnityAdapter.VectorNToUnityVector3(vectorN);
VectorN vn = UnityAdapter.UnityVector3ToVectorN(v3);
```

### Physics Synchronization

Interpolated physics state for smooth 60fps rendering at fixed physics rate:

```csharp
var world = new PhysicsWorld();
// ... add bodies ...
var sync = new PhysicsSync(world);

// In FixedUpdate:
sync.StepAndSync(dt: 0.02);   // 50 Hz physics

// In Update (rendering):
double alpha = 0.5; // interpolation factor
var pos = sync.GetInterpolatedPosition(bodyIndex: 0, alpha);
var rot = sync.GetInterpolatedRotation(bodyIndex: 0, alpha);

// Frame coherence check
bool ok = sync.CheckCoherence(maxDelta: 100);
```

### Fluid Rendering Bridge

Extract density/velocity textures for Unity VFX Graph or custom shaders:

```csharp
var renderer = new FluidRenderer(solver3D)
{
    DensityScale = 2.0,
    VelocityScale = 0.1,
    DensityThreshold = 0.01,
};

// Upload to Texture3D
float[] densityTex = renderer.UpdateDensityTexture();   // NX×NY×NZ floats [0,1]
float[] velocityTex = renderer.UpdateVelocityTexture(); // NX×NY×NZ×3 floats [-1,1]
float[] slice = renderer.GetDensitySlice(zIndex: 16);   // 2D slice
```

### Flight Controller (Unity MonoBehaviour Pattern)

```csharp
var controller = new FlightController(engine);

// In FixedUpdate — map Unity input axes
controller.SetInputAxis("throttle", 0.7);
controller.SetInputAxis("pitch", -0.3);
controller.SetInputAxis("roll", 0.1);
controller.SetGear(false);
controller.StepSimulation(dt: 0.02);

// In Update — read for rendering
UnityVector3 pos = controller.GetPosition();      // Y-up coordinates
UnityQuaternion rot = controller.GetRotation();
FlightController.HUDData hud = controller.GetHUDData();
// hud.Airspeed, hud.Altitude, hud.Heading, hud.Throttle, etc.
```

---

## 🚀 Rocket Simulation

Full 6DOF rocket launch simulation with multi-stage vehicles, orbital mechanics, GNC, and Unity integration.

**Namespace:** `CSharpNumerics.Engines.Game.Rocket`

### Architecture

```
RocketSimulationEngine (ISimulationEngine)
├── RocketVehicle (multi-stage stack)
│   ├── RocketStage[] (dry mass, engines, tanks, separation triggers)
│   └── Boosters (strap-on parallel stages)
├── GuidanceComputer (selects guidance law per phase)
│   ├── GravityTurnGuidance (atmospheric ascent)
│   ├── PEGGuidance (powered explicit guidance for orbit insertion)
│   └── AttitudeController (quaternion-feedback PID)
├── ThrustVectorControl (gimbal deflection → torque)
├── NavigationFilter (perfect or noisy state estimation)
├── MissionProfile (sequenced phases with exit conditions)
├── TelemetryRecorder (log states at configurable rate)
├── TelemetryStream (60fps HUD data delivery)
├── TrajectoryPredictor (fast Kepler propagation for orbit line)
├── TimeWarp (1x–1000x with fixed physics dt)
└── RocketUnityAdapter (NED/ECI → Unity Y-up coordinates)
```

### Quick Start

```csharp
// Build a two-stage rocket
var engine1 = new RocketEngine(thrustSL: 845000, thrustVac: 914000, ispSL: 282, ispVac: 311);
var tank1 = new PropellantTank(395700);
var stage1 = new RocketStage(dryMass: 25600, engines: new[] { engine1 }, tanks: new[] { tank1 });
stage1.SeparationTrigger = new StageSeparationTrigger(StageSeparationTriggerType.PropellantDepleted);

var engine2 = new RocketEngine(0, 934000, 0, 348);
var tank2 = new PropellantTank(92670);
var stage2 = new RocketStage(dryMass: 4000, engines: new[] { engine2 }, tanks: new[] { tank2 });

var vehicle = new RocketVehicle(new[] { stage1, stage2 }, payloadMass: 22800);
var sim = new RocketSimulationEngine(vehicle);

// Configure Earth model
sim.UseEarthModel = true;
sim.LaunchSite = GeoCoordinate.FromDegrees(28.5, -80.5, 0);
sim.SetStateFromLaunchSite(sim.LaunchSite.Value);
sim.Init();

// Run simulation
for (int i = 0; i < 10000; i++)
    sim.Step(0.1);

Console.WriteLine($"Altitude: {sim.State.Altitude:F0} m, Speed: {sim.State.Speed:F0} m/s");
```

### Guidance, Navigation & Control

```csharp
var gc = new GuidanceComputer();
gc.GravityTurn.KickAltitude = 500;
gc.GravityTurn.KickAngle = 3.0 * Math.PI / 180.0;
gc.GravityTurn.TargetInclination = 28.5 * Math.PI / 180.0;
gc.PEG.TargetSemiMajorAxis = 6378137 + 200000;
gc.PEG.Mu = 3.986004418e14;
gc.PEGActivationAltitude = 80000;
gc.Mode = GuidanceMode.Auto;

// Each simulation step:
gc.Update(state, time, altitude, thrustAccel, exhaustVelocity, dt);
Quaternion commanded = gc.CommandedAttitude;
Vector torque = gc.TorqueCommand;
```

### Thrust Vector Control

```csharp
var tvc = new ThrustVectorControl
{
    MaxGimbalAngle = 5.0 * Math.PI / 180.0,
    MaxGimbalRate = 10.0 * Math.PI / 180.0,
    MomentArm = 20.0
};

// From desired torque:
tvc.CommandFromTorque(desiredTorque, thrustMagnitude, dt);
Vector gimbalTorque = tvc.ComputeTorque(thrustMagnitude);
```

### Time Warp

```csharp
var warp = new TimeWarp
{
    PhysicsTimestep = 0.01,     // 100 Hz physics
    MaxWarpFactor = 1000
};

warp.WarpFactor = 10;           // 10x speed
warp.IncreaseWarp();            // jumps to next level (1→2→5→10→50→100→500→1000)

// Each render frame:
warp.StepEngine(sim, renderDeltaTime);
```

### Telemetry Stream (HUD)

```csharp
var stream = new TelemetryStream { TargetFrameRate = 60 };

// Each physics step:
stream.Push(state, simTime, thrustMagnitude, fuelFraction);

// Each render frame:
if (stream.TryDeliver(simTime))
{
    double alt = stream.Current.Altitude;
    double speed = stream.Current.Speed;
    double apo = stream.Current.ApoapsisAltitude;
    double peri = stream.Current.PeriapsisAltitude;
    double fuel = stream.Current.FuelPercent;
}
```

### Trajectory Prediction

```csharp
var predictor = new TrajectoryPredictor
{
    NumPoints = 200,
    PredictionHorizon = 5400,  // 90 minutes
    IncludeDrag = true
};

var points = predictor.Predict(position, velocity);
var (apo, peri) = predictor.PredictApsides(position, velocity);

// Render orbit line from points[i].Position
```

### Unity Adapter

```csharp
var adapter = new RocketUnityAdapter
{
    CoordinateMode = AdapterCoordinateMode.NEDToUnity,
    PositionScale = 1.0
};

adapter.Update(rocketState);
// Apply to Unity Transform:
// transform.position = ToUnityVector3(adapter.UnityPosition);
// transform.rotation = ToUnityQuaternion(adapter.UnityRotation);
```

### RL Environments

Two reinforcement learning environments for autonomous rocket control:

```csharp
// Propulsive landing (SpaceX-style)
var landingEnv = new RocketLandingEnv(dt: 0.1, maxSteps: 500);
var (obs, info) = landingEnv.Reset(seed: 42);
// obs: [altitude, vx, vy, vz, pitch, pitchRate, fuelFraction, speed]
// action: [throttle(0-1), gimbal(-1 to 1)]

// Ascent trajectory optimization
var ascentEnv = new AscentOptimizationEnv(dt: 1.0, maxSteps: 600);
var (obs2, info2) = ascentEnv.Reset();
// obs: [altitude, speed, gamma, downrange, fuel, Q, time, orbitalEnergy]
// action: [pitchRate(-1 to 1)]
```
```

### AI Bridge

Feed Unity game state to ML agents:

```csharp
var bridge = new AIBridge(trainedAgent);

// Each frame: set observation from Unity world state
bridge.SetObservationFromUnity(
    agentPosition: new UnityVector3(10, 5, 20),
    agentVelocity: new UnityVector3(1, 0, 2),
    targetPosition: new UnityVector3(50, 5, 30),
    extraFeatures: new double[] { health, ammo });

// Get AI action
VectorN action = bridge.GetAction();
UnityVector3 moveCmd = bridge.GetActionAsUnityVector3();  // with coord conversion
```

---

## Unity Examples

### Example 1 — Flight Simulator with HUD

Complete setup for a flight simulator scene: aircraft physics, input mapping, HUD overlay.

```csharp
// === FlightSimManager.cs (attach to empty GameObject) ===
using CSharpNumerics.Engines.Game.Flight;
using CSharpNumerics.Engines.Game.Unity;

// Setup
var config = AircraftConfig.GenericLightAircraft();
var engine = new FlightDynamicsEngine(config);
engine.Init();
engine.SetState(new AircraftState(
    new Vector(0, 0, -1000), new Vector(50, 0, 0),
    Quaternion.Identity, new Vector(0, 0, 0)));

var controller = new FlightController(engine);

// FixedUpdate (50 Hz physics)
controller.SetInputAxis("throttle", Input.GetAxis("Throttle"));
controller.SetInputAxis("pitch",    Input.GetAxis("Vertical"));
controller.SetInputAxis("roll",     Input.GetAxis("Horizontal"));
controller.SetInputAxis("yaw",      Input.GetAxis("Yaw"));
controller.StepSimulation(Time.fixedDeltaTime);

// Update (render at display rate)
var pos = controller.GetPosition();
var rot = controller.GetRotation();
transform.position = new Vector3(pos.x, pos.y, pos.z);
transform.rotation = new Quaternion(rot.x, rot.y, rot.z, rot.w);

// HUD overlay
var hud = controller.GetHUDData();
speedText.text  = $"{hud.Airspeed:F0} m/s";
altText.text    = $"{hud.Altitude:F0} m";
headingText.text = $"{hud.Heading:F0}°";
throttleBar.value = (float)hud.Throttle;
```

### Example 2 — Real-Time Smoke with Fluid-Body Interaction

Volumetric smoke driven by the 3D fluid solver, interacting with a moving sphere.

```csharp
// === SmokeSimManager.cs ===
using CSharpNumerics.Engines.Game.Fluids;
using CSharpNumerics.Engines.Game.Unity;

// Setup: 64³ fluid grid with buoyancy + vorticity confinement
var config = new FluidConfig
{
    GridX = 64, GridY = 64, GridZ = 64,
    Viscosity = 0.0001,
    VorticityConfinementStrength = 3.0,
    EnableBuoyancy = true,
};
var solver = new GameFluidSolver3D(config);
var renderer = new FluidRenderer(solver) { DensityScale = 5.0 };

// Ground-level smoke emitter
solver.AddEmitter(new FluidEmitter(new Vector(32, 32, 3), densityRate: 15)
{
    Velocity = new Vector(0, 0, 3),
    Temperature = 700,
});

// Obstacle: sphere at (32, 32, 20)
solver.AddObstacle(FluidObstacle.Box3D(28, 28, 16, 36, 36, 24));

// FixedUpdate
solver.Step(Time.fixedDeltaTime);

// Two-way coupling: sphere feels drag from fluid
var (fu, fv, fw) = solver.SampleVelocity(spherePos.x, spherePos.y, spherePos.z);
Vector fluidVel = new Vector(fu, fv, fw);
Vector drag = FluidBodyCoupling.ComputeDragForce(
    bodyVelocity, fluidVel, 0.47, Mathf.PI * radius * radius, 1.225);
sphereRb.AddForce(new Vector3((float)drag.x, (float)drag.z, (float)drag.y));

// Body displaces fluid
FluidBodyCoupling.DisplaceFluid3D(/* u, v, w arrays, body position/velocity */);

// Update VFX: upload density texture to Texture3D
float[] densityData = renderer.UpdateDensityTexture();
densityTexture3D.SetPixelData(densityData, 0);
densityTexture3D.Apply();
```

### Example 3 — ML-Trained Dogfight AI Opponents

Train AI pilots offline, then deploy as NPC opponents in a multiplayer dogfight.

```csharp
// === Training (offline, headless) ===
using CSharpNumerics.ML.ReinforcementLearning.Environments;
using CSharpNumerics.ML.ReinforcementLearning.Algorithms.PolicyGradient;
using CSharpNumerics.Engines.Game.AI;

var env = new DogfightEnv();
var ppo = new PPOAgent(env.ObservationSize, env.ActionSize);
var trainer = new AITrainer(env, ppo)
    .WithEpisodes(1000, maxStepsPerEpisode: 500)
    .WithSeed(42);
GameAIAgent redPilot = trainer.Train("RedLeader");

// === Runtime (in Unity scene) ===
using CSharpNumerics.Engines.Game.Unity;

var bridge = new AIBridge(redPilot);
var difficulty = new AdaptiveDifficulty(0.5, targetPerformance: 0.5);

// Each frame:
bridge.SetObservationFromUnity(
    agentPosition: aiAircraft.transform.ToUnityVector3(),
    agentVelocity: aiVelocity,
    targetPosition: playerAircraft.transform.ToUnityVector3(),
    extraFeatures: new double[] { aiHealth, aiAmmo, distToTarget });

VectorN action = bridge.GetAction();
aiController.SetInputAxis("throttle", action[0]);
aiController.SetInputAxis("pitch",    action[1]);
aiController.SetInputAxis("roll",     action[2]);
aiController.SetInputAxis("yaw",      action[3]);
aiController.StepSimulation(Time.fixedDeltaTime);

// After each round: adapt difficulty
difficulty.RecordPerformance(playerScore);
var aiParams = difficulty.GetAIParameters();
// Apply: add aiParams["explorationNoise"] to AI actions for sloppier flying
```

---

## API Reference

### Core Engine Types

| Type | Namespace | Description |
|------|-----------|-------------|
| `PhysicsWorld` | `Engines.Game` | Main simulation orchestrator — bodies, constraints, collisions |
| `CollisionDetection` | `Engines.Game` | AABB/sphere overlap tests and contact generation |
| `CollisionResponse` | `Engines.Game` | Impulse-based collision resolution with friction |
| `ContinuousCollisionDetection` | `Engines.Game` | Swept sphere CCD (sphere, AABB, plane) |
| `ConstraintSolver` | `Engines.Game` | Sequential impulse iterative solver |
| `RaycastExtensions` | `Engines.Game` | Ray vs AABB/sphere intersection |

### Bounding Volumes & Objects

| Type | Namespace | Description |
|------|-----------|-------------|
| `AABB` | `Engines.Game.Objects` | Axis-aligned bounding box |
| `BoundingSphere` | `Engines.Game.Objects` | Bounding sphere |
| `ContactPoint` | `Engines.Game.Objects` | Collision contact data |

### Constraints

| Type | Namespace | Description |
|------|-----------|-------------|
| `IConstraint` | `Engines.Game.Constraints` | Constraint interface |
| `DistanceConstraint` | `Engines.Game.Constraints` | Fixed-distance rod |
| `BallSocketJoint` | `Engines.Game.Constraints` | Free-rotation pivot |
| `HingeJoint` | `Engines.Game.Constraints` | Single-axis revolute joint |
| `SpringJoint` | `Engines.Game.Constraints` | Soft spring-damper |

### Broad Phase

| Type | Namespace | Description |
|------|-----------|-------------|
| `IBroadPhase` | `Engines.Game.BroadPhase` | Broad phase interface |
| `SweepAndPruneBroadPhase` | `Engines.Game.BroadPhase` | O(n log n) sweep on X axis |
| `BruteForceBroadPhase` | `Engines.Game.BroadPhase` | O(n²) for small scenes |
| `BVHBroadPhase` | `Engines.Game.BroadPhase` | Bounding volume hierarchy |

### Flight Dynamics

| Type | Namespace | Description |
|------|-----------|-------------|
| `FlightDynamicsEngine` | `Engines.Game.Flight` | 6DOF aircraft with RK4 integration |
| `AircraftState` | `Engines.Game.Flight` | 12-state vector (pos, vel, quat, ω) |
| `AircraftConfig` | `Engines.Game.Flight` | Aircraft geometry, mass, aero models |
| `ControlInput` | `Engines.Game.Flight` | Pilot input (throttle, pitch, roll, yaw) |

### Fluid Simulation

| Type | Namespace | Description |
|------|-----------|-------------|
| `GameFluidSolver2D` | `Engines.Game.Fluids` | 2D Stable Fluids solver |
| `GameFluidSolver3D` | `Engines.Game.Fluids` | 3D Stable Fluids solver |
| `FluidConfig` | `Engines.Game.Fluids` | Grid size, viscosity, quality settings |
| `FluidEmitter` | `Engines.Game.Fluids` | Inject density/velocity (smoke, exhaust) |
| `FluidObstacle` | `Engines.Game.Fluids` | Block flow (buildings, terrain) |
| `FluidBodyCoupling` | `Engines.Game.Fluids` | Two-way fluid-body interaction |
| `VorticityConfinement` | `Engines.Game.Fluids` | Re-inject lost vorticity for sharp curls |

### Game AI

| Type | Namespace | Description |
|------|-----------|-------------|
| `GameAIAgent` | `Engines.Game.AI` | Wraps trained RL policy for runtime use |
| `AITrainer` | `Engines.Game.AI` | Offline training loop |
| `BehaviorTree` | `Engines.Game.AI` | Behavior tree executor (Selector, Sequence, Action, Condition) |
| `FormationController` | `Engines.Game.AI` | Multi-agent formation flying |
| `AdaptiveDifficulty` | `Engines.Game.AI` | Dynamic difficulty adjustment |

### Particles

| Type | Namespace | Description |
|------|-----------|-------------|
| `ParticleSystem` | `Engines.Game.Particles` | Emission, lifetime, forces, ground collision |
| `ParticleEmitter` | `Engines.Game.Particles` | Cone emission, speed/lifetime variance |

### Terrain

| Type | Namespace | Description |
|------|-----------|-------------|
| `TerrainCollider` | `Engines.Game.Terrain` | Heightmap collision, surface normals |
| `WindOverTerrain` | `Engines.Game.Terrain` | Terrain-deflected wind (Jackson-Hunt) |

### Unity Integration

| Type | Namespace | Description |
|------|-----------|-------------|
| `UnityAdapter` | `Engines.Game.Unity` | Surrogate types + Z-up ↔ Y-up conversion |
| `PhysicsSync` | `Engines.Game.Unity` | Interpolated physics snapshots for rendering |
| `FluidRenderer` | `Engines.Game.Unity` | Density/velocity textures for VFX Graph |
| `FlightController` | `Engines.Game.Unity` | Input axis mapping + HUD data |
| `AIBridge` | `Engines.Game.Unity` | Unity game state → AI observation/action |

### Physics (used by Game Engine)

| Type | Namespace | Description |
|------|-----------|-------------|
| `DeformableMesh` | `Physics.Mechanics.SoftBody` | Mass-spring Verlet mesh |
| `ClothSimulation` | `Physics.Mechanics.SoftBody` | Cloth with wind + self-collision |
| `SPHSolver` | `Physics.FluidDynamics.SPH` | Smoothed Particle Hydrodynamics |
| `VOFTracker` | `Physics.FluidDynamics.FreeSurface` | Volume-of-Fluid surface tracking |
| `AtmosphereModel` | `Physics.FluidDynamics.Aerodynamics` | ISA standard atmosphere |
| `AirfoilModel` | `Physics.FluidDynamics.Aerodynamics` | Lift/drag vs angle of attack |
| `ControlSurface` | `Physics.FluidDynamics.Aerodynamics` | Elevator/aileron/rudder ΔCl/ΔCd |
| `PropulsionModel` | `Physics.FluidDynamics.Aerodynamics` | Jet + propeller thrust |
| `TurbulenceModel` | `Physics.FluidDynamics.Turbulence` | Smagorinsky SGS turbulence |
| `BuoyancyForce` | `Physics.FluidDynamics.Buoyancy` | Thermal buoyancy |
