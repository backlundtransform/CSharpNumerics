# 🗺️ CSharpNumerics — Physics Roadmap

> **Context:** This project was built from a mathematics and physics background. Recent community interest — particularly from game developers — has highlighted the value of the kinematics module. This roadmap outlines a path from the current foundation toward a complete real-time physics toolkit, with **rigid body dynamics** as the next major milestone.

---

## ✅ What We Have Today

### Kinematics (`KinematicsExtensions.cs`)
- Free fall (scalar + vector)
- Constant velocity & acceleration (SUVAT, scalar + vector)
- Circular motion (centripetal acceleration, angular speed/velocity, period, frequency, tangential velocity)
- Projectile motion (position, velocity, time of flight, max height, range — scalar + vector)
- Orbital mechanics (gravitational field, force, orbital speed/period/position/velocity/acceleration, escape velocity)
- Relative motion (relative velocity/position/acceleration, closing speed, closest approach)

### ODE Solvers (`DifferentialEquationExtensions.cs`)
- Scalar RK4, Euler, Trapezoidal, custom Butcher tableau
- **Vector (3D)** RK4, Euler + trajectory
- **VectorN (N-dimensional)** RK4, Euler + trajectory — ready for 6D+ dynamics state

### Core Math Objects
- `Vector` (3D), `VectorN` (N-dimensional), `Matrix`, `Tensor`
- All with operator overloads (`+`, `-`, `*`, `/`, dot, cross, norm, etc.)

### Physics Constants
- Gravitation, electromagnetism, thermodynamics, quantum, astronomy

---

## 🎯 Phase 1 — Dynamics Foundation

**Goal:** Newton's laws applied to particles and rigid bodies. Enable force-driven simulation.

### 1.1 Particle Dynamics (`DynamicsExtensions.cs`) ✅ DONE
Forces acting on point masses — the bridge from kinematics to dynamics.

| Feature | Description | Builds On |
|---|---|---|
| `ApplyForce(mass, force) → Vector` | a = F/m (Newton's 2nd law) | Vector |
| `NetForce(params Vector[] forces) → Vector` | Sum of forces | Vector |
| `Momentum(mass, velocity) → Vector` | p = mv | Vector |
| `Impulse(force, dt) → Vector` | J = F·Δt | Vector |
| `KineticEnergy(mass, velocity) → double` | ½mv² | Vector.Dot |
| `PotentialEnergy(mass, height) → double` | mgh | PhysicsConstants |
| `GravitationalPotentialEnergy(m1, m2, r) → double` | -Gm₁m₂/r | PhysicsConstants |
| `Work(force, displacement) → double` | W = F·d | Vector.Dot |
| `Power(force, velocity) → double` | P = F·v | Vector.Dot |

### 1.2 Rigid Body Primitives ✅ DONE

| Feature | Description | Builds On |
|---|---|---|
| `struct RigidBody` | Mass, position, velocity, orientation, angular velocity, inertia tensor | Vector, Matrix |
| `MomentOfInertia` helpers | Sphere, box, cylinder, rod (solid + hollow) | double |
| `InertiaMatrix(shape) → Matrix` | 3×3 inertia tensor for standard shapes | Matrix |

### 1.3 Torque & Rotational Dynamics ✅ DONE

| Feature | Description | Builds On |
|---|---|---|
| `Torque(r, F) → Vector` | τ = r × F | Vector.Cross |
| `AngularMomentum(I, ω) → Vector` | L = Iω | Matrix * Vector |
| `AngularAcceleration(I, τ) → Vector` | α = I⁻¹τ | Matrix.Inverse |
| `RotationalKineticEnergy(I, ω) → double` | ½ωᵀIω | Matrix, Vector.Dot |

---

## 🎯 Phase 2 — Simulation Loop

**Goal:** Time-stepping a scene with multiple bodies. This is where game devs start to integrate.

### 2.1 State Integration

| Feature | Description | Builds On |
|---|---|---|
| `IntegrateState(RigidBody, forces, torques, dt)` | Update position, velocity, orientation, ω using RK4/Euler | VectorN ODE solvers |
| Semi-implicit Euler | `v += a·dt; x += v·dt` — stable for games | New solver method |
| Verlet integration | `x(t+dt) = 2x(t) - x(t-dt) + a·dt²` — energy-preserving | New solver method |

### 2.2 Orientation & Rotation

| Feature | Description | Builds On |
|---|---|---|
| `struct Quaternion` | (w, x, y, z), multiply, normalize, conjugate, slerp | New type |
| `Quaternion.ToMatrix() → Matrix` | 3×3 rotation matrix | Matrix |
| `Quaternion.Rotate(Vector) → Vector` | Apply rotation to a vector | Vector |
| `IntegrateOrientation(q, ω, dt) → Quaternion` | q += ½·dt·ω·q | Quaternion, Vector |

### 2.3 Common Force Models

| Feature | Description | Builds On |
|---|---|---|
| `SpringForce(k, restLength, p1, p2) → Vector` | Hooke's law: F = -k(|Δr| - L₀)·r̂ | Vector |
| `DampingForce(c, velocity) → Vector` | F = -cv | Vector |
| `DragForce(Cd, ρ, A, velocity) → Vector` | F = -½CdρA|v|v | Vector |
| `FrictionForce(μ, normal, velocity) → Vector` | Static/kinetic friction | Vector |

---

## 🎯 Phase 3 — Collision Detection & Response

**Goal:** Bodies interact. Critical for any game physics scenario.

### 3.1 Bounding Volumes

| Feature | Description |
|---|---|
| `struct AABB` | Axis-Aligned Bounding Box (min, max) |
| `struct BoundingSphere` | Center + radius |
| `Intersects(AABB, AABB) → bool` | Box-box overlap test |
| `Intersects(Sphere, Sphere) → bool` | Sphere-sphere overlap test |
| `Intersects(AABB, Sphere) → bool` | Box-sphere overlap test |

### 3.2 Narrow Phase

| Feature | Description |
|---|---|
| `struct ContactPoint` | Position, normal, penetration depth |
| `SphereSphereContact → ContactPoint?` | Exact sphere collision |
| `SphereAABBContact → ContactPoint?` | Sphere-box collision |
| `GJK / SAT` (stretch) | Convex-convex collision (advanced) |

### 3.3 Collision Response

| Feature | Description | Builds On |
|---|---|---|
| `ResolveCollision(body1, body2, contact)` | Impulse-based resolution | RigidBody, Vector.Cross |
| Coefficient of restitution (e) | Elastic/inelastic control | double |
| Friction impulse | Tangential impulse at contact | Vector |
| Positional correction | Prevent sinking (Baumgarte stabilization) | Vector |

---

## 🎯 Phase 4 — Constraints & Joints

**Goal:** Connected bodies — pendulums, chains, ragdolls, vehicles.

| Feature | Description |
|---|---|
| `DistanceConstraint` | Maintains fixed distance between two points |
| `HingeJoint` | Single-axis rotation (doors, elbows) |
| `BallSocketJoint` | Free rotation, fixed position (shoulders) |
| `SpringJoint` | Elastic connection between bodies |
| Iterative constraint solver | Sequential impulse method |

---

## 🎯 Phase 5 — Spatial Optimization & Scene

**Goal:** Handle many objects efficiently.

| Feature | Description |
|---|---|
| `BroadPhase` | Sweep-and-prune or spatial hashing for O(n log n) pair reduction |
| `PhysicsWorld` | Container: list of bodies, forces, constraints, step function |
| `PhysicsWorld.Step(dt)` | Integrate → broad phase → narrow phase → resolve → update |
| Fixed timestep accumulator | Deterministic, framerate-independent simulation |

---

## 💡 Bonus Modules (Community-Driven)

These are not on the critical path but could drive significant interest:

| Module | Description | Depends On |
|---|---|---|
| **Raycasting** | `Ray.Intersect(AABB)`, `Ray.Intersect(Sphere)`, `Ray.Intersect(Plane)` | Vector, new types |
| **Soft body basics** | Spring-mass systems using existing `SpringForce` + ODE solvers | Phase 2 |
| **Fluid particles (SPH)** | Smoothed-particle hydrodynamics using VectorN ODE + spatial hashing | Phase 5 |
| **Character controller** | Grounded check, slope limits, step height | Phase 3 |

---

## 🔗 Dependency Graph

```
Phase 1 ─── Dynamics Foundation
  │           (forces, RigidBody, torque)
  │
  ├── uses ── Existing: Vector, Matrix, VectorN, ODE solvers, PhysicsConstants
  │
  ▼
Phase 2 ─── Simulation Loop
  │           (Quaternion, integration, force models)
  │
  ▼
Phase 3 ─── Collision Detection & Response
  │           (AABB, spheres, contact points, impulses)
  │
  ▼
Phase 4 ─── Constraints & Joints
  │           (distance, hinge, ball-socket, spring)
  │
  ▼
Phase 5 ─── Scene & Optimization
              (broad phase, PhysicsWorld, fixed timestep)
```

---

## 📐 How Current Code Maps to Game Physics

| Game Physics Need | Already Exists | To Build |
|---|---|---|
| Move objects | ✅ Position, Velocity, Acceleration | — |
| Gravity | ✅ Free fall, Gravitational force | — |
| Projectiles | ✅ Full projectile motion suite | — |
| Orbits | ✅ Circular orbit full suite | Elliptical orbits (Kepler) |
| Spin & rotate | ✅ Angular velocity, Torque (cross product) | Quaternion, inertia tensor |
| Collisions | ❌ | Phase 3 |
| Springs & ropes | ❌ | Phase 2 force models |
| Ragdolls | ❌ | Phase 4 constraints |
| Scene management | ❌ | Phase 5 |
| Numerical integration | ✅ RK4, Euler, VectorN trajectories | Verlet, semi-implicit Euler |
| Energy conservation check | ✅ VectorN.Dot, trajectory energy audit | — |

---

## 🚀 Suggested First Implementation

Start **Phase 1.1** — it requires no new types, only new extension methods on existing `Vector` and `double`:

```csharp
// Example: what Phase 1.1 looks like
var force = new Vector(0, 0, -9.8 * mass);
var acceleration = force.ApplyForce(mass);    // a = F/m
double ke = mass.KineticEnergy(velocity);     // ½mv²
double work = force.Work(displacement);       // F·d

// Then immediately simulate with existing ODE solver:
Func<(double t, VectorN y), VectorN> dynamics = v =>
{
    var pos = new Vector(v.y[0], v.y[1], v.y[2]);
    var vel = new Vector(v.y[3], v.y[4], v.y[5]);
    var netForce = gravity + SpringForce(k, L0, pos, anchor) + DampingForce(c, vel);
    var acc = netForce / mass;
    return new VectorN([vel.x, vel.y, vel.z, acc.x, acc.y, acc.z]);
};
```

No new solver code needed — the `VectorN` RK4 already handles 6D+ state integration.
