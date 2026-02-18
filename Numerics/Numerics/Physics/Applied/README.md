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
