## 🛞 Kinematics Extensions

The `KinematicsExtensions` class provides a set of extension methods for performing common kinematic calculations in both **scalar** and **vector** form. It covers free fall, constant velocity, constant acceleration, time-independent SUVAT equations, and circular motion.

**Free Fall**

Compute the velocity or time of a freely falling object:

```csharp
double v = 10.0.FreeFallVelocity();       // scalar: v = sqrt(2*g*h)
double t = 10.0.FreeFallTime();           // scalar: t = sqrt(2*h/g)

// Vector version with optional direction (default downward)
Vector dir = new Vector(0,0,-1);
Vector vVec = 10.0.FreeFallVelocity(dir);
```

**Constant Velocity**

Compute position given constant velocity:

```csharp
double s = 3.0.PositionFromConstantVelocity(time: 5, initialPosition: 2); 

Vector velocity = new Vector(2, 0, 0);
Vector initialPosition = new Vector(1, 0, 0);
Vector position = velocity.PositionFromConstantVelocity(4, initialPosition);
```

**Constant Acceleration**

Compute velocity or position under constant acceleration:

```csharp
// Scalar
double v = 2.0.VelocityFromConstantAcceleration(time: 4, initialVelocity: 1);
double s = 2.0.PositionFromConstantAcceleration(time: 3, initialVelocity: 1, initialPosition: 0);

// Vector
Vector a = new Vector(1, 0, 0);
Vector v0 = new Vector(0, 0, 0);
Vector s0 = new Vector(0, 0, 0);

Vector vVec = a.VelocityFromConstantAcceleration(3, v0);
Vector sVec = a.PositionFromConstantAcceleration(2, v0, s0);
```

**Time-Independent (SUVAT) Equations**

Compute kinematics without explicit time:

```csharp
double finalV = 2.0.VelocityFromDisplacement(5, initialVelocity: 1);
double displacement = 2.0.DisplacementFromVelocities(finalVelocity: 5, initialVelocity: 1);
double t = 2.0.TimeToReachVelocity(finalVelocity: 5, initialVelocity: 1);
double sAvg = 3.0.DisplacementFromAverageVelocity(initialVelocity: 2, finalVelocity: 4);

// Vector versions are also supported
Vector a = new Vector(2, 0, 0);
Vector s = new Vector(3, 0, 0);
Vector v0 = new Vector(1, 0, 0);

Vector finalVVec = a.VelocityFromDisplacement(s, v0);
Vector displacementVec = a.DisplacementFromVelocities(finalVVec, v0);
```

**Circular Motion**

Centripetal acceleration:

```csharp
double a = 3.0.CentripetalAcceleration(radius: 2); // scalar

Vector velocity = new Vector(2, 0, 0);
Vector radius = new Vector(0, 3, 0);
Vector ac = velocity.CentripetalAcceleration(radius); // vector, towards center
```

Angular speed, period, and frequency:

```csharp
double omega = 10.0.AngularSpeed(radius: 5);  // ω = v/r = 2 rad/s
double T     = 10.0.Period(radius: 5);         // T = 2πr/v
double f     = 10.0.Frequency(radius: 5);      // f = v/(2πr) = 1/T
```

Angular velocity vector and tangential velocity:

```csharp
// Object at (5,0,0) moving at (0,10,0) → angular velocity along +Z
Vector omega = vel.AngularVelocity(radius);    // ω = (r × v) / |r|²

// Reverse: angular velocity → tangential velocity
Vector v = omega.TangentialVelocity(radius);   // v = ω × r
```

**Projectile Motion**

Create an initial velocity vector from speed and launch angle:

```csharp
// 20 m/s at 45° → v₀ = (v·cos(θ), 0, v·sin(θ))
Vector v0 = 20.0.ProjectileVelocityFromAngle(Math.PI / 4);
```

Position and velocity at any time:

```csharp
Vector pos = v0.ProjectilePosition(time: 1.5);
Vector vel = v0.ProjectileVelocity(time: 1.5);

// With initial height (e.g. launched from a 10m cliff)
Vector pos2 = v0.ProjectilePosition(time: 1.5, initialHeight: 10);
```

Time of flight, maximum height, and range:

```csharp
// Vector versions (support initial height)
double T = v0.ProjectileTimeOfFlight();
double H = v0.ProjectileMaxHeight();
double R = v0.ProjectileRange();

// With elevated launch
double T2 = v0.ProjectileTimeOfFlight(initialHeight: 10);
double R2 = v0.ProjectileRange(initialHeight: 10);

// Scalar versions (speed + angle, same-height launch)
double T3 = 20.0.ProjectileTimeOfFlight(Math.PI / 4);   // T = 2v₀sin(θ)/g
double H3 = 20.0.ProjectileMaxHeight(Math.PI / 4);      // H = v₀²sin²(θ)/(2g)
double R3 = 20.0.ProjectileRange(Math.PI / 4);           // R = v₀²sin(2θ)/g
```

**Orbital Mechanics**

Gravitational helpers:

```csharp
// Gravitational field strength at distance r from a mass: g = GM/r²
double g = PhysicsConstants.EarthMass.GravitationalFieldStrength(PhysicsConstants.EarthRadius);

// Gravitational force between two masses: F = G·m₁·m₂/r²
double F = PhysicsConstants.EarthMass.GravitationalForce(PhysicsConstants.MoonMass, 3.844e8);

// Escape velocity: v = √(2GM/r)
double vEsc = PhysicsConstants.EarthMass.EscapeVelocity(PhysicsConstants.EarthRadius);
```

Circular orbit scalars:

```csharp
double r = PhysicsConstants.EarthRadius + 408000; // ISS altitude

double speed  = PhysicsConstants.EarthMass.OrbitalSpeed(r);  // v = √(GM/r)  ≈ 7660 m/s
double period = PhysicsConstants.EarthMass.OrbitalPeriod(r);  // T = 2π√(r³/GM) ≈ 92 min
```

Position, velocity, and acceleration on a circular orbit at time t:

```csharp
double M = PhysicsConstants.EarthMass;
double r = 1e7; // 10 000 km radius

Vector pos = M.OrbitalPosition(r, time);      // R·(cos ωt, sin ωt, 0)
Vector vel = M.OrbitalVelocity(r, time);      // Rω·(-sin ωt, cos ωt, 0)
Vector acc = M.OrbitalAcceleration(r, time);   // -ω²R·(cos ωt, sin ωt, 0)
```

**Relative Motion**

Basic relative kinematics between two objects or reference frames:

```csharp
var vA = new Vector(30, 0, 0);
var vB = new Vector(-20, 0, 0);

Vector vRel = vA.RelativeVelocity(vB);        // v_A - v_B = (50, 0, 0)
Vector rRel = posA.RelativePosition(posB);     // r_A - r_B
Vector aRel = accA.RelativeAcceleration(accB); // a_A - a_B
```

Closing speed (positive = approaching, negative = separating):

```csharp
double cs = vA.ClosingSpeed(vB, posA, posB);
```

Position in a moving reference frame over time:

```csharp
Vector relPos = vObj.RelativePositionAtTime(vRef, time: 5.0, r0Obj, r0Ref);
```

Closest approach between two objects at constant velocity:

```csharp
double t   = vA.TimeOfClosestApproach(vB, posA, posB);
double d   = vA.MinimumDistance(vB, posA, posB);
```

---

## ⚡ Dynamics Extensions

The `DynamicsExtensions` class provides extension methods for particle dynamics — **forces, momentum, energy, work, power, and collisions**. This is Phase 1.1 of the [physics roadmap](../../ROADMAP.md).

**Newton's Laws**

```csharp
// Newton's 2nd law: a = F/m
var force = new Vector(10, 0, 0);
Vector acceleration = force.Acceleration(mass: 5);    // (2, 0, 0)

// F = ma
Vector F = 10.0.Force(new Vector(0, 0, -9.8));        // (0, 0, -98)

// Sum of forces
Vector net = f1.NetForce(f2, f3);

// Weight (default direction: -Z)
Vector W = 80.0.Weight();                              // (0, 0, -784.5)
```

**Momentum & Impulse**

```csharp
Vector p = 5.0.Momentum(new Vector(3, 4, 0));         // (15, 20, 0)
Vector J = force.Impulse(duration: 0.5);               // F·Δt
Vector J2 = mass.ImpulseFromVelocityChange(vBefore, vAfter); // m·Δv
Vector vNew = impulse.ApplyImpulse(mass: 5, v0);       // v + J/m
```

**Energy**

```csharp
double ke = 4.0.KineticEnergy(new Vector(3, 4, 0));   // ½mv² = 50 J
double pe = 10.0.PotentialEnergy(height: 5);           // mgh
double U  = m1.GravitationalPotentialEnergy(m2, r);    // -Gm₁m₂/r
double E  = mass.MechanicalEnergy(velocity, height);   // KE + PE
double v  = mass.SpeedFromKineticEnergy(ke);            // √(2·KE/m)
```

**Work & Power**

```csharp
double W = force.Work(displacement);                   // F·d (dot product)
double W2 = 10.0.Work(5, angleRadians: Math.PI / 3);  // F·d·cos(θ) = 25 J
double P = force.Power(velocity);                      // F·v (instantaneous)
double P2 = 100.0.AveragePower(duration: 5);           // W/Δt = 20 W
double ΔKE = mass.WorkEnergyTheorem(vBefore, vAfter);  // ½m(v₂² - v₁²)
```

**Elastic & Inelastic Collisions**

```csharp
// 1D elastic collision — both momentum and energy conserved
var (v1f, v2f) = m1.ElasticCollision(v1, m2, v2);

// Perfectly inelastic (sticky) collision
double vf = m1.InelasticCollisionVelocity(v1, m2, v2);
Vector vf3d = m1.InelasticCollisionVelocity(v1Vec, m2, v2Vec);   // 3D version
double loss = m1.InelasticCollisionEnergyLoss(v1, m2, v2);

// Coefficient of restitution: e = 1 (elastic), e = 0 (perfectly inelastic)
double e = v1Before.CoefficientOfRestitution(v2Before, v1After, v2After);
```

**Rigid Body** (`Physics/Objects/RigidBody.cs`)

Create rigid bodies from standard shapes with automatic inertia tensors:

```csharp
var sphere  = RigidBody.CreateSolidSphere(mass: 10, radius: 2);
var box     = RigidBody.CreateSolidBox(mass: 12, width: 2, height: 3, depth: 4);
var cyl     = RigidBody.CreateSolidCylinder(mass: 6, radius: 1, height: 4);
var hollow  = RigidBody.CreateHollowSphere(mass: 5, radius: 3);
var tube    = RigidBody.CreateHollowCylinder(mass: 8, innerRadius: 1, outerRadius: 2, height: 3);
var rod     = RigidBody.CreateThinRod(mass: 4, length: 5);
var wall    = RigidBody.CreateStatic(new Vector(0, 0, 0));  // immovable
```

Apply forces, torques, and query state:

```csharp
var body = RigidBody.CreateSolidSphere(10, 1);
body.Position = new Vector(0, 5, 0);
body.Velocity = new Vector(3, 0, 0);

body.ApplyForce(new Vector(0, 0, -98));           // gravity
body.ApplyForceAtPoint(                            // generates torque too
    new Vector(10, 0, 0), worldPoint: new Vector(0, 1, 0));
body.ApplyTorque(new Vector(0, 0, 5));

Vector a     = body.LinearAcceleration;            // F/m
Vector alpha = body.AngularAcceleration;           // I⁻¹τ
Vector p     = body.LinearMomentum;                // mv
Vector L     = body.AngularMomentum;               // Iω
double KE    = body.KineticEnergy;                 // ½mv² + ½ωᵀIω

body.ClearForces();  // reset accumulators after integration step
```

**Moment of Inertia (Scalar)**

```csharp
double I = 10.0.MomentOfInertiaSolidSphere(radius: 2);       // 2/5·mr²
double I2 = 10.0.MomentOfInertiaHollowSphere(radius: 2);     // 2/3·mr²
double I3 = 8.0.MomentOfInertiaSolidCylinder(radius: 3);      // ½mr²
double I4 = 6.0.MomentOfInertiaThinRod(length: 4);            // mL²/12
double I5 = 6.0.MomentOfInertiaThinRodEnd(length: 4);         // mL²/3
double I6 = 12.0.MomentOfInertiaSolidBox(sideA: 3, sideB: 4); // m/12·(a²+b²)

// Parallel axis theorem: I_new = I_cm + m·d²
double Inew = I.ParallelAxis(mass: 10, distance: 3);
```

**Inertia Tensor (3×3 Matrix)**

```csharp
Matrix I  = 10.0.InertiaTensorSolidSphere(radius: 2);
Matrix Ib = 12.0.InertiaTensorSolidBox(width: 2, height: 3, depth: 4);
Matrix Ic = 6.0.InertiaTensorSolidCylinder(radius: 1, height: 4);

// Parallel axis theorem for 3×3 tensor: I_new = I_cm + m·(d²E - d⊗d)
Matrix Inew = I.ParallelAxis(mass: 10, offset: new Vector(3, 0, 0));
```

**Torque & Rotational Dynamics**

```csharp
// Vector: τ = r × F
Vector tau = momentArm.Torque(force);

// Scalar: τ = F·r·sin(θ)
double tau2 = 20.0.Torque(momentArm: 3, angleRadians: Math.PI / 6);

// Angular momentum: L = Iω
Vector L = inertiaTensor.AngularMomentum(omega);
double Ls = momentOfInertia.AngularMomentum(omega);

// Angular acceleration: α = I⁻¹τ
Vector alpha = inverseInertiaTensor.AngularAcceleration(torque);

// Rotational kinetic energy: KE = ½ωᵀIω
double KE = inertiaTensor.RotationalKineticEnergy(omega);
double KEs = momentOfInertia.RotationalKineticEnergy(omega);
```

**RigidBody Integration**

Three time-stepping methods advance a `RigidBody` by one `dt`. Each delegates to the corresponding ODE solver in `DifferentialEquationExtensions`, packing body state into `VectorN` internally — no duplicated stepping logic.

```csharp
// Semi-implicit (symplectic) Euler — stable for games, first-order
// Apply forces before calling; accumulators are cleared afterwards
var body = RigidBody.CreateSolidSphere(10, 1);
body.Position = new Vector(0, 0, 100);
body.ApplyForce(new Vector(0, 0, -9.80665 * body.Mass));
body.IntegrateSemiImplicitEuler(dt: 0.001);

// Velocity Verlet — O(dt²) accuracy, excellent energy conservation
// Forces are evaluated by forceFunc (no need to pre-apply)
Func<RigidBody, (Vector force, Vector torque)> gravity = b =>
    (new Vector(0, 0, -9.80665 * b.Mass), new Vector(0, 0, 0));
body.IntegrateVelocityVerlet(gravity, dt: 0.01);

// Explicit Euler — simplest, least stable
body.ApplyForce(new Vector(0, 0, -9.80665 * body.Mass));
body.IntegrateEuler(dt: 0.001);
```

**Common Force Models**

Ready-made force functions for building simulations — combine them in a `forceFunc` for the integrators.

```csharp
// Spring (Hooke's law): F = -k·(|Δr| - L₀)·r̂
var spring = 50.0.SpringForce(restLength: 2.0, position, anchor);

// Viscous damping: F = -c·v
var damping = 0.5.DampingForce(velocity);

// Aerodynamic drag: F = -½·Cd·ρ·A·|v|·v
var drag = 0.47.DragForce(fluidDensity: 1.225, crossSectionArea: 0.01, velocity);

// Friction (kinetic — moving object)
var kinetic = 0.3.FrictionForce(normalForceMagnitude: 98, velocity);

// Friction (static — stationary object, opposes applied force up to μN)
var friction = 0.5.FrictionForce(normalForceMagnitude: 100, velocity: new Vector(0, 0, 0),
    appliedTangentialForce: new Vector(20, 0, 0));  // returns (-20, 0, 0)

// Compose forces in a Verlet simulation
Func<RigidBody, (Vector, Vector)> forces = b =>
{
    var F = 50.0.SpringForce(2.0, b.Position, new Vector(0, 0, 0))
          + 0.5.DampingForce(b.Velocity);
    return (F, new Vector(0, 0, 0));
};
body.IntegrateVelocityVerlet(forces, dt: 0.01);
```

---



## 🔭 Astronomy Extensions

The `AstronomyExtensions` class provides extension methods for astronomical calculations: **distance conversions**, **Julian date**, **sidereal time**, and **horizontal ↔ equatorial coordinate transforms**. No external libraries are used.

**Distance Conversions**

Convert between light-years, parsecs, and astronomical units:

```csharp
double pc = 4.37.LightYearsToParsecs();      // Proxima Centauri ≈ 1.34 pc
double ly = 1.0.ParsecsToLightYears();        // 1 pc ≈ 3.26 ly
double au = 1.0.LightYearsToAU();             // 1 ly ≈ 63241 AU
double au2 = 1.0.ParsecsToAU();               // 1 pc ≈ 206265 AU
```

**Julian Date**

Compute Julian Date and Julian centuries since J2000.0:

```csharp
var utc = new DateTime(2024, 6, 15, 21, 0, 0, DateTimeKind.Utc);

double jd = utc.ToJulianDate();               // Julian Date
double T  = utc.JulianCenturiesSinceJ2000();  // centuries since J2000.0
```

**Sidereal Time**

Compute Greenwich and local sidereal time from UTC, or from local time with time zone:

```csharp
var utc = new DateTime(2024, 6, 15, 21, 0, 0, DateTimeKind.Utc);

double gmst = utc.GreenwichMeanSiderealTimeHours();       // GMST in hours
double gmstDeg = utc.GreenwichMeanSiderealTimeDegrees();   // GMST in degrees

// Local sidereal time (Stockholm: 18.07° E)
double lmst = utc.LocalMeanSiderealTimeHours(18.07);

// From local time + time zone + longitude
var local = new DateTime(2024, 6, 15, 23, 0, 0);
double lst = local.LocalSiderealTimeFromLocal(
    utcOffsetHours: 2.0,       // UTC+2 (CEST)
    longitudeDegrees: 18.07);  // Stockholm
```

**Horizontal → Equatorial (Altitude/Azimuth → RA/Dec)**

Convert what you observe (altitude, azimuth) to sky coordinates (right ascension, declination):

```csharp
// Using explicit local sidereal time
var (ra, dec) = AstronomyExtensions.HorizontalToEquatorial(
    altitudeDegrees: 45.0,
    azimuthDegrees: 180.0,     // due south
    latitudeDegrees: 51.48,    // London
    localSiderealTimeDegrees: 120.0);

// Using UTC time and longitude (LST computed automatically)
var utc = new DateTime(2024, 3, 20, 22, 0, 0, DateTimeKind.Utc);
var (ra2, dec2) = AstronomyExtensions.HorizontalToEquatorial(
    altitudeDegrees: 45.0,
    azimuthDegrees: 180.0,
    latitudeDegrees: 51.48,
    longitudeDegrees: 0.0,     // Greenwich
    utc: utc);
```

**Equatorial → Horizontal (RA/Dec → Altitude/Azimuth)**

Find where a star appears in the sky from your location:

```csharp
var (alt, az) = AstronomyExtensions.EquatorialToHorizontal(
    rightAscensionDegrees: 200.0,
    declinationDegrees: 30.0,
    latitudeDegrees: 59.33,    // Stockholm
    localSiderealTimeDegrees: 120.0);

// Or with UTC + position
var (alt2, az2) = AstronomyExtensions.EquatorialToHorizontal(
    rightAscensionDegrees: 200.0,
    declinationDegrees: 30.0,
    latitudeDegrees: 59.33,
    longitudeDegrees: 18.07,
    utc: utc);
```

**Angle Helpers**

Convert between common astronomical angle formats:

```csharp
// Right ascension: hours/min/sec ↔ degrees
double deg = AstronomyExtensions.RightAscensionToDegrees(6, 30, 0);  // → 97.5°
var (h, m, s) = AstronomyExtensions.DegreesToRightAscension(97.5);    // → (6, 30, 0.0)

// Declination: degrees/arcmin/arcsec → decimal degrees
double dec = AstronomyExtensions.DeclinationToDegrees(-16, 42, 58);   // → -16.7161°
```

### Transit Geometry

Compute geometric properties of planetary transits:

```csharp
using CSharpNumerics.Physics.Astro;

// Impact parameter: b = (a/R★) · cos(i)
double b = TransitGeometry.ImpactParameter(aOverRstar: 15.0, inclination: Math.PI / 2 - 0.01);

// Geometric transit probability: P = (R★ + Rp) / a
double prob = TransitGeometry.TransitProbability(a: 0.05, rStar: 0.01, rPlanet: 0.001);

// Total transit duration T14 (days)
double t14 = TransitGeometry.TransitDuration(
    period: 3.0, aOverRstar: 15.0, radiusRatio: 0.1, inclination: Math.PI / 2);

// Ingress/egress duration T12
double t12 = TransitGeometry.IngressDuration(
    period: 3.0, aOverRstar: 15.0, radiusRatio: 0.1, inclination: Math.PI / 2);

// Contact times T1–T4
var (t1, t2, t3, t4) = TransitGeometry.ContactTimes(
    period: 3.0, aOverRstar: 15.0, radiusRatio: 0.1, inclination: Math.PI / 2, epoch: 100.0);
```

### Limb Darkening

Stellar limb darkening laws — intensity as a function of $\mu = \cos\theta$ (angle from disk centre):

```csharp
using CSharpNumerics.Physics.Astro;

double mu = 0.5;

// Linear: I(μ) = 1 − u1·(1 − μ)
double iLinear = LimbDarkening.Linear(mu, u1: 0.6);

// Quadratic: I(μ) = 1 − u1·(1−μ) − u2·(1−μ)²
double iQuad = LimbDarkening.Quadratic(mu, u1: 0.4, u2: 0.2);

// Nonlinear four-parameter (Claret 2000)
double iNl = LimbDarkening.NonlinearFourParam(mu, c1: 0.5, c2: -0.2, c3: 0.3, c4: -0.1);

// Intensity profile for an array of μ values
double[] muArray = { 0.0, 0.25, 0.5, 0.75, 1.0 };
double[] profile = LimbDarkening.IntensityProfile(
    LimbDarkeningModel.Quadratic, new[] { 0.4, 0.2 }, muArray);
```

| Model | Enum | Coefficients |
|-------|------|-------------|
| Uniform | `LimbDarkeningModel.Uniform` | — |
| Linear | `LimbDarkeningModel.Linear` | u1 |
| Quadratic | `LimbDarkeningModel.Quadratic` | u1, u2 |
| Nonlinear 4-param | `LimbDarkeningModel.NonlinearFourParam` | c1, c2, c3, c4 |

### Transit Model (Mandel & Agol 2002)

Analytical transit light curve — computes fractional flux drop as a planet crosses a limb-darkened stellar disk. Supports circular orbits with all limb darkening models.

```csharp
using CSharpNumerics.Physics.Astro;
using CSharpNumerics.Engines.Exoplanet.Data;

var model = new TransitModel();

// With stellar properties (a/R★ computed from Kepler's third law)
var p = new TransitParameters
{
    Period = 3.0,
    Epoch = 100.5,
    RadiusRatio = 0.1,
    ImpactParameter = 0.3,
    Duration = 0.15
};
var star = new StellarProperties { Radius = 1.0, Mass = 1.0 };

double[] times = /* observation times in days (BJD) */;
double[] flux = model.Evaluate(times, p, LimbDarkeningModel.Quadratic, new[] { 0.4, 0.2 }, star);
// flux[i] ≈ 1.0 out of transit, < 1.0 during transit

// Or specify orbital parameters directly (no stellar properties needed)
double[] flux2 = model.Evaluate(times,
    period: 3.0, epoch: 100.5, radiusRatio: 0.1,
    aOverRstar: 15.0, inclination: Math.PI / 2,
    LimbDarkeningModel.Quadratic, new[] { 0.4, 0.2 });
```

### Kepler Orbit

Orbital mechanics: Kepler's equation and third law.

```csharp
using CSharpNumerics.Physics.Astro;

// Solve Kepler's equation: M = E − e·sin(E) → true anomaly ν
double nu = KeplerOrbit.TrueAnomaly(meanAnomaly: 1.0, eccentricity: 0.3);

// Kepler's third law: a = (G·M★·P²/(4π²))^(1/3)
double a = KeplerOrbit.SemiMajorAxis(period: 259200, stellarMass: 1.989e30); // SI units
double aAU = KeplerOrbit.SemiMajorAxisAU(periodDays: 365.25, stellarMassSolar: 1.0); // ≈ 1 AU

// Mean orbital velocity: v = 2πa/P
double v = KeplerOrbit.OrbitalVelocity(a: 1.496e11, period: 3.156e7); // ≈ 29.8 km/s
```

### Exoplanet Classification

Classify stars by temperature, compute habitable zones, and measure how Earth-like a planet is — all in `AstronomyExtensions`.

**Spectral Classification**

Map a stellar effective temperature to its Harvard spectral type (O B A F G K M L T Y):

```csharp
using CSharpNumerics.Physics.Astro;
using CSharpNumerics.Physics.Astro.Enums;

SpectralType sun   = AstronomyExtensions.GetSpectralFromTemp(5778);  // G
SpectralType hot   = AstronomyExtensions.GetSpectralFromTemp(30000); // O
SpectralType cool  = AstronomyExtensions.GetSpectralFromTemp(3000);  // M
SpectralType brown = AstronomyExtensions.GetSpectralFromTemp(1800);  // L
```

| Temperature (K) | Spectral Type |
|------------------|---------------|
| ≥ 30 000 | O |
| 10 000 – 29 999 | B |
| 7 500 – 9 999 | A |
| 6 000 – 7 499 | F |
| 5 200 – 5 999 | G |
| 3 700 – 5 199 | K |
| 2 400 – 3 699 | M |
| 1 300 – 2 399 | L |
| 550 – 1 299 | T |
| < 550 | Y |

**Habitable Zone (Goldilocks Zone)**

Compute the conservative habitable zone boundaries using the Kopparapu et al. (2013) parameterisation:

```csharp
// From luminosity (solar units) + effective temperature
var (inner, outer) = AstronomyExtensions.CalculateGoldilocksZone(
    stellarLuminositySolar: 1.0,
    effectiveTemperatureK: 5778);
// inner ≈ 0.99 AU, outer ≈ 1.69 AU — Earth sits comfortably inside

// From radius (solar radii) + temperature (luminosity derived via Stefan–Boltzmann)
var (inner2, outer2) = AstronomyExtensions.CalculateGoldilocksZoneFromRadius(
    stellarRadiusSolar: 1.0,
    effectiveTemperatureK: 5778);

// Red dwarf (Proxima Centauri-like)
var (innerM, outerM) = AstronomyExtensions.CalculateGoldilocksZone(0.04, 3200);
// innerM ≈ 0.19 AU, outerM ≈ 0.35 AU
```

**Earth Similarity Index (ESI)**

Quantify how Earth-like a planet is on a 0–1 scale (Schulze-Makuch et al. 2011). Uses four parameters normalised to Earth values:

```csharp
// Earth (reference) → ESI = 1.0
double esiEarth = AstronomyExtensions.CalculateEsi(
    radiusEarth: 1.0,
    densityEarth: 1.0,
    escapeVelocityEarth: 1.0,
    surfaceTemperatureK: 288);   // 1.0

// Mars: R≈0.53, ρ≈0.71, vesc≈0.45, T≈210 K
double esiMars = AstronomyExtensions.CalculateEsi(0.53, 0.71, 0.45, 210);  // ≈ 0.73

// Venus: R≈0.95, ρ≈0.95, vesc≈0.93, T≈737 K
double esiVenus = AstronomyExtensions.CalculateEsi(0.95, 0.95, 0.93, 737); // ≈ 0.44
```

| Parameter | Earth value | Weight |
|-----------|-------------|--------|
| Radius | 1.0 R⊕ | 0.57 |
| Bulk density | 1.0 ρ⊕ | 1.07 |
| Escape velocity | 1.0 v⊕ | 0.70 |
| Surface temperature | 288 K | 5.58 |

---

## 🔄 Oscillations

The `Physics.Oscillations` namespace provides one-dimensional oscillator models with both **analytic** and **numerical** solutions. All oscillators implement `IOscillator` for a consistent API.

### Simple Harmonic Oscillator

Models the undamped system $\ddot{x} + \omega_0^2 x = 0$ where $\omega_0 = \sqrt{k/m}$.

Integration uses **Velocity Verlet** (symplectic, O(dt²)) — excellent energy conservation for undamped systems.

**Create and inspect properties**

```csharp
using CSharpNumerics.Physics.Oscillations;

// mass = 2 kg, spring constant = 50 N/m, initial displacement = 0.5 m
var sho = new SimpleHarmonicOscillator(mass: 2, stiffness: 50, initialPosition: 0.5);

double w = sho.AngularFrequency;   // ω₀ = √(k/m) = 5 rad/s
double f = sho.Frequency;          // f₀ = ω₀/2π ≈ 0.796 Hz
double T = sho.Period;             // T = 1/f₀ ≈ 1.257 s
double A = sho.Amplitude;          // A = √(x₀² + (v₀/ω₀)²) = 0.5 m
double phi = sho.PhaseOffset;      // φ = atan2(-v₀/ω₀, x₀)
```

**Analytic solution** — exact reference for verification

```csharp
double x = sho.AnalyticPosition(t: 1.0);   // x(t) = A·cos(ω₀t + φ)
double v = sho.AnalyticVelocity(t: 1.0);   // v(t) = -Aω₀·sin(ω₀t + φ)
```

**Step-by-step simulation**

```csharp
double dt = 0.001;
for (int i = 0; i < 10000; i++)
{
    sho.Step(dt);
    Console.WriteLine($"t={sho.Time:F3}  x={sho.Position:F6}  v={sho.Velocity:F6}");
}

sho.Reset();  // restore initial conditions
```

**Trajectory** — displacement vs time as a `List<Serie>`

```csharp
// Simulate 5 full periods, returns (Index = time, Value = position)
List<Serie> trajectory = sho.Trajectory(tEnd: sho.Period * 5, dt: 0.001);
```

**Phase portrait** — (position, velocity) plot

```csharp
// Returns (Index = x, Value = v) — traces an ellipse for SHM
List<Serie> phase = sho.PhasePortrait(tEnd: sho.Period * 2, dt: 0.001);
```

**Energy**

```csharp
double KE = sho.KineticEnergy;     // ½mv²
double PE = sho.PotentialEnergy;    // ½kx²
double E  = sho.TotalEnergy;       // KE + PE (conserved for undamped SHM)

// Energy over time — verify conservation
List<Serie> energy = sho.EnergyOverTime(tEnd: sho.Period * 10, dt: 0.001);
```

**Frequency spectrum** — FFT of the displacement signal

```csharp
// Returns frequency (Hz) vs normalised magnitude
// Peak appears at the natural frequency f₀
List<Serie> spectrum = sho.FrequencySpectrum(tEnd: sho.Period * 20, dt: 0.01);
```

**With initial velocity**

```csharp
// Released from equilibrium with an impulse
var sho2 = new SimpleHarmonicOscillator(mass: 1, stiffness: 25,
    initialPosition: 0, initialVelocity: 5);

// A = √(0² + (5/5)²) = 1 m
Console.WriteLine(sho2.Amplitude);  // 1.0
```

### Damped Oscillator

Models the damped system $\ddot{x} + 2\gamma\dot{x} + \omega_0^2 x = 0$ where $\gamma = c/(2m)$.

Three damping regimes are automatically detected:
- **Underdamped** ($\gamma < \omega_0$): oscillates with exponentially decaying amplitude
- **Critically damped** ($\gamma = \omega_0$): returns to equilibrium as fast as possible without oscillating
- **Overdamped** ($\gamma > \omega_0$): decays without oscillation, slower than critical

Integration uses **RK4** (4th-order Runge-Kutta) — appropriate for dissipative systems where energy is not conserved.

**Create and inspect properties**

```csharp
using CSharpNumerics.Physics.Oscillations;

// mass = 1 kg, k = 100 N/m, damping coefficient c = 4 N·s/m
var osc = new DampedOscillator(mass: 1, stiffness: 100, damping: 4, initialPosition: 1.0);

double w0 = osc.NaturalFrequency;    // ω₀ = √(k/m) = 10 rad/s
double g  = osc.Gamma;               // γ = c/(2m) = 2 rad/s
double wd = osc.DampedFrequency;     // ω_d = √(ω₀² − γ²) = √96 ≈ 9.80 rad/s
double Td = osc.DampedPeriod;        // T_d = 2π/ω_d ≈ 0.641 s

DampingRegime regime = osc.Regime;   // Underdamped
double Q     = osc.QualityFactor;    // Q = ω₀/(2γ) = 2.5
double delta = osc.LogarithmicDecrement; // δ = γ·T_d ≈ 1.283
```

**Damping regime detection**

```csharp
var under   = new DampedOscillator(1, 100,  4);  // γ=2  < ω₀=10 → Underdamped
var crit    = new DampedOscillator(1, 100, 20);  // γ=10 = ω₀=10 → CriticallyDamped
var over    = new DampedOscillator(1, 100, 40);  // γ=20 > ω₀=10 → Overdamped
```

**Analytic solution** — exact reference for all three regimes

```csharp
// Underdamped:  x(t) = A·e^(−γt)·cos(ω_d·t + φ)
// Critical:     x(t) = (C₁ + C₂t)·e^(−γt)
// Overdamped:   x(t) = C₁·e^(r₁t) + C₂·e^(r₂t)
double x = osc.AnalyticPosition(t: 1.0);
double v = osc.AnalyticVelocity(t: 1.0);
```

**Step-by-step simulation**

```csharp
double dt = 0.001;
for (int i = 0; i < 5000; i++)
{
    osc.Step(dt);
    Console.WriteLine($"t={osc.Time:F3}  x={osc.Position:F6}  v={osc.Velocity:F6}");
}

osc.Reset();  // restore initial conditions
```

**Trajectory and phase portrait**

```csharp
// Displacement vs time
List<Serie> trajectory = osc.Trajectory(tEnd: 5.0, dt: 0.001);

// Phase portrait — spirals to origin for underdamped
List<Serie> phase = osc.PhasePortrait(tEnd: 5.0, dt: 0.001);
```

**Energy** — monotonically decreasing for damped systems

```csharp
double E = osc.TotalEnergy;         // KE + PE (decreases over time)

// Energy over time
List<Serie> energy = osc.EnergyOverTime(tEnd: 5.0, dt: 0.001);

// Energy dissipated by damping: E₀ − E(t), monotonically increasing
List<Serie> dissipated = osc.EnergyDissipation(tEnd: 5.0, dt: 0.001);
```

**Envelope** — the exponential decay bounding the oscillation

```csharp
double env = osc.Envelope(t: 1.0);  // A·e^(−γt) at a given time

// Envelope curve over time
List<Serie> envCurve = osc.EnvelopeCurve(tEnd: 5.0, dt: 0.01);
```

**Frequency spectrum** — peak at ω_d with spectral broadening proportional to γ

```csharp
List<Serie> spectrum = osc.FrequencySpectrum(tEnd: 20.0, dt: 0.005);
```

**Zero damping** — reduces to `SimpleHarmonicOscillator` behaviour

```csharp
var undamped = new DampedOscillator(mass: 1, stiffness: 25, damping: 0);
// Regime = Underdamped, QualityFactor = ∞, DampedFrequency = ω₀
```

### Driven Oscillator

Extends the damped oscillator with a harmonic driving force:

$$\ddot{x} + 2\gamma\dot{x} + \omega_0^2 x = \frac{F_0}{m}\cos(\omega_d t)$$

**Constructor**

```csharp
var osc = new DrivenOscillator(
    mass:           1.0,
    stiffness:      100.0,
    damping:        4.0,
    driveAmplitude: 10.0,
    driveFrequency: 7.0,
    initialPosition: 0.0,   // optional, default 0
    initialVelocity: 0.0    // optional, default 0
);
```

| Parameter | Description |
|---|---|
| `mass` | Mass $m$ (must be > 0) |
| `stiffness` | Spring constant $k$ (must be ≥ 0) |
| `damping` | Damping coefficient $c$ (must be ≥ 0) |
| `driveAmplitude` | Force amplitude $F_0$ (must be ≥ 0) |
| `driveFrequency` | Angular frequency $\omega_d$ of the driving force (must be ≥ 0) |

**Physical properties**

| Property | Formula | Description |
|---|---|---|
| `NaturalFrequency` | $\omega_0 = \sqrt{k/m}$ | Undamped natural frequency |
| `Gamma` | $\gamma = c/(2m)$ | Damping rate |
| `DampedFrequency` | $\omega_d = \sqrt{\omega_0^2 - \gamma^2}$ | Damped frequency |
| `Regime` | auto-detected | `Underdamped`, `CriticallyDamped`, or `Overdamped` |
| `QualityFactor` | $Q = \omega_0/(2\gamma)$ | Sharpness of resonance |
| `ResonanceFrequency` | $\omega_r = \sqrt{\omega_0^2 - 2\gamma^2}$ | Frequency of maximum amplitude response (0 if $\omega_0^2 < 2\gamma^2$) |
| `Bandwidth` | $\Delta\omega = 2\gamma$ | Width of resonance peak |

**Steady-state response**

The forced response after transients decay:

```csharp
// Amplitude & phase at a specific driving frequency
double A  = osc.SteadyStateAmplitude(wd);    // A(ωd) = (F₀/m) / √((ω₀²−ωd²)² + (2γωd)²)
double φ  = osc.SteadyStatePhase(wd);        // φ(ωd) = −atan2(2γωd, ω₀²−ωd²)

// At the configured drive frequency
double Adrive = osc.SteadyStateAmplitudeAtDrive;
double φdrive = osc.SteadyStatePhaseAtDrive;
```

**Resonance curve & phase response**

```csharp
// Sweep amplitude A(ω) over a frequency range
List<Serie> resonance = osc.ResonanceCurve(wMin: 0.1, wMax: 20, steps: 500);

// Sweep phase φ(ω) over a frequency range
List<Serie> phase = osc.PhaseResponse(wMin: 0.1, wMax: 20, steps: 500);
```

**Transfer function & frequency response**

```csharp
// H(s) = 1 / (s² + 2γs + ω₀²) — Laplace domain
ComplexNumber H = osc.TransferFunction(s);

// H(iω) — evaluate on imaginary axis
ComplexNumber Hiw = osc.FrequencyResponse(w);
// |H(iω)| = A(ω) / (F₀/m),   arg(H(iω)) = φ(ω)
```

**Simulation (RK4)**

```csharp
osc.Step(dt: 0.001);                        // advance one step
List<Serie> traj  = osc.Trajectory(5.0, 0.001);
List<Serie> phase = osc.PhasePortrait(5.0, 0.001);
```

**Steady-state extraction**

```csharp
// Simulates past the transient, returns only the steady-state window
List<Serie> steady = osc.SteadyStateTrajectory(
    steadyDuration: 10.0,
    dt: 0.001,
    transientDuration: null  // auto = 5/γ
);

// Measure numerical peak-to-peak amplitude in steady state
double Ameas = osc.MeasuredSteadyStateAmplitude(dt: 0.001);
```

**Energy & power**

```csharp
List<Serie> energy = osc.EnergyOverTime(tEnd: 5.0, dt: 0.001);
List<Serie> power  = osc.PowerInput(tEnd: 5.0, dt: 0.001);
// PowerInput = F(t) · v(t) — instantaneous power from the drive
```

**Frequency spectrum**

```csharp
List<Serie> spectrum = osc.FrequencySpectrum(tEnd: 50.0, dt: 0.005);
// Peak appears at the drive frequency f_d = ω_d / (2π)
```

**Zero driving force** — reduces to `DampedOscillator` behaviour

```csharp
var undriven = new DrivenOscillator(mass: 1, stiffness: 100, damping: 4,
                                     driveAmplitude: 0, driveFrequency: 0, initialPosition: 1);
// Identical trajectory to DampedOscillator with same parameters
```

### Coupled Oscillators

Models an N-mass spring chain with fixed walls on both ends:

$$M\ddot{\mathbf{x}} + C\dot{\mathbf{x}} + K\mathbf{x} = 0$$

where $M$ is a diagonal mass matrix, $K$ is a tridiagonal stiffness matrix, and $C$ is an optional diagonal damping matrix.

```
wall —k₀— m₁ —k₁— m₂ —k₂— … —mₙ —kₙ— wall
```

**General constructor**

```csharp
var osc = new CoupledOscillators(
    masses:      new[] { 1.0, 2.0, 1.5 },
    stiffnesses: new[] { 5.0, 10.0, 8.0, 5.0 },  // N+1 springs
    dampings:    new[] { 0.1, 0.2, 0.1 },          // optional
    initialPositions:  new[] { 1.0, 0.0, -0.5 },   // optional
    initialVelocities: new[] { 0.0, 0.0, 0.0 }     // optional
);
```

**Uniform constructor**

```csharp
var osc = new CoupledOscillators(
    count: 5, mass: 1.0, stiffness: 4.0, damping: 0.0,
    initialPositions: new[] { 1.0, 0.0, 0.0, 0.0, 0.0 }
);
```

| Property | Description |
|---|---|
| `Count` | Number of masses $N$ |
| `Position(i)` / `Velocity(i)` | State of mass $i$ |
| `Positions` / `Velocities` | Copy of all current states |
| `KineticEnergy` | $\frac{1}{2}\sum m_i v_i^2$ |
| `PotentialEnergy` | $\frac{1}{2}\sum k_j (\Delta x_j)^2$ (all springs) |
| `TotalEnergy` | Kinetic + Potential |

**Matrices**

```csharp
Matrix K = osc.StiffnessMatrix();   // N×N tridiagonal symmetric
Matrix M = osc.MassMatrix();        // N×N diagonal
```

**Normal modes** — computed via the Jacobi eigenvalue algorithm on the symmetrised dynamical matrix $D = L^{-1}KL^{-1}$ where $L = \text{diag}(\sqrt{m_i})$

```csharp
List<double> frequencies = osc.NormalModes();   // ω₁ ≤ ω₂ ≤ … ≤ ωₙ
List<VectorN> shapes     = osc.ModeShapes();    // normalised eigenvectors
```

For a uniform 2-mass chain ($m=1$, $k=1$): $\omega_1 = 1$, $\omega_2 = \sqrt{3}$ with modes $[1,1]/\sqrt{2}$ (symmetric) and $[1,-1]/\sqrt{2}$ (antisymmetric).

For a uniform N-mass chain: $\omega_r = 2\sqrt{k/m}\,\sin\!\bigl(\frac{r\pi}{2(N+1)}\bigr)$

**Simulation (RK4)**

```csharp
osc.Step(dt: 0.001);
osc.Reset();
List<List<Serie>> trajectories = osc.Trajectory(tEnd: 5.0, dt: 0.001);
// trajectories[i] = time series for mass i
List<Serie> phase = osc.PhasePortrait(massIndex: 0, tEnd: 5.0, dt: 0.001);
```

**Energy analysis**

```csharp
List<Serie> energy = osc.EnergyOverTime(tEnd: 5.0, dt: 0.001);

// Modal energy decomposition (sum over modes ≈ total energy)
List<Serie> modeEnergy = osc.ModalEnergy(modeIndex: 0, tEnd: 5.0, dt: 0.001);
```

**Dispersion relation** — theoretical for uniform periodic chain: $\omega(k) = 2\sqrt{k/m}\,|\sin(ka/2)|$

```csharp
double[] kValues = Enumerable.Range(0, 100)
    .Select(i => i * Math.PI / 100.0).ToArray();
List<Serie> dispersion = osc.DispersionRelation(kValues, latticeSpacing: 1.0);

double vp = osc.PhaseVelocity(k: 1.0);   // ω/k
double vg = osc.GroupVelocity(k: 1.0);    // dω/dk (numerical derivative)
```

**Frequency spectrum**

```csharp
List<Serie> spectrum = osc.FrequencySpectrum(massIndex: 0, tEnd: 50.0, dt: 0.01);
// Peaks appear at the excited normal mode frequencies
```

### IOscillator Interface

All oscillators share this contract:

```csharp
public interface IOscillator
{
    double Position { get; }
    double Velocity { get; }
    double Time { get; }
    double TotalEnergy { get; }
    double KineticEnergy { get; }
    double PotentialEnergy { get; }

    void Step(double dt);
    void Reset();
    List<Serie> Trajectory(double tEnd, double dt);
    List<Serie> PhasePortrait(double tEnd, double dt);
}
```

This enables polymorphic use — the same analysis code works for `SimpleHarmonicOscillator`, `DampedOscillator`, and `DrivenOscillator`. `CoupledOscillators` has a similar API but operates on N masses simultaneously.

---

## 🌊 Waves

The `Physics.Waves` namespace provides PDE-based wave simulation using **Method of Lines** (spatial discretisation via `GridOperators`, temporal integration via `ITimeStepper`), as well as analytical tools for superposition, wave packets, and Fourier analysis. All wave field classes implement `IWaveField`.

### IWaveField Interface

```csharp
public interface IWaveField
{
    double Time { get; }
    double TotalEnergy { get; }
    void Step(double dt);
    void Reset();
    List<Serie> Snapshot();
}
```

### BoundaryType

Wave-specific boundary conditions that map internally to the finite-difference `BoundaryCondition` enum:

```csharp
public enum BoundaryType { Fixed, Free, Periodic, Absorbing }
```

### WaveEquation1D

Solves the 1D wave equation $\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2}$ via MOL.

State vector: $\mathbf{y} = [u_0 \dots u_{N-1} \mid v_0 \dots v_{N-1}]$. Default stepper: Velocity Verlet (symplectic).

**Create and set initial conditions**

```csharp
using CSharpNumerics.Physics.Waves;

int n = 201;
double dx = 0.005;  // 1 m domain
double c = 340.0;   // speed of sound in air

var wave = new WaveEquation1D(n, dx, c, BoundaryType.Fixed);

// Pluck the first mode
double L = wave.Length;
wave.SetInitialCondition(
    u0: x => Math.Sin(Math.PI * x / L),
    v0: null);
```

**Simulation**

```csharp
// Single step
wave.Step(dt: 0.00001);

// Run to end time, optionally recording trajectory
TimeStepResult result = wave.Simulate(tEnd: 0.01, dt: 0.00001, recordTrajectory: true);

wave.Reset();  // restore initial conditions
```

**Properties and diagnostics**

```csharp
double t  = wave.Time;          // current simulation time
double E  = wave.TotalEnergy;   // conserved for undamped wave
double cfl = wave.CFL(dt);      // Courant number c·dt/dx — must be ≤ 1

double u3 = wave.Displacement(3);  // u at grid point 3
double v3 = wave.Velocity(3);      // ∂u/∂t at grid point 3
```

**Snapshots & field output**

```csharp
// Current displacement as (x, u) pairs
List<Serie> snap = wave.Snapshot();

// Full u(x,t) matrix — rows = spatial, columns = time steps
Matrix field = wave.SpaceTimeField(tEnd: 0.05, dt: 0.0001);

// Energy density at each grid point
List<Serie> eDensity = wave.EnergyDensity();
```

**Frequency analysis & standing wave reference**

```csharp
// FFT of u(x_i, t) at a single spatial point
List<Serie> spectrum = wave.FrequencyContent(spatialIndex: 50, tEnd: 0.1, dt: 0.0001);

// Analytic standing-wave mode shape: sin(nπx/L)
List<Serie> mode = wave.StandingWaveMode(modeNumber: 1);

// Analytic frequency: f_n = n·c/(2L)
double f1 = wave.StandingWaveFrequency(modeNumber: 1);
```

### WaveEquation2D

Solves $\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u$ on a `Grid2D` using the 5-point Laplacian stencil.

```csharp
using CSharpNumerics.Numerics.FiniteDifference;

var grid = new Grid2D(nx: 51, ny: 51, d: 0.02);  // 1 m × 1 m
var wave2d = new WaveEquation2D(grid, c: 1.0, BoundaryType.Fixed);

double Lx = (grid.Nx - 1) * grid.Dx;
double Ly = (grid.Ny - 1) * grid.Dy;

wave2d.SetInitialCondition(
    u0: (x, y) => Math.Sin(Math.PI * x / Lx) * Math.Sin(Math.PI * y / Ly));

wave2d.Simulate(tEnd: 0.5, dt: 0.005);

// 2D snapshot as array[ix, iy]
double[,] snap = wave2d.SnapshotArray();

// Energy density per cell
double[,] eDensity = wave2d.EnergyDensityArray();

double E = wave2d.TotalEnergy;
```

### WaveSuperposition

Analytical superposition of travelling sinusoidal waves — no PDE solver needed:

$$u(x,t) = \sum_i A_i \sin(\omega_i t - k_i x + \varphi_i)$$

```csharp
var ws = new WaveSuperposition();

// Two counter-propagating waves → standing wave
ws.AddHarmonic(amplitude: 1, angularFrequency: 10, waveNumber: 5);
ws.AddHarmonic(amplitude: 1, angularFrequency: 10, waveNumber: -5);

double u = ws.Evaluate(x: 0.5, t: 0.1);

// Field at many x-values
double[] xVals = Enumerable.Range(0, 200).Select(i => i * 0.01).ToArray();
List<Serie> field = ws.EvaluateRange(xVals, t: 0.1);

// Beat frequency for two-component superposition
double fBeat = ws.BeatFrequency();  // |ω₁ − ω₂| / 2π

// Spatial Fourier coefficients via FFT
List<Serie> coeffs = ws.FourierCoefficients(xVals, t: 0);
```

### WavePacket

FFT-based Gaussian wave packet propagation with arbitrary dispersion relation $\omega(k)$:

```csharp
int nPts = 512;
double dxp = 0.05;
double[] x = Enumerable.Range(0, nPts).Select(i => i * dxp).ToArray();

// Non-dispersive: ω = c·k → packet translates without deformation
var packet = new WavePacket(x, centerK: 10.0, sigma: 2.0,
    dispersion: k => 3.0 * k);

double vp = packet.PhaseVelocity;   // ω(k₀)/k₀ = c
double vg = packet.GroupVelocity;   // dω/dk = c  (same for linear)

// Propagate to time t and get displacement
List<Serie> field = packet.Propagate(t: 1.0);

// Measure spatial width (std dev of |u|²)
double w0 = packet.Width(0);
double w1 = packet.Width(1.0);

// Dispersive example: ω = k² → packet spreads
var dispersive = new WavePacket(x, centerK: 10.0, sigma: 2.0,
    dispersion: k => k * k);

double rate = dispersive.SpreadRate(t1: 0, t2: 2.0);
```

### DampedDrivenWaveEquation1D

Solves the damped/driven 1D wave equation:

$$\frac{\partial^2 u}{\partial t^2} + 2\alpha \frac{\partial u}{\partial t} = c^2 \frac{\partial^2 u}{\partial x^2} + S(x, t)$$

When $\alpha = 0$ and $S = 0$ this reduces to the standard wave equation. Default stepper: RK4 (appropriate for dissipative systems).

**Damped wave — energy decay**

```csharp
var damped = new DampedDrivenWaveEquation1D(
    n: 101, dx: 0.01, c: 1.0,
    alpha: 5.0,          // damping coefficient
    boundaryType: BoundaryType.Fixed);

double L = (101 - 1) * 0.01;
damped.SetInitialCondition(u0: x => Math.Sin(Math.PI * x / L));

double e0 = damped.TotalEnergy;
damped.Simulate(tEnd: 1.0, dt: 0.005);
double eFinal = damped.TotalEnergy;   // eFinal << e0
```

**Driven wave — source injection**

```csharp
// Point-like sinusoidal source at the centre
Func<double, double, double> source = (x, t) =>
{
    double xc = 0.5;
    return 10.0 * Math.Sin(20 * t) * Math.Exp(-100 * (x - xc) * (x - xc));
};

var driven = new DampedDrivenWaveEquation1D(
    n: 101, dx: 0.01, c: 1.0,
    alpha: 0,
    source: source);

driven.SetInitialCondition();  // start from rest
driven.Simulate(tEnd: 2.0, dt: 0.005);
// Energy grows as the source injects power
```

**Damped + driven — steady state**

```csharp
var system = new DampedDrivenWaveEquation1D(
    n: 101, dx: 0.01, c: 1.0,
    alpha: 2.0,
    source: source);

system.SetInitialCondition();
system.Simulate(tEnd: 10.0, dt: 0.005);
// Energy reaches bounded steady state — damping balances input
```

**Replace source at runtime**

```csharp
system.SetSource((x, t) => 5.0 * Math.Cos(30 * t) * Math.Exp(-200 * x * x));
```

---

## ⚡ Electromagnetic Field Extensions

The `ElectroMagneticFieldExtensions` class bridges `VectorField` (∇·, ∇×) and `VectorFieldExtensions` (∇, ∇²) to classical electrodynamics — Coulomb's law, Lorentz force, Maxwell's equations, Poynting vector, potentials, and Biot–Savart sources.

### Point Charges

Electric field, potential, and Coulomb force from point charges:

```csharp
double q = 1e-6; // 1 μC
var pos = new Vector(0, 0, 0);
var fieldPoint = new Vector(1, 0, 0);

Vector E = q.ElectricField(pos, fieldPoint);        // E = kq/r² · r̂
double V = q.ElectricPotential(pos, fieldPoint);     // V = kq/r

// Coulomb force between two charges
double q2 = -2e-6;
Vector F = q.CoulombForce(q2, new Vector(0, 0, 0), new Vector(1, 0, 0));
```

**Superposition** — multiple charges:

```csharp
var charges = new[]
{
    (1e-6,  new Vector(-1, 0, 0)),
    (-1e-6, new Vector( 1, 0, 0))
};

Vector E = charges.ElectricFieldSuperposition(new Vector(0, 1, 0));
double V = charges.ElectricPotentialSuperposition(new Vector(0, 1, 0));
```

**VectorField** — create a field usable with ∇· and ∇×:

```csharp
// Single charge → VectorField
VectorField E = q.ElectricVectorField(pos);

// Verify Gauss's law: ρ = ε₀ ∇·E
double rho = E.ChargeDensity((2, 0, 0));

// Multiple charges → VectorField
VectorField Edipole = charges.ElectricVectorField();
```

### Lorentz Force

```csharp
double q = 1.6e-19;   // proton
var v = new Vector(1e6, 0, 0);
var E = new Vector(0, 100, 0);
var B = new Vector(0, 0, 0.5);

Vector F = q.LorentzForce(v, E, B);       // F = q(E + v × B)
Vector Fe = q.ElectricForce(E);            // F = qE
Vector Fm = q.MagneticForce(v, B);         // F = q(v × B)
```

### Maxwell's Equations (differential form)

All four Maxwell equations are accessible through `VectorField`:

```csharp
VectorField E = /* electric field */;
VectorField B = /* magnetic field */;
var point = (1.0, 2.0, 0.0);

// (1) Gauss (electric): ρ = ε₀ ∇·E
double rho = E.ChargeDensity(point);

// (2) Gauss (magnetic): ∇·B = 0  (verification)
double divB = B.GaussMagnetic(point);    // should be ≈ 0

// (3) Faraday: ∂B/∂t = −∇×E
Vector dBdt = E.FaradayLaw(point);

// (4) Ampère–Maxwell: J + ε₀∂E/∂t = (1/μ₀)∇×B
Vector Jeff = B.AmpereLaw(point);
```

**Wave equation** — verify a time-dependent scalar field satisfies ∇²f = μ₀ε₀ ∂²f/∂t²:

```csharp
// Plane wave: f(t)(r) = cos(kx − ωt)
double k = 2 * Math.PI;
double omega = k * PhysicsConstants.SpeedOfLight;

Func<double, Func<Vector, double>> wave = t =>
    r => Math.Cos(k * r.x - omega * t);

double residual = wave.WaveEquationResidual(t: 0, point: (0.5, 0, 0));
// residual ≈ 0 for a valid electromagnetic wave
```

### Energy & Momentum

```csharp
var E = new Vector(100, 0, 0);
var B = new Vector(0, 0, 0.001);

// Electromagnetic energy density: u = ½(ε₀|E|² + |B|²/μ₀)
double u = E.EnergyDensity(B);

// Poynting vector: S = (1/μ₀)(E × B)
Vector S = E.PoyntingVector(B);

// Radiation pressure (absorption or reflection)
double P_abs = E.RadiationPressure(B);
double P_ref = E.RadiationPressure(B, reflected: true);  // 2×
```

### Potentials

```csharp
// Electric field from scalar potential: E = −∇V
Func<Vector, double> V = r =>
    ElectroMagneticFieldExtensions.CoulombConstant * 1e-6 / r.GetMagnitude();

VectorField E = V.ElectricFieldFromPotential();
// E is now a full VectorField — use E.Divergence, E.Curl, etc.

// Magnetic field from vector potential: B = ∇×A
VectorField A = new VectorField(
    r => 0,
    r => 0,
    r => r.x * r.y);
Vector B = A.MagneticFieldFromVectorPotential((1, 1, 0));
```

### Biot–Savart Sources

Magnetic field from common current configurations:

```csharp
// Infinite straight wire: B = μ₀I/(2πd)
double B = 10.0.MagneticFieldFromWire(distance: 0.05);  // 10 A, 5 cm away

// Wire in 3D — full vector field
Vector Bvec = 10.0.MagneticFieldFromWire(
    wireDirection: new Vector(0, 0, 1),
    wirePoint:     new Vector(0, 0, 0),
    fieldPoint:    new Vector(0.05, 0, 0));

// Circular loop center: B = μ₀I/(2R) along normal
Vector Bloop = 5.0.MagneticFieldFromLoop(
    radius: 0.1,
    normal: new Vector(0, 0, 1));

// Solenoid: B = μ₀nI along axis
Vector Bsol = 2.0.MagneticFieldFromSolenoid(
    turnsPerUnitLength: 1000,
    direction: new Vector(0, 0, 1));
```

### Dipoles

```csharp
// Electric dipole moment: p = qd
var p = new Vector(0, 0, 1e-12); // 1 pC·m along Z
var dipolePos = new Vector(0, 0, 0);

double V = p.DipolePotential(dipolePos, fieldPoint: new Vector(1, 0, 0));
Vector E = p.DipoleElectricField(dipolePos, fieldPoint: new Vector(1, 0, 0));
// E = (1/4πε₀) · [3(p·r̂)r̂ − p] / r³
```

### Electrostatic Model (`IElectrostaticModel`)

The `IElectrostaticModel` interface provides an abstraction for electrostatic physics used by simulation engines. The default implementation `ElectrostaticModel` uses `PhysicsConstants.VacuumPermittivity`.

| Method | Description |
|--------|-------------|
| `VacuumPermittivity` | Returns $\varepsilon_0$ in F/m |
| `EffectivePermittivity(εᵣ)` | Computes $\varepsilon = \varepsilon_0 \cdot \varepsilon_r$ |
| `PoissonRhs(ρ, ε)` | Computes $-\rho / \varepsilon$ for the Poisson equation |

```csharp
using CSharpNumerics.Physics.Electromagnetism;
using CSharpNumerics.Physics.Electromagnetism.Interfaces;

IElectrostaticModel em = new ElectrostaticModel();

double eps  = em.EffectivePermittivity(relativePermittivity: 4.0);
double rhs  = em.PoissonRhs(chargeDensity: 1e-6, permittivity: eps);
```

---

## 🌿 Environmental Extensions

The `EnvironmentalExtensions` class bridges `ScalarField` (∇, ∇²) and `VectorField` to atmospheric and aquatic transport physics — Gaussian plume dispersion, Fickian diffusion, advection, and the advection–diffusion equation.

### Gaussian Plume

Steady-state atmospheric dispersion from a point source with ground reflection:

```csharp
using CSharpNumerics.Physics.Enums;

double Q = 10.0; // emission rate (kg/s)
var source = new Vector(0, 0, 50);     // 50 m stack
var wind   = new Vector(1, 0, 0);      // wind along +X

// Using Pasquill–Gifford stability class
ScalarField C = Q.GaussianPlume(
    windSpeed: 5.0,
    stackHeight: 50,
    sourcePosition: source,
    windDirection: wind,
    stability: StabilityClass.D);        // neutral conditions

double conc = C.Evaluate(new Vector(500, 0, 0));  // concentration 500 m downwind
Vector grad = C.Gradient((500, 0, 0));             // concentration gradient
double lap  = C.Laplacian((500, 0, 0));            // ∇²C

// Quick ground-level centerline concentration
double Cgl = Q.GaussianPlumeGroundLevel(
    windSpeed: 5.0, stackHeight: 50,
    downwindDistance: 1000,
    stability: StabilityClass.C);
```

**Custom dispersion** — supply your own σy(x), σz(x):

```csharp
ScalarField C = Q.GaussianPlume(
    windSpeed: 5.0,
    stackHeight: 50,
    sourcePosition: source,
    windDirection: wind,
    sigmaY: x => 0.22 * x / Math.Sqrt(1 + 0.0001 * x),
    sigmaZ: x => 0.20 * x);
```

**Briggs dispersion parameters** are also available directly:

```csharp
double sy = EnvironmentalExtensions.BriggsSigmaY(1000, StabilityClass.D);
double sz = EnvironmentalExtensions.BriggsSigmaZ(1000, StabilityClass.D);
```

### Transient Gaussian Puff

Time-dependent dispersion of an instantaneous release that advects downwind — the puff centre moves at wind speed while the cloud expands according to Briggs σ(u·t):

```csharp
double Q = 5.0;  // emission rate (kg/s)
double releaseSeconds = 10.0;  // release duration → mass = Q·Δt

ScalarField C = Q.GaussianPuff(
    releaseSeconds: releaseSeconds,
    windSpeed: 10.0,
    stackHeight: 50,
    sourcePosition: new Vector(0, 0, 50),
    windDirection: new Vector(1, 0, 0),
    time: 30.0,                          // seconds since release
    stability: StabilityClass.D);

// At t=30s, puff centre is at x = u·t = 300 m downwind
double conc = C.Evaluate(new Vector(300, 0, 50));  // peak region

// Evaluate at different times for animation
for (double t = 10; t <= 120; t += 10)
{
    var puff = Q.GaussianPuff(releaseSeconds, 10, 50,
        new Vector(0, 0, 50), new Vector(1, 0, 0), t, StabilityClass.D);
    double peak = puff.Evaluate(new Vector(10 * t, 0, 50));
}
```

### Diffusion (Fick's Laws)

```csharp
var C = new ScalarField(r => Math.Exp(-(r.x * r.x + r.y * r.y)));
double D = 0.01; // m²/s

// Fick's first law: diffusive flux J = −D∇C → VectorField
VectorField J = C.DiffusionFlux(D);

// Fick's second law: ∂C/∂t = D∇²C → ScalarField
ScalarField dCdt = C.DiffusionRate(D);
double rate = dCdt.Evaluate((0.5, 0, 0));
```

**Analytical point-source diffusion** — 3D Gaussian spreading:

```csharp
// Mass M released at origin, evaluated at time t
ScalarField C = 1.0.DiffusionPointSource(
    diffusionCoefficient: 0.01,
    time: 100,
    sourcePosition: new Vector(0, 0, 0));
// C(r,t) = M / (4πDt)^(3/2) · exp(−|r|² / (4Dt))
```

### Advection

```csharp
var C = new ScalarField(r => Math.Exp(-r.x * r.x));
var wind = new VectorField(r => 2.0, r => 0, r => 0);  // uniform 2 m/s along X

// ∂C/∂t = −v·∇C
ScalarField dCdt = C.AdvectionRate(wind);
```

### Advection–Diffusion

Combined transport with optional source term:

```csharp
var C = new ScalarField(r => Math.Exp(-(r.x * r.x + r.y * r.y)));
var wind = new VectorField(r => 1.0, r => 0.5, r => 0);
double D = 0.1;

// ∂C/∂t = D∇²C − v·∇C
ScalarField dCdt = C.AdvectionDiffusionRate(wind, D);

// With source term S(r)
var source = new ScalarField(r => r.x > 0 && r.x < 1 ? 0.5 : 0);
ScalarField dCdt2 = C.AdvectionDiffusionRate(wind, D, source);
```

**Péclet number** — determines transport regime:

```csharp
double Pe = 2.0.PecletNumber(characteristicLength: 100, diffusionCoefficient: 0.1);
// Pe = 2000 → advection-dominated
```

### Rothermel Surface Fire Spread

The `RothermelModel` static class implements the **Rothermel (1972) surface fire spread model** — the standard used by FARSITE, FlamMap, and BehavePlus. All parameters are in **SI units**.

**Rate of spread:**

$$R = \frac{I_R \cdot \xi \cdot (1 + \phi_w + \phi_s)}{\rho_b \cdot \varepsilon \cdot Q_{ig}}$$

```csharp
using CSharpNumerics.Physics.Environmental.Fire;
using CSharpNumerics.Physics.Materials.Fire;
using CSharpNumerics.Physics.Materials.Fire.Enums;

// Get a standard Anderson 13 fuel model
FuelModel grass = FuelLibrary.Get(FuelModelType.ShortGrass);

// Rate of spread on flat terrain, moderate wind
double ros = RothermelModel.RateOfSpread(
    fuel: grass,
    moistureContent: 0.05,     // 5% dead fuel moisture
    windSpeed: 2.24,           // m/s (≈ 5 mph mid-flame)
    slopeRadians: 0);          // flat
// ros ≈ 23–25 m/min for Short Grass

// Individual sub-equations
double IR   = RothermelModel.ReactionIntensity(grass, moistureContent: 0.05);
double phiW = RothermelModel.WindFactor(grass, midflameWindSpeed: 2.24);
double phiS = RothermelModel.SlopeFactor(
    RothermelModel.PackingRatio(grass), slopeRadians: Math.PI / 6);
double xi   = RothermelModel.PropagatingFluxRatio(grass);
double Qig  = RothermelModel.HeatOfPreignition(moistureContent: 0.05);
double eps  = RothermelModel.EffectiveHeatingNumber(grass);

// Flame length from Byram's fireline intensity correlation
double fl = RothermelModel.FlameLength(IR, ros);  // metres
```

**Sub-equation reference:**

| Method | Symbol | Unit |
|--------|--------|------|
| `ReactionIntensity` | $I_R$ | kJ/m²·min |
| `WindFactor` | $\phi_w$ | — |
| `SlopeFactor` | $\phi_s$ | — |
| `PropagatingFluxRatio` | $\xi$ | — |
| `HeatOfPreignition` | $Q_{ig}$ | kJ/kg |
| `EffectiveHeatingNumber` | $\varepsilon$ | — |
| `PackingRatio` | $\beta$ | — |
| `OptimalPackingRatio` | $\beta_{op}$ | — |
| `FlameLength` | $L$ | m |

### Fuel Models (Anderson 13)

The `FuelModel` readonly struct holds Rothermel fuel parameters. All 13 standard Anderson models are pre-loaded in `FuelLibrary`:

```csharp
// All models are in SI units
FuelModel chaparral = FuelLibrary.Get(FuelModelType.Chaparral);
// chaparral.SurfaceAreaToVolumeRatio = 4921 (1/m)
// chaparral.FuelBedDepth             = 1.829 (m)
// chaparral.OvendryFuelLoad          = 3.663 (kg/m²)
// chaparral.MoistureOfExtinction     = 0.20

// Iterate all models
foreach (var fuel in FuelLibrary.All)
    Console.WriteLine($"{fuel.Type}: δ={fuel.FuelBedDepth:F3} m");

// Register a custom fuel model at runtime
FuelLibrary.Register(new FuelModel(
    (FuelModelType)99, "Custom Sage", sigma: 5000,
    delta: 0.5, w0: 1.2, Mx: 0.25, h: 18000));
```

| Model | Type | σ (1/m) | δ (m) | w₀ (kg/m²) | Mₓ |
|-------|------|---------|-------|-------------|-----|
| 1 | Short Grass | 11,483 | 0.305 | 0.166 | 0.12 |
| 2 | Timber/Grass Understory | 3,281 | 0.305 | 0.897 | 0.15 |
| 3 | Tall Grass | 4,921 | 0.762 | 0.675 | 0.25 |
| 4 | Chaparral | 6,562 | 1.829 | 3.663 | 0.20 |
| 5 | Brush | 6,562 | 0.610 | 0.784 | 0.20 |
| 6 | Dormant Brush | 5,741 | 0.762 | 0.672 | 0.25 |
| 7 | Southern Rough | 5,741 | 0.762 | 0.529 | 0.40 |
| 8 | Closed Timber Litter | 6,562 | 0.061 | 1.121 | 0.30 |
| 9 | Hardwood Litter | 8,202 | 0.061 | 0.327 | 0.25 |
| 10 | Timber/Litter Understory | 6,562 | 0.305 | 1.121 | 0.25 |
| 11 | Light Logging Slash | 4,921 | 0.305 | 1.345 | 0.15 |
| 12 | Medium Logging Slash | 4,921 | 0.701 | 3.363 | 0.20 |
| 13 | Heavy Logging Slash | 4,921 | 0.914 | 5.604 | 0.25 |

### Water Hydraulics & Contaminant Data

The `Physics.Environmental.Water` namespace provides open-channel hydraulics (Manning's equation) and longitudinal dispersion for river transport modelling. The `Physics.Materials.Water` namespace defines aquatic contaminant descriptors with decay, adsorption, and toxicity properties.

#### Manning's Equation — Open-Channel Velocity

$$u = \frac{1}{n} R_h^{2/3} S_0^{1/2}$$

```csharp
using CSharpNumerics.Physics.Environmental.Water;

// Rectangular channel: W=10 m, H=2 m, Manning's n=0.035, slope=0.001
double Rh = ManningEquation.RectangularHydraulicRadius(10, 2);
double u  = ManningEquation.Velocity(0.035, Rh, 0.001);     // ≈ 1.14 m/s
double Q  = ManningEquation.Discharge(0.035, Rh, 0.001, 20); // 20 m² area

// Trapezoidal channel
double RhTrap = ManningEquation.TrapezoidalHydraulicRadius(8, 2, 1.5);
double uTrap  = ManningEquation.Velocity(0.03, RhTrap, 0.0005);
```

#### Fischer Longitudinal Dispersion

$$E_L = 0.011 \frac{u^2 W^2}{H \cdot u_*}$$

```csharp
double uStar = LongitudinalDispersion.ShearVelocity(Rh, 0.001);
double EL = LongitudinalDispersion.FischerCoefficient(u, 10, 2, uStar);

// Decay & retardation helpers
double lambda = LongitudinalDispersion.DecayConstant(halfLifeSeconds: 9.5e8); // Cs-137
double Rf = LongitudinalDispersion.RetardationFactor(
    bulkDensity: 1600, Kd: 50, porosity: 0.4); // Rf > 1 → retarded transport
```

#### Tributary Mixing

```csharp
// Two rivers merge: main at 80 mg/L (Q=20 m³/s) + tributary at 0 mg/L (Q=10 m³/s)
double downstreamConc = MixingZoneModel.TributaryMixing(80, 20, 0, 10);
// → 53.3 mg/L (mass-balance dilution)

double mixLength = MixingZoneModel.MixingLength(1.0, 10, 0.6); // ≈ 116.7 m
```

#### Aquatic Contaminants

Pre-loaded contaminant library with decay, adsorption, and toxicity data:

```csharp
using CSharpNumerics.Physics.Materials.Water;

// Built-in contaminants: Cs137, Sr90, I131, Benzene, Toluene, Cyanide,
// Mercury, Lead, Arsenic, EColi, Enterococcus, GenericHeat
var cs137 = ContaminantLibrary.Get("Cs-137");
// cs137.HalfLifeSeconds ≈ 9.5×10⁸, cs137.PartitionCoefficient = 1000, etc.

var benzene = AquaticContaminant.Benzene;
bool isConservative = benzene.IsConservative;  // false (has half-life)
double lambda = benzene.DecayConstant;          // ln(2) / t½

// Create a custom contaminant
var custom = new AquaticContaminant("PFOS",
    ContaminantType.Chemical,
    halfLifeSeconds: 0,       // persistent (conservative tracer)
    partitionCoefficient: 10,
    toxicityThresholdMgL: 0.00007,
    lethalThresholdMgL: 50);

ContaminantLibrary.Register(custom);
```

| Contaminant | Type | Half-life | Kd (L/kg) | Toxicity (mg/L) |
|-------------|------|-----------|-----------|-----------------|
| Cs-137 | Radioactive | 30.17 yr | 1000 | 0.002 |
| Sr-90 | Radioactive | 28.8 yr | 200 | 0.008 |
| I-131 | Radioactive | 8.02 d | 10 | 0.003 |
| Benzene | Chemical | 180 d | 1.8 | 0.005 |
| Cyanide | Chemical | ∞ | 0 | 0.07 |
| Mercury | Chemical | ∞ | 5000 | 0.001 |
| E. coli | Biological | 2 d | 0 | 0.0001 |

---

## 🔥 Heat Extensions

The `HeatExtensions` class provides extension methods for heat transfer calculations — **conduction, convection, radiation, the heat equation, dimensionless numbers, lumped-capacitance transient analysis, and fins**.

### Conduction (Fourier's Law)

```csharp
// 1-D heat flux through a slab: q = k·ΔT / L
double q = 50.0.ConductiveHeatFlux(temperatureDifference: 80, thickness: 0.2);  // W/m²

// Heat rate through a slab of area A: Q̇ = k·A·ΔT / L
double Q = 50.0.ConductiveHeatRate(area: 2.0, temperatureDifference: 80, thickness: 0.2);

// 3-D heat flux vector field: q = −k∇T
var T = new ScalarField(r => 300 + 100 * Math.Exp(-(r.x * r.x + r.y * r.y)));
VectorField heatFlux = T.ConductiveHeatFlux(thermalConductivity: 50);
```

**Thermal resistance** for common geometries:

```csharp
// Flat slab: R = L / (kA)
double Rslab = 0.1.SlabThermalResistance(thermalConductivity: 50, area: 2.0);

// Cylindrical shell: R = ln(r₂/r₁) / (2πkL)
double Rcyl = 0.05.CylindricalThermalResistance(outerRadius: 0.10, thermalConductivity: 16, length: 1.0);

// Spherical shell: R = (1/r₁ − 1/r₂) / (4πk)
double Rsph = 0.05.SphericalThermalResistance(outerRadius: 0.10, thermalConductivity: 16);

// Series and parallel resistance networks
double Rseries  = new[] { 0.01, 0.05, 0.02 }.SeriesThermalResistance();
double Rparallel = new[] { 0.05, 0.10 }.ParallelThermalResistance();
```

### Convection (Newton's Law of Cooling)

```csharp
// Heat flux: q = h·(T_s − T_∞)
double q = 25.0.ConvectiveHeatFlux(surfaceTemperature: 400, fluidTemperature: 300);

// Heat rate: Q̇ = h·A·(T_s − T_∞)
double Q = 25.0.ConvectiveHeatRate(area: 0.5, surfaceTemperature: 400, fluidTemperature: 300);

// Convective thermal resistance: R = 1/(hA)
double Rconv = 25.0.ConvectiveThermalResistance(area: 0.5);
```

### Radiation (Stefan–Boltzmann)

```csharp
// Emissive power: q = εσT⁴
double q = 500.0.RadiativeHeatFlux(emissivity: 0.8);

// Net radiative flux: q_net = εσ(T_s⁴ − T_sur⁴)
double qNet = 500.0.NetRadiativeHeatFlux(surroundingTemperature: 300, emissivity: 0.8);

// Radiative heat rate: Q̇ = εσA(T_s⁴ − T_sur⁴)
double Q = 500.0.RadiativeHeatRate(surroundingTemperature: 300, area: 2.0, emissivity: 0.8);

// Linearised radiative thermal resistance
double Rrad = 500.0.RadiativeThermalResistance(surroundingTemperature: 300, area: 2.0, emissivity: 0.8);
```

### Heat Equation

```csharp
// Thermal diffusivity: α = k / (ρ·c_p)
double alpha = 50.0.ThermalDiffusivity(density: 7800, specificHeat: 500);

// ∂T/∂t = α∇²T
var T = new ScalarField(r => 300 + 200 * Math.Exp(-(r.x * r.x)));
ScalarField dTdt = T.HeatEquationRate(thermalDiffusivity: alpha);

// With volumetric source: ∂T/∂t = α∇²T + q̇/(ρc_p)
var source = new ScalarField(r => 1e6);  // 1 MW/m³
ScalarField dTdt2 = T.HeatEquationRate(alpha, source, density: 7800, specificHeat: 500);

// Semi-infinite solid analytical solution: T(x,t) = T_i + (T_s − T_i)·erfc(x/(2√(αt)))
double Tx = 20.0.SemiInfiniteTemperature(surfaceTemperature: 100, thermalDiffusivity: alpha,
    depth: 0.05, time: 60);
```

### Dimensionless Numbers

```csharp
// Biot number: Bi = hL_c / k   (Bi < 0.1 → lumped model valid)
double Bi = 25.0.BiotNumber(characteristicLength: 0.01, thermalConductivity: 50);

// Nusselt number: Nu = hL / k_f
double Nu = 25.0.NusseltNumber(characteristicLength: 0.5, fluidConductivity: 0.6);

// Fourier number: Fo = αt / L²
double Fo = 1.28e-5.FourierNumber(time: 300, characteristicLength: 0.01);

// Prandtl number: Pr = c_p·μ / k
double Pr = 4184.0.PrandtlNumber(dynamicViscosity: 8.9e-4, thermalConductivity: 0.6);

// Grashof number: Gr = gβΔTL³ / ν²
double Gr = 50.0.GrashofNumber(characteristicLength: 0.3,
    thermalExpansionCoefficient: 3.4e-3, kinematicViscosity: 1.5e-5);

// Rayleigh number: Ra = Gr·Pr = gβΔTL³ / (να)
double Ra = 50.0.RayleighNumber(characteristicLength: 0.3,
    thermalExpansionCoefficient: 3.4e-3, kinematicViscosity: 1.5e-5, thermalDiffusivity: 2.2e-5);
```

### Lumped Capacitance

```csharp
// Time constant: τ = ρVc_p / (hA)
double tau = 7800.0.ThermalTimeConstant(volume: 1e-4, specificHeat: 500,
    heatTransferCoefficient: 25, surfaceArea: 0.02);

// Temperature at time t: T(t) = T_∞ + (T₀ − T_∞)·exp(−t/τ)
double T = 500.0.LumpedCapacitanceTemperature(ambientTemperature: 300, time: 60, timeConstant: tau);

// Time to reach target temperature
double t = 500.0.LumpedCapacitanceTime(targetTemperature: 350, ambientTemperature: 300, timeConstant: tau);
```

### Fins

```csharp
// Fin parameter: m = √(hP / (kA_c))
double m = 25.0.FinParameter(perimeter: 0.1, thermalConductivity: 200, crossSectionalArea: 5e-4);

// Infinite fin heat rate: Q̇ = √(hPkA_c)·(T_b − T_∞)
double Qinf = 400.0.InfiniteFinHeatRate(ambientTemperature: 300,
    heatTransferCoefficient: 25, perimeter: 0.1, thermalConductivity: 200, crossSectionalArea: 5e-4);

// Finite fin with insulated tip: Q̇ = √(hPkA_c)·(T_b − T_∞)·tanh(mL)
double Qfin = 400.0.InsulatedTipFinHeatRate(ambientTemperature: 300,
    heatTransferCoefficient: 25, perimeter: 0.1, thermalConductivity: 200,
    crossSectionalArea: 5e-4, finLength: 0.05);

// Fin efficiency: η = tanh(mL) / (mL)
double eta = m.InsulatedTipFinEfficiency(finLength: 0.05);
```

### Heat Transfer Model (`IHeatTransferModel`)

The `IHeatTransferModel` interface provides an abstraction for thermal physics used by simulation engines. The default implementation `HeatTransferModel` delegates to `HeatExtensions`.

| Method | Description |
|--------|-------------|
| `ThermalDiffusivity(k, ρ, cₚ)` | Computes $\alpha = k / (\rho c_p)$ |
| `HeatSourceRate(Q, ρ, cₚ)` | Converts volumetric power to temperature rate $Q / (\rho c_p)$ |
| `ConvectiveBoundaryRate(h, ρ, cₚ, dx, T, T∞)` | Robin BC correction: $-h / (\rho c_p \cdot dx) \cdot (T - T_\infty)$ |

```csharp
using CSharpNumerics.Physics.Thermodynamics;
using CSharpNumerics.Physics.Thermodynamics.Interfaces;

IHeatTransferModel heat = new HeatTransferModel();

double alpha = heat.ThermalDiffusivity(conductivity: 50, density: 7800, specificHeat: 500);
// α ≈ 1.28e-5 m²/s (steel)

double rate = heat.HeatSourceRate(power: 1e6, density: 7800, specificHeat: 500);
// temperature rate from 1 MW/m³ source

// Robin / convection boundary correction (used by HeatPlate / HeatBlock3D solvers)
double correction = heat.ConvectiveBoundaryRate(
    h: 50, density: 2700, specificHeat: 900, dx: 0.01,
    cellTemperature: 100, ambientTemperature: 20);
// < 0 — boundary cools towards ambient
```

---

## 💧 Fluid Extensions

The `FluidExtensions` class bridges `VectorField` (∇·, ∇×) and `ScalarField` (∇, ∇²) to classical fluid dynamics — Navier–Stokes equations, Bernoulli's principle, continuity, vorticity, drag/lift, dimensionless numbers, viscous pipe flow, and hydrostatics.

### Navier–Stokes Equations (Incompressible)

$$\rho\!\left(\frac{\partial \mathbf{v}}{\partial t} + (\mathbf{v}\cdot\nabla)\mathbf{v}\right) = -\nabla p + \mu\nabla^2\mathbf{v} + \mathbf{f}, \qquad \nabla\cdot\mathbf{v} = 0$$

Individual terms and the full residual:

```csharp
var velocity = new VectorField(
    r => -r.y, r => r.x, r => 0);  // rigid-body rotation
var pressure = new ScalarField(
    r => 0.5 * 1000 * (r.x * r.x + r.y * r.y));
var point = (1.0, 0.0, 0.0);

// Convective acceleration: (v·∇)v
Vector conv = velocity.ConvectiveAcceleration(point);

// Viscous term: μ∇²v
Vector visc = velocity.ViscousTerm(dynamicViscosity: 1e-3, point);

// Pressure gradient force: −∇p
Vector gradP = pressure.PressureGradientForce(point);

// Full residual: ∂v/∂t = (−∇p + μ∇²v + f)/ρ − (v·∇)v
// For steady-state flow this should be ≈ 0
Vector residual = velocity.NavierStokesResidual(
    pressure, density: 1000, dynamicViscosity: 1e-3, point);

// With body force (e.g. gravity)
Vector residualG = velocity.NavierStokesResidual(
    pressure, density: 1000, dynamicViscosity: 1e-3, point,
    bodyForce: new Vector(0, 0, -9810));

// Incompressibility check: ∇·v = 0
double divV = velocity.IncompressibilityResidual(point);  // should be ≈ 0
```

### Euler Equations (Inviscid Flow)

Navier–Stokes with $\mu = 0$:

```csharp
Vector euler = velocity.EulerEquationResidual(
    pressure, density: 1000, point);

// With body force
Vector eulerG = velocity.EulerEquationResidual(
    pressure, density: 1000, point,
    bodyForce: new Vector(0, 0, -9810));
```

### Bernoulli's Principle

$$p + \tfrac{1}{2}\rho v^2 + \rho g h = \text{const}$$

Along a streamline for steady, inviscid, incompressible flow:

```csharp
// Bernoulli constant
double B = 101325.0.BernoulliConstant(density: 1000, velocity: 5.0, height: 10);

// Pressure at a second point
double p2 = 101325.0.BernoulliPressure(
    density: 1000, v1: 2.0, v2: 8.0, h1: 0, h2: 3.0);

// Dynamic pressure: q = ½ρv²
double q = 1.225.DynamicPressure(speed: 50);

// Stagnation (total) pressure: p₀ = p + ½ρv²
double p0 = 101325.0.StagnationPressure(density: 1.225, speed: 50);
```

### Continuity Equation

$$\frac{\partial\rho}{\partial t} + \nabla\cdot(\rho\mathbf{v}) = 0$$

```csharp
// Mass flux as a vector field: ρv (constant density)
VectorField massFlux = velocity.MassFlux(density: 1000);

// Mass flux with spatially varying density
var rhoField = new ScalarField(r => 1000 + 10 * r.x);
VectorField massFlux2 = velocity.MassFlux(rhoField);

// Continuity residual: ∇·(ρv) — should be 0 for incompressible flow
double cont = velocity.ContinuityResidual(density: 1000, point);

// Volume flow rate: Q = A·v
double Q = 0.01.VolumeFlowRate(speed: 3.0);   // 0.03 m³/s

// Speed from continuity: A₁v₁ = A₂v₂
double v2 = 0.01.ContinuitySpeed(speed1: 3.0, area2: 0.005);  // 6.0 m/s
```

### Vorticity & Circulation

$$\boldsymbol{\omega} = \nabla\times\mathbf{v}$$

```csharp
// Vorticity: ω = ∇×v
Vector omega = velocity.Vorticity(point);

// Enstrophy density: ½|ω|²
double enstrophy = velocity.Enstrophy(point);

// Helicity density: v·ω
double helicity = velocity.HelicityDensity(point);
```

### Stream Function (2-D)

For incompressible 2-D flow, $v_x = \partial\psi/\partial y$, $v_y = -\partial\psi/\partial x$:

```csharp
// Create velocity field from a stream function
var psi = new ScalarField(r => r.x * r.y);  // ψ = xy → saddle flow
VectorField v = psi.VelocityFromStreamFunction();
// vx = ∂ψ/∂y = x,  vy = −∂ψ/∂x = −y → strain flow

// The resulting field is automatically divergence-free
double div = v.Divergence((1, 1, 0));  // ≈ 0
```

### Drag & Lift

```csharp
// Drag force: F_D = ½ρv²C_D A
double drag = 0.47.DragForce(density: 1.225, speed: 30, referenceArea: 0.01);

// Lift force: F_L = ½ρv²C_L A
double lift = 1.2.LiftForce(density: 1.225, speed: 60, referenceArea: 20);

// Terminal velocity: v_t = √(2mg / (ρC_D A))
double vTerm = 0.5.TerminalVelocity(
    dragCoefficient: 0.47, density: 1.225, referenceArea: 0.01);
```

### Dimensionless Numbers

```csharp
// Reynolds number: Re = ρvL/μ
double Re = 1000.0.ReynoldsNumber(speed: 2.0, characteristicLength: 0.1,
    dynamicViscosity: 1e-3);  // Re = 200 000

// Mach number: Ma = v/c (defaults to speed of sound in air)
double Ma = 300.0.MachNumber();           // Ma ≈ 0.875
double Ma2 = 1500.0.MachNumber(1480);    // in water

// Froude number: Fr = v/√(gL)
double Fr = 5.0.FroudeNumber(characteristicLength: 2.0);

// Strouhal number: St = fL/v
double St = 50.0.StrouhalNumber(characteristicLength: 0.1, speed: 25);

// Weber number: We = ρv²L/σ
double We = 1000.0.WeberNumber(speed: 1.0, characteristicLength: 0.001,
    surfaceTension: 0.072);
```

### Viscous / Pipe Flow

**Hagen–Poiseuille** flow in circular pipes:

```csharp
// Volume flow rate: Q = πR⁴ΔP / (8μL)
double Q = 0.005.PoiseuilleFlowRate(
    pressureDrop: 1000, dynamicViscosity: 1e-3, length: 1.0);

// Velocity profile: v(r) = (ΔP/4μL)(R² − r²)
double vCenter = 0.005.PoiseuilleVelocity(

---

## 🔦 Optics

The `CSharpNumerics.Physics.Optics` namespace provides geometric optics: rays, reflection/refraction (Snell's law, Fresnel equations), optical elements (mirrors, lenses, prisms), apertures, sensors, and a recursive ray tracer.

### Core Types

| Type | Description |
|------|-------------|
| `Ray` | Light ray: origin, direction, wavelength (nm), intensity |
| `OpticalMedium` | Refractive index, absorption coefficient, Abbe number for dispersion |
| `RayHit` | Intersection result: point, normal, distance, medium |
| `IOpticalSurface` | Interface for any surface interacting with rays |

```csharp
var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1), wavelengthNm: 550);
var glass = OpticalMaterialLibrary.CrownGlass;  // n=1.5168, V=64.17
double nBlue = glass.RefractiveIndexAt(450);     // Cauchy dispersion
```

### Reflection & Refraction

```csharp
// Snell's law
double? theta2 = OpticsExtensions.RefractionAngle(n1: 1.0, n2: 1.5, thetaIncident: Math.PI / 6);

// Total internal reflection
double critical = OpticsExtensions.CriticalAngle(n1: 1.5, n2: 1.0);
bool tir = OpticsExtensions.IsTotalInternalReflection(1.5, 1.0, thetaIncident);

// Fresnel reflectance (unpolarised)
double R = OpticsExtensions.FresnelReflectance(1.0, 1.5, thetaIncident);

// Schlick approximation (fast, game-friendly)
double Rs = OpticsExtensions.SchlickReflectance(1.0, 1.5, cosTheta);

// Beer–Lambert absorption
double transmittance = OpticsExtensions.BeerLambertTransmittance(alpha: 0.5, distance: 2.0);

// Vector-level reflection and refraction
Vector reflected = OpticsExtensions.Reflect(incident, normal);
Vector? refracted = OpticsExtensions.Refract(incident, normal, n1, n2);
```

### Mirrors

```csharp
// Plane mirror
var mirror = new PlaneMirror(center: new Vector(0, 0, 10), normal: new Vector(0, 0, -1));
RayHit? hit = mirror.Intersect(ray);

// Spherical mirror (concave/convex)
var concave = new SphericalMirror(
    centreOfCurvature: new Vector(0, 0, 20),
    radiusOfCurvature: 20, MirrorType.Concave);
double imageDistance = concave.ImageDistance(objectDistance: 30); // 1/f = 1/do + 1/di
double magnification = concave.Magnification(objectDistance: 30);
```

### Thin Lenses

```csharp
var lens = new ThinLens(
    center: new Vector(0, 0, 5),
    axis: new Vector(0, 0, 1),
    focalLength: 10.0, LensType.Converging);

double di = lens.ImageDistance(objectDistance: 20);   // thin-lens equation
double m  = lens.Magnification(objectDistance: 20);

// Lensmaker's equation: 1/f = (n-1)(1/R₁ - 1/R₂)
double f = ThinLens.LensmakerFocalLength(n: 1.5168, r1: 0.2, r2: -0.2);

// Refract a ray through the lens (paraxial model)
Ray outRay = lens.RefractRay(incomingRay);
```

### Prisms & Dispersion

```csharp
var prism = new Prism(apexAngleRadians: Math.PI / 3, OpticalMaterialLibrary.CrownGlass);

double? deviation = prism.DeviationAngle(thetaIncident: 0.5, wavelengthNm: 550);
double dMin = prism.MinimumDeviation(wavelengthNm: 589.3);

// Angular dispersion between blue and red
double spread = prism.AngularDispersion(lambda1Nm: 450, lambda2Nm: 650, thetaIncident: 0.5);
```

### Apertures & Sensors

```csharp
// Circular aperture (iris)
var aperture = new CircularAperture(center, normal, radius: 2.0);

// Rectangular aperture (slit)
var slit = new RectangularAperture(center, normal, right, width: 4.0, height: 1.0);

// Image sensor (CCD / film plane)
var sensor = new ImageSensor(center, normal, right, width: 4.0, height: 4.0,
    resolutionX: 256, resolutionY: 256);
// After tracing: sensor.HitCount, sensor.GetPixelIntensity(x, y), sensor.GetImage()
```

### Ray Tracer

```csharp
var scene = new OpticalScene();
scene.Add(new ThinLens(lensCenter, axis, focalLength: 10.0));
scene.Add(sensor);

var tracer = new RayTracer(scene) { MaxBounces = 16, IntensityCutoff = 0.001 };
TraceResult result = tracer.Trace(ray);

foreach (var segment in result.Segments)
    Console.WriteLine($"{segment.Start} → {segment.End}  I={segment.Intensity:F3}");
```

### Optical Materials

Pre-defined media are available via `OpticalMaterialLibrary` or the `Materials.Optical()` factory:

```csharp
var diamond = OpticalMaterialLibrary.Diamond;       // n=2.417
var bk7 = Materials.Optical("CrownGlass");          // factory lookup
var flint = Materials.Optical("SF11");               // n=1.7847, low Abbe → high dispersion
```

### Game Engine Integration

`RaycastExtensions` casts optics rays against game-engine bounding volumes:

```csharp
using CSharpNumerics.Engines.Game;

var ray = new Ray(origin, direction);
double? tAABB   = ray.IntersectAABB(aabb);
double? tSphere = ray.IntersectSphere(sphere);
var detailed    = ray.IntersectSphereDetailed(sphere); // (point, normal)
```
    pressureDrop: 1000, dynamicViscosity: 1e-3, length: 1.0,
    radialDistance: 0);    // maximum at centre

double vWall = 0.005.PoiseuilleVelocity(
    pressureDrop: 1000, dynamicViscosity: 1e-3, length: 1.0,
    radialDistance: 0.005);  // zero at wall
```

**Stokes drag** on a sphere in creeping flow:

```csharp
// F = 6πμRv
double F = 1e-3.StokesDrag(radius: 0.01, speed: 0.1);
```

### Hydrostatics

```csharp
// Hydrostatic pressure at depth: p = p₀ + ρgh
double p = 10.0.HydrostaticPressure(density: 1025);  // 10 m depth in seawater

// Buoyant force (Archimedes): F_b = ρ_fluid · V · g
double Fb = 1025.0.BuoyantForce(displacedVolume: 0.5);  // ≈ 5025 N
```

### Fluid Energy & Momentum

```csharp
var v = new Vector(3, 4, 0);

// Kinetic energy density: e_k = ½ρ|v|²
double ek = v.KineticEnergyDensity(density: 1000);  // 12 500 J/m³

// Momentum density: ρv
Vector mom = v.MomentumDensity(density: 1000);  // (3000, 4000, 0) kg/(m²·s)
```

---

### ☢️ Nuclear Physics & Radioactive Fallout

Model radioactive isotopes, decay chains, radiation dose, and couple them to the dispersion pipeline for fallout scenarios that produce activity (Bq) and dose (Sv) maps alongside concentration.

#### Isotopes & Library

```csharp
using CSharpNumerics.Physics.Materials.Nuclear.Isotopes;

// Built-in isotopes
Isotope cs = Isotope.Cs137;           // Cs-137: t½ = 30.17 years
Isotope iodine = Isotope.I131;        // I-131: t½ = 8.02 days

// Lookup by name
Isotope sr = IsotopeLibrary.Get("Sr90");
bool found = IsotopeLibrary.TryGet("Co60", out Isotope co);

// Properties
double halfLife = cs.HalfLifeSeconds;       // ~9.51e8 s
double lambda = cs.Lambda;                   // decay constant (ln2/t½)
double activity = cs.SpecificActivityBqKg;   // Bq per kg
bool stable = Isotope.Ba137.IsStable;        // true

// Filter by element
var caesiumIsotopes = IsotopeLibrary.ByElement(55);  // Z = 55

// Register custom isotope at runtime
IsotopeLibrary.Register(new Isotope("Pu239", 94, 239, 7.6e11, 2.3e9,
    DecayMode.Alpha, 0.0, 4.4e-5, 0));
```

#### Radioactive Decay

```csharp
using CSharpNumerics.Physics.Materials.Nuclear.Decay;
using CSharpNumerics.Physics.Materials.Nuclear.DecayChains;

// Simple decay
double a0 = Decay.Activity(massKg: 0.001, Isotope.Cs137);      // initial Bq
double aT = Decay.ActivityAtTime(a0, timeSeconds: 3600, Isotope.Cs137);
double remaining = Decay.RemainingMass(0.001, 3600, Isotope.Cs137);

// Decay chain (Bateman equations)
var chain = DecayChain.Cs137Chain();    // Cs137 → Ba137m → Ba137
double[] masses = chain.Evolve(
    initialMasses: new[] { 1.0, 0.0, 0.0 },
    timeSeconds: 3600);

double[] activities = chain.EvolveActivity(
    new[] { 1.0, 0.0, 0.0 }, 3600);  // Bq for each isotope

// I-131 chain
var iChain = DecayChain.I131Chain();   // I131 → Xe131
```

#### Radiation Dose

```csharp
using CSharpNumerics.Physics.Materials.Nuclear.RadiationDose;

// Point-source dose rate (inverse-square law)
double svPerHour = RadiationDose.DoseRate(
    activityBq: 1e9, distanceM: 1.0, Isotope.Cs137);

// Ground-shine dose (deposition on surface)
double sv = RadiationDose.GroundShineDose(
    depositionBqM2: 1e6, Isotope.Cs137, timeSeconds: 3600);

// Inhalation dose
double inhDose = RadiationDose.InhalationDose(
    concentrationBqM3: 100,
    breathingRateM3s: RadiationDose.DefaultBreathingRate,
    exposureTimeSeconds: 3600,
    Isotope.Cs137);
```

---

### 🧪 Chemical Materials — Toxic & Flammable Gases

**Namespace:** `CSharpNumerics.Physics.Materials.Chemical`

Model hazardous chemical substances with toxicological thresholds (IDLH, ERPG, LC50, TLV) and unit conversions (kg/m³ ↔ ppm). Integrates with the GIS dispersion pipeline via `.WithMaterial(Materials.Chemical("Cl2"))`.

#### Built-in substances

| Formula | Name | MW (g/mol) | IDLH (ppm) | ERPG-2 (ppm) | ERPG-3 (ppm) | Phase |
|---------|------|-----------|------------|---------------|---------------|-------|
| Cl₂ | Chlorine | 70.9 | 10 | 3 | 20 | Gas |
| NH₃ | Ammonia | 17.0 | 300 | 200 | 1000 | Gas |
| H₂S | Hydrogen Sulfide | 34.1 | 50 | 30 | 100 | Gas |
| CH₄ | Methane | 16.0 | N/A* | 5000 | 50000 | Gas |
| C₃H₈ | Propane | 44.1 | 2100 | 5000 | 17000 | Liquefied gas |

*Simple asphyxiant — no toxicity-based IDLH.

#### Substance lookup & custom registration

```csharp
using CSharpNumerics.Physics.Materials.Chemical;

// Static instances
ChemicalSubstance cl = ChemicalSubstance.Chlorine;
ChemicalSubstance nh3 = ChemicalSubstance.Ammonia;

// Library lookup (case-insensitive)
ChemicalSubstance h2s = ChemicalLibrary.Get("H2S");
bool found = ChemicalLibrary.TryGet("CH4", out ChemicalSubstance methane);

// Properties
double mw = cl.MolarMass;          // 70.906 g/mol
double idlh = cl.IDLH;             // 10 ppm
double erpg2 = cl.ERPG2;           // 3 ppm
double vapDens = cl.VapourDensity;  // 2.49 (heavier than air)

// Register custom substance
ChemicalLibrary.Register(new ChemicalSubstance(
    "COCl2", "Phosgene", "75-44-5",
    98.92, 3.4, 8.3, PhaseAtSTP.Gas,
    idlh: 2, erpg2: 0.5, erpg3: 1.5, lc50: 5,
    tlvTwa: 0.1, tlvStel: 0.3));
```

#### Unit conversion (kg/m³ ↔ ppm)

```csharp
// Convert mass concentration to ppm at 20°C, 1 atm
double ppm = cl.KgM3ToPpm(2.95e-6);   // ≈ 1 ppm
double kgm3 = cl.PpmToKgM3(10);        // IDLH in kg/m³
```

#### GIS pipeline integration

```csharp
using CSharpNumerics.Physics.Materials;

// Attach chemical material to a dispersion scenario
var result = RiskScenario
    .ForGaussianPlume(5.0)
    .FromSource(new Vector(0, 0, 10))
    .WithWind(5, new Vector(1, 0, 0))
    .WithStability(StabilityClass.D)
    .WithMaterial(Materials.Chemical("Cl2"))    // ← chemical
    .OverGrid(new GeoGrid(-500, 500, -500, 500, 0, 100, 50))
    .OverTime(0, 3600, 60)
    .RunSingle();

// Snapshot layers: "ppm" and "toxicDose" (ppm·s)
double[] ppmLayer = result.Snapshots[0].GetLayer("ppm");
double[] toxDose  = result.Snapshots[0].GetLayer("toxicDose");

// IDLH exceedance polygon (Cl2 IDLH = 10 ppm)
var idlhZone = result.GeneratePeakExposurePolygon(10, "ppm");

// ERPG-2 exceedance polygon (3 ppm)
var erpg2Zone = result.GeneratePeakExposurePolygon(3, "ppm");

// Integrated toxic dose polygon
var integratedZone = result.GenerateIntegratedExposurePolygon(
    threshold: 1000, layerName: "toxicDose");  // 1000 ppm·s
```

---

### 🦠 Biological Materials — Bioaerosol Agents

**Namespace:** `CSharpNumerics.Physics.Materials.Biological`

Model biological aerosol agents (viruses, bacteria, spores) with viability decay, unit-mass conversion (kg/m³ ↔ units/m³), and screening-level infectious dose layers. Integrates with the GIS dispersion pipeline via `Materials.Biological("Virus")`.

#### Agent properties

| Property | Description |
|----------|-------------|
| `Code` | Short identifier for lookup (e.g. `"Virus"`) |
| `Name` | Human-readable name |
| `Classification` | `BiologicalAgentClass` enum: `Virus`, `Bacteria`, `Spore` |
| `TypicalDiameterMicrons` | Representative aerodynamic diameter (µm) |
| `UnitMassKg` | Mass of one biological unit (kg) — converts kg/m³ → units/m³ |
| `ViabilityHalfLifeSeconds` | Screening half-life for viability in air (seconds) |
| `IsPersistent` | `true` when viability is treated as non-decaying |

#### Built-in agents

| Code | Name | Class | Diameter (µm) | Unit mass (kg) | Viability t½ |
|------|------|-------|---------------|----------------|-------------|
| Virus | Generic Viral Aerosol | Virus | 0.12 | 1 × 10⁻¹⁸ | 6 h |
| Bacteria | Generic Bacterial Aerosol | Bacteria | 1.5 | 1 × 10⁻¹⁵ | 12 h |
| Spore | Generic Biological Spore | Spore | 2.5 | 3 × 10⁻¹⁵ | 7 d |

#### Agent lookup & custom registration

```csharp
using CSharpNumerics.Physics.Materials.Biological;

// Static instances
BiologicalAgent virus = BiologicalAgent.GenericVirus;
BiologicalAgent bacteria = BiologicalAgent.GenericBacteria;
BiologicalAgent spore = BiologicalAgent.GenericSpore;

// Library lookup (case-insensitive, supports aliases)
BiologicalAgent v = BiologicalLibrary.Get("virus");
BiologicalAgent s = BiologicalLibrary.Get("Spores");   // alias
bool found = BiologicalLibrary.TryGet("bacteria", out BiologicalAgent b);

// All registered agents (de-duplicated)
IReadOnlyCollection<BiologicalAgent> all = BiologicalLibrary.All();

// Register custom agent with aliases
BiologicalLibrary.Register(
    new BiologicalAgent(
        code: "Anthrax",
        name: "Bacillus anthracis spore",
        classification: BiologicalAgentClass.Spore,
        typicalDiameterMicrons: 1.5,
        unitMassKg: 1e-15,
        viabilityHalfLifeSeconds: 30 * 24 * 3600),  // 30 days
    "anthrax", "b.anthracis");
```

#### Unit conversion (kg/m³ ↔ units/m³)

```csharp
BiologicalAgent bacteria = BiologicalAgent.GenericBacteria;

// Convert mass concentration to biological units per m³
double units = bacteria.KgM3ToUnitsPerM3(2.5e-12);  // 2500 units/m³
double kgm3  = bacteria.UnitsPerM3ToKgM3(units);    // round-trip

// Viability fraction at a given time
double fraction = bacteria.ViabilityFractionAt(6 * 3600);  // after 6 h
```

#### GIS pipeline integration

```csharp
using CSharpNumerics.Physics.Materials;

// Attach biological material to a dispersion scenario
var sim = new PlumeSimulator(
    2.0, 8, new Vector(1, 0, 0), 50,
    new Vector(0, 0, 50), StabilityClass.D);
sim.Material = Materials.Biological("Bacteria");

var snaps = sim.Run(grid, tf);

// Snapshot layers: "bioUnits", "viableBioUnits", and "infectiousDose"
double[] bioUnits     = snaps[0].GetLayer("bioUnits");
double[] viable       = snaps[0].GetLayer("viableBioUnits");
double[] infectious   = snaps[0].GetLayer("infectiousDose");
```

The `viableBioUnits` layer applies exponential viability decay based on the agent's half-life, so viable counts decrease over successive time steps while raw `bioUnits` reflect total transported mass only.

---

### Kármán Vortex Street

The `KarmanVortexStreetExtensions` class provides formulas for analysing periodic vortex shedding behind bluff bodies (cylinders). It covers the Roshko empirical Strouhal–Reynolds correlation, vortex street geometry (von Kármán stability ratio), regime classification, and the idealised wake-drag formula.

| Formula | Method | Description |
|---------|--------|-------------|
| $f = \mathrm{St} \cdot U / D$ | `VortexSheddingFrequency` | Shedding frequency from Strouhal number |
| $\mathrm{St} \approx 0.198(1 - 19.7/\mathrm{Re})$ | `RoshkoStrouhalNumber` | Roshko (1954) empirical correlation |
| $f = \mathrm{St}(\mathrm{Re}) \cdot U / D$ | `VortexSheddingFrequencyFromReynolds` | Shedding frequency from Reynolds number |
| $\mathrm{Ro} = f D^2 / \nu$ | `RoshkoNumber` | Roshko number (= St · Re) |
| $a = U / f$ | `VortexStreamwiseSpacing` | Longitudinal spacing between vortices |
| $h/a = \tfrac{1}{\pi}\mathrm{arccosh}(\sqrt{2})$ | `StableSpacingRatio` | Von Kármán stability criterion ≈ 0.281 |
| $h = (h/a) \cdot a$ | `VortexLateralSpacing` | Cross-stream spacing between rows |
| $U_v \approx 0.875\,U$ | `VortexConvectionVelocity` | Downstream vortex convection speed |
| — | `IsPeriodicSheddingRegime` | True if 47 < Re < 2×10⁵ |
| — | `CylinderWakeRegime` | Wake regime string from Re |
| $C_d$ (wake momentum) | `VortexStreetDragCoefficient` | Idealised wake-drag from vortex spacing |

```csharp
using CSharpNumerics.Physics.FluidDynamics;

// Roshko's Strouhal–Reynolds correlation for a circular cylinder
double Re = 1000;
double St = Re.RoshkoStrouhalNumber();           // ≈ 0.194

// Vortex shedding frequency
double U = 5.0, D = 0.05;
double f = St.VortexSheddingFrequency(U, D);     // ≈ 19.4 Hz

// Or directly from Reynolds number
double f2 = Re.VortexSheddingFrequencyFromReynolds(U, D);

// Roshko number Ro = f·D²/ν
double nu = 1e-5;
double Ro = f.RoshkoNumber(D, nu);

// Vortex street geometry
double a = U.VortexStreamwiseSpacing(f);         // streamwise spacing
double h = a.VortexLateralSpacing();             // lateral spacing ≈ 0.281·a

// Vortex convection velocity
double Uv = U.VortexConvectionVelocity();        // ≈ 0.875·U

// Regime classification
bool shedding = Re.IsPeriodicSheddingRegime();   // true
string regime = Re.CylinderWakeRegime();         // "Turbulent transition"

// Wake-drag coefficient from vortex street
double Cd = h.VortexStreetDragCoefficient(a, U, Uv, D);
```

---

### Viscous Flow Model (`IViscousFlowModel`)

The `IViscousFlowModel` interface provides an abstraction for pipe flow physics used by simulation engines. The default implementation `ViscousFlowModel` encapsulates Hagen–Poiseuille physics.

| Method | Description |
|--------|-------------|
| `DrivingForce(dP/dx, ρ)` | Body force per unit mass $f = -(dP/dx) / \rho$ |
| `CylindricalDiffusion(ν, d²v/dr², dv/dr, r)` | $\nu [d^2v/dr^2 + (1/r)(dv/dr)]$ |
| `SymmetryAxisDiffusion(ν, d²v/dr²)` | $2\nu \, d^2v/dr^2$ (L'Hôpital at $r=0$) |

```csharp
using CSharpNumerics.Physics.FluidDynamics;
using CSharpNumerics.Physics.FluidDynamics.Interfaces;

IViscousFlowModel flow = new ViscousFlowModel();

double driving = flow.DrivingForce(pressureGradient: -100, density: 1000);
// 0.1 m/s² driving force

double diff = flow.CylindricalDiffusion(nu: 1e-6, d2vdr2: 50, dvdr: 10, r: 0.01);
```

---

## 🔩 Solid Mechanics Extensions

The `SolidExtensions` class provides solid mechanics calculations centred on the **Euler–Bernoulli beam equation** $EI u^{\prime\prime\prime\prime} = q$, plus Hooke's law, second moment of area, the flexure formula, and analytical beam deflections.

### Stress & Strain (Hooke's Law)

| Formula | Method | Description |
|---------|--------|-------------|
| $\sigma = F/A$ | `NormalStress` | Axial stress |
| $\varepsilon = \sigma/E$ | `NormalStrain` | Axial strain |
| $\sigma = E\varepsilon$ | `HookesLaw` | Hooke's law |
| $\tau = V/A$ | `ShearStress` | Average shear stress |
| $G = E / 2(1+\nu)$ | `ShearModulus` | Shear modulus |

```csharp
double sigma = force.NormalStress(area);        // σ = F/A
double eps   = sigma.NormalStrain(200e9);        // ε = σ/E  (steel)
double G     = (200e9).ShearModulus(0.3);         // ≈ 76.9 GPa
```

### Second Moment of Area

| Cross-Section | Formula | Method |
|---------------|---------|--------|
| Solid rectangle | $I = bh^3/12$ | `RectangularSecondMoment` |
| Solid circle | $I = \pi r^4/4$ | `CircularSecondMoment` |
| Hollow tube | $I = \pi(R^4 - r^4)/4$ | `TubularSecondMoment` |

```csharp
double I_rect = (0.10).RectangularSecondMoment(0.10);  // 8.33e-6 m⁴
double I_circ = (0.05).CircularSecondMoment();          // 4.91e-7 m⁴
double I_tube = (0.05).TubularSecondMoment(0.04);       // hollow tube
```

### Euler–Bernoulli Beam Equation

$$EI \frac{d^4 u}{dx^4} = q(x)$$

| Relation | Method | Description |
|----------|--------|-------------|
| $M = EI\kappa$ | `BendingMoment` | Moment from curvature |
| $\sigma = My/I$ | `BendingStress` | Flexure formula |
| $V = EI u'''$ | `BeamShearForce` | Shear from 3rd derivative |
| $q = EI u''''$ | `BeamLoadIntensity` | Load from 4th derivative |
| $r = EI d^4u/dx^4 - q$ | `EulerBernoulliResidual` | FD residual via `Biharmonic1D` |

```csharp
double M     = E.BendingMoment(I, curvature);          // M = EIκ
double sigma = M.BendingStress(y: 0.05, I);            // σ = My/I
double q     = E.BeamLoadIntensity(I, d4u);             // q = EIu⁗

// Discrete residual (should ≈ 0 for correct deflection)
var residual = u.EulerBernoulliResidual(E * I, dx, qVector);
```

### Analytical Beam Deflections

| Case | Max Deflection | Method |
|------|---------------|--------|
| Cantilever + point load P | $\delta = PL^3 / 3EI$ | `CantileverPointLoadMaxDeflection` |
| Cantilever + uniform q | $\delta = qL^4 / 8EI$ | `CantileverUniformLoadMaxDeflection` |
| Simply supported + uniform q | $\delta = 5qL^4 / 384EI$ | `SimplySupportedUniformLoadMaxDeflection` |
| Simply supported + midpoint P | $\delta = PL^3 / 48EI$ | `SimplySupportedPointLoadMaxDeflection` |

Full deflection curves are also available: `CantileverPointLoadDeflection(P, L, EI, x)`, etc.

```csharp
// Steel cantilever: 10×10 cm, 1 m, 1 kN tip load
double EI  = 200e9 * (0.10).RectangularSecondMoment(0.10);
double max = (1000.0).CantileverPointLoadMaxDeflection(1.0, EI);  // PL³/(3EI)

// Deflection curve
for (double x = 0; x <= 1.0; x += 0.1)
    Console.WriteLine($"x={x:F1}  u={1000.0.CantileverPointLoadDeflection(1.0, EI, x):E3} m");
```

### Beam Support Enum

The `BeamSupport` enum lives in `CSharpNumerics.Physics.SolidMechanics.Enums` and defines the three standard support conditions:

| Value | Description |
|-------|-------------|
| `Cantilever` | Fixed at one end, free at the other |
| `SimplySupported` | Supported at both ends with free rotation |
| `FixedFixed` | Fixed at both ends (no rotation or displacement) |

### Beam Model (`IBeamModel`)

The `IBeamModel` interface provides a structured abstraction over the analytical beam formulas, returning deflection, bending moment, and shear force as a tuple for any support condition and load type. The default implementation `BeamModel` delegates to `BeamExtensions`.

| Method | Description |
|--------|-------------|
| `CantileverPointLoad(P, a, L, EI, x)` | Point load at distance $a$ from fixed end |
| `CantileverUniformLoad(q, L, EI, x)` | Uniform distributed load |
| `SimplySupportedPointLoad(P, a, L, EI, x)` | Point load at distance $a$ from left support |
| `SimplySupportedUniformLoad(q, L, EI, x)` | Uniform distributed load |
| `FixedFixedPointLoad(P, a, L, EI, x)` | Point load on fixed-fixed beam |
| `FixedFixedUniformLoad(q, L, EI, x)` | Uniform load on fixed-fixed beam |
| `BendingStress(M, y, I)` | $\sigma = |M| \cdot y_{\max} / I$ |

All methods return `(double deflection, double moment, double shear)` — a complete set of beam responses at position $x$.

```csharp
using CSharpNumerics.Physics.SolidMechanics;
using CSharpNumerics.Physics.SolidMechanics.Interfaces;

IBeamModel beam = new BeamModel();

double EI = 200e9 * (0.10).RectangularSecondMoment(0.10);

// Cantilever with 1 kN tip load — query at midspan
var (d, m, s) = beam.CantileverPointLoad(P: 1000, a: 1.0, L: 1.0, EI: EI, x: 0.5);

// Simply supported with uniform 5 kN/m
var (d2, m2, s2) = beam.SimplySupportedUniformLoad(q: 5000, L: 2.0, EI: EI, x: 1.0);

// Bending stress at max fibre distance
double sigma = beam.BendingStress(moment: m2, halfHeight: 0.05, secondMoment: 8.33e-6);
```

---

## 🏭 Engineering Materials

The `EngineeringMaterial` immutable struct bundles thermo-mechanical-electrical-magnetic properties for multiphysics simulations. The `EngineeringLibrary` provides common pre-defined materials.

**Namespace:** `CSharpNumerics.Physics.Materials.Engineering`

### Material Properties

| Property | Unit | Description |
|----------|------|-------------|
| `ThermalConductivity` | W/(m·K) | Heat conduction coefficient $k$ |
| `SpecificHeat` | J/(kg·K) | Specific heat capacity $c_p$ |
| `Density` | kg/m³ | Mass density $\rho$ |
| `DynamicViscosity` | Pa·s | Dynamic viscosity $\mu$ |
| `ElectricPermittivity` | F/m | Dielectric permittivity $\varepsilon$ |
| `MagneticPermeability` | — | Relative magnetic permeability $\mu_r$ |
| `YoungsModulus` | Pa | Young's modulus $E$ |
| `PoissonsRatio` | — | Poisson's ratio $\nu$ |

**Computed properties:**

| Computed | Formula | Description |
|----------|---------|-------------|
| `ThermalDiffusivity` | $\alpha = k/(\rho c_p)$ | Heat diffusion rate |
| `KinematicViscosity` | $\nu = \mu/\rho$ | Flow diffusion rate |

### Pre-defined Materials

```csharp
using CSharpNumerics.Physics.Materials.Engineering;

// Metals
var steel = EngineeringLibrary.Steel;            // E=200 GPa, k=50 W/(m·K), μ_r=100
var al    = EngineeringLibrary.Aluminum;         // E=69 GPa, k=237 W/(m·K)
var cu    = EngineeringLibrary.Copper;           // E=120 GPa, k=401 W/(m·K)
var ti    = EngineeringLibrary.Titanium;         // E=116 GPa, k=21.9 W/(m·K)
var brass = EngineeringLibrary.Brass;            // E=100 GPa, k=109 W/(m·K)
var ss    = EngineeringLibrary.StainlessSteel;   // E=193 GPa, k=16 W/(m·K)

// Fluids
var water = EngineeringLibrary.Water;            // μ=1e-3 Pa·s, ρ=998 kg/m³
var air   = EngineeringLibrary.Air;              // μ=1.81e-5 Pa·s
var oil   = EngineeringLibrary.Oil;              // μ=0.03 Pa·s, ρ=870 kg/m³
var glyc  = EngineeringLibrary.Glycerin;         // μ=1.412 Pa·s, ρ=1261 kg/m³

// Construction / Structural
var conc  = EngineeringLibrary.Concrete;         // E=30 GPa
var glass = EngineeringLibrary.Glass;            // E=70 GPa
var wood  = EngineeringLibrary.Wood;             // E=12 GPa, k=0.15 W/(m·K)
var rub   = EngineeringLibrary.Rubber;           // E=0.01 GPa, ν=0.49
var hdpe  = EngineeringLibrary.Plastic;          // E=1.1 GPa (HDPE)

double alpha = steel.ThermalDiffusivity;         // k/(ρ·cp)
double nu    = water.KinematicViscosity;         // μ/ρ
double muR   = steel.MagneticPermeability;       // 100.0 (ferromagnetic)
```

### Custom Materials

```csharp
var iron = new EngineeringMaterial(
    name: "Iron",
    thermalConductivity: 80,          // W/(m·K)
    specificHeat: 450,                // J/(kg·K)
    density: 7874,                    // kg/m³
    dynamicViscosity: 0,              // not a fluid
    electricPermittivity: 0,
    youngsModulus: 211e9,             // Pa
    poissonsRatio: 0.29,
    magneticPermeability: 5000.0);    // ferromagnetic
```

`magneticPermeability` defaults to 1.0 (vacuum/non-magnetic) when omitted.

---

## ⚛️ Quantum Gates

**Namespace:** `CSharpNumerics.Physics.Quantum`

Quantum gates are unitary operators acting on one or more qubits. Each gate is represented by a `ComplexMatrix` and knows how to apply itself to a `ComplexVectorN` state vector.

All gates extend the abstract `QuantumGate` base class:

| Class | Qubits | Matrix | Description |
|---|---|---|---|
| `HadamardGate` | 1 | $\frac{1}{\sqrt{2}}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$ | Creates equal superposition |
| `PauliXGate` | 1 | $\begin{pmatrix}0&1\\1&0\end{pmatrix}$ | Quantum NOT — flips \|0⟩ ↔ \|1⟩ |
| `PauliYGate` | 1 | $\begin{pmatrix}0&-i\\i&0\end{pmatrix}$ | Bit + phase flip, $Y^2 = I$ |
| `PauliZGate` | 1 | $\begin{pmatrix}1&0\\0&-1\end{pmatrix}$ | Phase flip — \|1⟩ → −\|1⟩ |
| `SGate` | 1 | $\begin{pmatrix}1&0\\0&i\end{pmatrix}$ | Phase gate (π/2), $S^2 = Z$ |
| `TGate` | 1 | $\begin{pmatrix}1&0\\0&e^{i\pi/4}\end{pmatrix}$ | π/8 gate, $T^2 = S$ |
| `PhaseGate(θ)` | 1 | $\begin{pmatrix}1&0\\0&e^{i\theta}\end{pmatrix}$ | General phase gate, $P(\pi) = Z$ |
| `RxGate(θ)` | 1 | $\begin{pmatrix}\cos\frac{\theta}{2}&-i\sin\frac{\theta}{2}\\-i\sin\frac{\theta}{2}&\cos\frac{\theta}{2}\end{pmatrix}$ | Rotation about X-axis |
| `RyGate(θ)` | 1 | $\begin{pmatrix}\cos\frac{\theta}{2}&-\sin\frac{\theta}{2}\\\sin\frac{\theta}{2}&\cos\frac{\theta}{2}\end{pmatrix}$ | Rotation about Y-axis |
| `RzGate(θ)` | 1 | $\begin{pmatrix}e^{-i\theta/2}&0\\0&e^{i\theta/2}\end{pmatrix}$ | Rotation about Z-axis |
| `CNOTGate` | 2 | 4×4 permutation | Flips target qubit when control is \|1⟩ |
| `CZGate` | 2 | $\text{diag}(1,1,1,-1)$ | Phase flip on \|11⟩ |
| `CPhaseGate(θ)` | 2 | $\text{diag}(1,1,1,e^{i\theta})$ | Controlled phase, $CP(\pi) = CZ$ |
| `SWAPGate` | 2 | 4×4 permutation | Swaps two qubit states |
| `ToffoliGate` | 3 | 8×8 permutation | CCNOT — flips target when both controls are \|1⟩ |
| `FredkinGate` | 3 | 8×8 permutation | CSWAP — swaps targets when control is \|1⟩ |
| `PhaseOracle(n, states)` | n | 2ⁿ×2ⁿ diagonal | Flips phase of marked basis states: \|w⟩ → −\|w⟩ |
| `ControlledGate(U)` | n+1 | 2ⁿ⁺¹×2ⁿ⁺¹ block | Applies U when control qubit is \|1⟩ |
| `ModularMultiplyGate(a,N,n)` | n | 2ⁿ×2ⁿ permutation | \|y⟩ → \|ay mod N⟩, used in Shor's algorithm |

### QuantumGate (abstract)

```csharp
public abstract class QuantumGate
{
    public abstract int QubitCount { get; }
    public abstract ComplexMatrix GetMatrix();
    public ComplexVectorN Apply(ComplexVectorN amplitudes, int[] qubitIndices, int totalQubits);
}
```

`Apply` uses the gate's unitary matrix to transform the relevant amplitudes in-place (tensor-product expansion) — works for arbitrary qubit counts and target indices.

### Usage

Gates are consumed by `QuantumInstruction` and `QuantumCircuit` in the Engines Quantum section:

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;

var circuit = new QuantumCircuit(2);
circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
circuit.AddInstruction(new QuantumInstruction(new CNOTGate(), new List<int> { 0, 1 }));

var state = new QuantumSimulator().Run(circuit);   // Bell state (|00⟩ + |11⟩)/√2
```

### BlochVector

Represents a single-qubit pure state on the Bloch sphere. Given $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$:

$$x = 2\,\text{Re}(\alpha^*\beta), \quad y = 2\,\text{Im}(\alpha^*\beta), \quad z = |\alpha|^2 - |\beta|^2$$

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;

// From a circuit result
var circuit = new QuantumCircuit(1);
circuit.AddInstruction(new QuantumInstruction(new HadamardGate(), new List<int> { 0 }));
var state = new QuantumSimulator().Run(circuit);

BlochVector bloch = state.GetBlochVector();  // (1, 0, 0) — +X axis
double theta = bloch.Theta;                  // polar angle θ
double phi   = bloch.Phi;                    // azimuthal angle φ
double r     = bloch.Radius;                 // 1.0 for pure states
Vector v     = bloch.ToVector();             // 3D Vector for rendering

// Directly from amplitudes
var b = BlochVector.FromAmplitudes(
    new ComplexNumber(1, 0), new ComplexNumber(0, 0));  // |0⟩ → (0, 0, 1)
```

| Property | Description |
|---|---|
| `X`, `Y`, `Z` | Cartesian Bloch coordinates |
| `Theta` | Polar angle $\theta \in [0, \pi]$ from +Z |
| `Phi` | Azimuthal angle $\varphi \in (-\pi, \pi]$ |
| `Radius` | Vector length (1 for pure states) |
| `ToVector()` | Returns a 3D `Vector` for visualization |

**Canonical states on the sphere:**

| State | Bloch |
|---|---|
| $|0\rangle$ | $(0, 0, 1)$ — north pole |
| $|1\rangle$ | $(0, 0, -1)$ — south pole |
| $(|0\rangle+|1\rangle)/\sqrt{2}$ | $(1, 0, 0)$ — +X |
| $(|0\rangle-|1\rangle)/\sqrt{2}$ | $(-1, 0, 0)$ — −X |
| $(|0\rangle+i|1\rangle)/\sqrt{2}$ | $(0, 1, 0)$ — +Y |
| $(|0\rangle-i|1\rangle)/\sqrt{2}$ | $(0, -1, 0)$ — −Y |

### QuantumFidelity

Static methods for computing the overlap between quantum states. Fidelity ranges from 0 (orthogonal) to 1 (identical).

**State fidelity** — $F = |\langle\psi|\phi\rangle|^2$

```csharp
using CSharpNumerics.Physics.Quantum;
using CSharpNumerics.Engines.Quantum;

var sim = new QuantumSimulator();
var s0    = sim.Run(QuantumCircuitBuilder.New(1).Build());        // |0⟩
var sPlus = sim.Run(QuantumCircuitBuilder.New(1).H(0).Build());  // |+⟩
var s1    = sim.Run(QuantumCircuitBuilder.New(1).X(0).Build());  // |1⟩

double f1 = QuantumFidelity.Fidelity(s0, s0);     // 1.0  — identical
double f2 = QuantumFidelity.Fidelity(s0, sPlus);  // 0.5  — overlap
double f3 = QuantumFidelity.Fidelity(s0, s1);     // 0.0  — orthogonal
```

**Vector fidelity** — works directly on `ComplexVectorN`

```csharp
double f = QuantumFidelity.Fidelity(psi, phi);   // |⟨ψ|φ⟩|²
```

**Bloch fidelity** — geometric formula for single-qubit states: $F = \frac{1}{2}(1 + \hat{n}_1 \cdot \hat{n}_2)$

```csharp
var b1 = s1State.GetBlochVector();
var b2 = s2State.GetBlochVector();
double f = QuantumFidelity.BlochFidelity(b1, b2);
// Matches state fidelity for pure single-qubit states
```

### Noise Models

The `Physics.Quantum.NoiseModels` namespace provides quantum noise channels using the **Kraus operator** formalism. Each channel implements `INoiseChannel` and satisfies the completeness relation $\sum E_k^\dagger E_k = I$.

| Channel | Parameter | Kraus Operators | Physical Model |
|---|---|---|---|
| `DepolarizingNoise(p)` | $p \in [0,1]$ | $E_0 = \sqrt{1 - 3p/4}\,I$, $E_{1,2,3} = \sqrt{p/4}\,\{X, Y, Z\}$ | Random Pauli error — qubit replaced by maximally mixed state with probability $p$ |
| `DephasingNoise(p)` | $p \in [0,1]$ | $E_0 = \sqrt{1-p}\,I$, $E_1 = \sqrt{p}\,Z$ | Phase-flip error — $T_2$ decoherence |
| `AmplitudeDampingNoise(\gamma)` | $\gamma \in [0,1]$ | $E_0 = [[1,0],[0,\sqrt{1-\gamma}]]$, $E_1 = [[0,\sqrt{\gamma}],[0,0]]$ | Energy dissipation — $\|1\rangle \to \|0\rangle$ decay ($T_1$) |

**Usage with NoisyQuantumSimulator** (from `CSharpNumerics.Engines.Quantum`):

```csharp
using CSharpNumerics.Physics.Quantum.NoiseModels;
using CSharpNumerics.Engines.Quantum;

var circuit = QuantumCircuitBuilder.New(2).H(0).CNOT(0, 1).Build();

// Ideal simulation
var ideal = new QuantumSimulator().Run(circuit);

// Noisy simulation — channels stacked, applied after each gate
var noisy = new NoisyQuantumSimulator(new Random(42))
    .WithNoise(new DepolarizingNoise(0.01))
    .WithNoise(new AmplitudeDampingNoise(0.005))
    .Run(circuit);

double fidelity = QuantumFidelity.Fidelity(ideal, noisy);
```

**INoiseChannel interface**

```csharp
public interface INoiseChannel
{
    int QubitCount { get; }
    ComplexMatrix[] GetKrausOperators();
}
```

### Module Structure

```
Physics/Quantum/
├── QuantumGate.cs       Abstract base with Apply logic
├── QuantumFidelity.cs   State fidelity metrics
├── BlochVector.cs       Bloch sphere representation
├── Interfaces/
│   └── INoiseChannel.cs Noise channel contract
├── NoiseModels/
│   ├── DepolarizingNoise.cs      Depolarizing channel
│   ├── DephasingNoise.cs         Dephasing (phase-flip) channel
│   └── AmplitudeDampingNoise.cs  Amplitude damping (T₁ decay)
├── ErrorCorrection/
│   ├── Interfaces/
│   │   └── IQuantumErrorCorrectionCode.cs  QEC code contract
│   └── Codes/
│       ├── BitFlipCode3.cs       [[3,1,1]] bit-flip code
│       ├── PhaseFlipCode3.cs     [[3,1,1]] phase-flip code
│       ├── ShorCode9.cs          [[9,1,3]] Shor code
│       └── SteaneCode7.cs        [[7,1,3]] Steane (CSS) code
├── HadamardGate.cs      H gate
├── PauliXGate.cs        X gate (Pauli-X)
├── PauliZGate.cs        Z gate (Pauli-Z)
├── SGate.cs             S gate (Phase, π/2)
├── TGate.cs             T gate (π/8)
├── RxGate.cs            Rotation about X-axis
├── RyGate.cs            Rotation about Y-axis
├── RzGate.cs            Rotation about Z-axis
├── CNOTGate.cs          Controlled-NOT gate
├── CZGate.cs            Controlled-Z gate
└── SWAPGate.cs          SWAP gate
```

### Quantum Error Correction — Code Definitions

**Namespace:** `CSharpNumerics.Physics.Quantum.ErrorCorrection`

The `ErrorCorrection` sub-namespace defines quantum error-correcting codes as mathematical objects — stabilizers, correction maps, and logical operators. The simulation/orchestration layer lives in `Engines.Quantum.ErrorCorrection`.

#### IQuantumErrorCorrectionCode

Interface for any [[n, k, d]] stabilizer code:

| Property / Method | Description |
|---|---|
| `PhysicalQubits` | Number of physical qubits (n) |
| `LogicalQubits` | Number of logical qubits (k) |
| `Distance` | Code distance (d) |
| `SyndromeQubits` | Number of ancilla qubits for syndrome extraction |
| `GetStabilizers()` | Returns stabilizer generators as (qubit, Pauli) lists |
| `GetCorrectionMap()` | Syndrome integer → corrective (qubit, Pauli) operations |
| `GetLogicalX(i)` | Logical X̄ operator for logical qubit i |
| `GetLogicalZ(i)` | Logical Z̄ operator for logical qubit i |

#### BitFlipCode3 — [[3,1,1]]

Encodes $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$ into $\alpha|000\rangle + \beta|111\rangle$. Corrects any single-qubit X (bit-flip) error.

**Stabilizers:** $Z_0 Z_1$, $Z_1 Z_2$

| Syndrome | Error | Correction |
|---|---|---|
| 00 | None | — |
| 01 | $X_0$ | Apply X to qubit 0 |
| 10 | $X_2$ | Apply X to qubit 2 |
| 11 | $X_1$ | Apply X to qubit 1 |

**Logical operators:** $\bar{X} = X_0 X_1 X_2$, $\bar{Z} = Z_0$

#### PhaseFlipCode3 — [[3,1,1]]

Encodes $|\psi\rangle$ into $\alpha|{+}{+}{+}\rangle + \beta|{-}{-}{-}\rangle$. Corrects any single-qubit Z (phase-flip) error.

**Stabilizers:** $X_0 X_1$, $X_1 X_2$

| Syndrome | Error | Correction |
|---|---|---|
| 00 | None | — |
| 01 | $Z_0$ | Apply Z to qubit 0 |
| 10 | $Z_2$ | Apply Z to qubit 2 |
| 11 | $Z_1$ | Apply Z to qubit 1 |

**Logical operators:** $\bar{X} = X_0$, $\bar{Z} = Z_0 Z_1 Z_2$

```csharp
using CSharpNumerics.Physics.Quantum.ErrorCorrection;

var code = new BitFlipCode3();
var stabs = code.GetStabilizers();       // [{(0,'Z'),(1,'Z')}, {(1,'Z'),(2,'Z')}]
var map   = code.GetCorrectionMap();     // 0→[], 1→[(0,'X')], 2→[(2,'X')], 3→[(1,'X')]
var logX  = code.GetLogicalX();          // [(0,'X'),(1,'X'),(2,'X')]
```

#### ShorCode9 — [[9,1,3]]

Shor's 9-qubit code — the first code to correct **any** single-qubit error (X, Z, or Y). Uses three blocks of three qubits: inner bit-flip repetition within blocks, outer phase-flip repetition across blocks.

$$|0\rangle_L = \frac{(|000\rangle+|111\rangle)(|000\rangle+|111\rangle)(|000\rangle+|111\rangle)}{2\sqrt{2}}$$
$$|1\rangle_L = \frac{(|000\rangle-|111\rangle)(|000\rangle-|111\rangle)(|000\rangle-|111\rangle)}{2\sqrt{2}}$$

**8 stabilizer generators:**

| Generator | Operator | Type |
|---|---|---|
| $g_1 \ldots g_6$ | $Z_i Z_{i+1}$ within each block | Bit-flip detection |
| $g_7$ | $X_0 X_1 X_2 X_3 X_4 X_5$ | Phase-flip detection (blocks 0–1) |
| $g_8$ | $X_3 X_4 X_5 X_6 X_7 X_8$ | Phase-flip detection (blocks 1–2) |

**Correction:** The 8-bit syndrome (256 values) maps to at most one X correction + one Z correction. Each block's 2-bit sub-syndrome identifies bit-flip errors exactly as in `BitFlipCode3`. The 2-bit phase-flip sub-syndrome identifies which block suffered a Z error.

**Logical operators:** $\bar{X} = X_0 X_1 \ldots X_8$, $\bar{Z} = Z_0 Z_3 Z_6$

```csharp
var shor = new ShorCode9();
var stabs = shor.GetStabilizers();  // 8 generators
var map   = shor.GetCorrectionMap(); // 256 syndrome entries
```

#### SteaneCode7 — [[7,1,3]]

The Steane code — the smallest CSS (Calderbank-Shor-Steane) code, correcting **any** single-qubit error using only 7 physical qubits. Built from the classical [7,4,3] Hamming code and its dual [7,3,4] code.

$$|0\rangle_L = \frac{1}{\sqrt{8}} \sum_{x \in C_2} |x\rangle \qquad |1\rangle_L = \frac{1}{\sqrt{8}} \sum_{x \in C_2 + v} |x\rangle$$

where $C_2$ is the [7,3,4] dual Hamming code and $v = (1110000)$.

**6 stabilizer generators:**

| Generator | Operator | Type |
|---|---|---|
| $g_1$ | $Z_3 Z_4 Z_5 Z_6$ | X-error detection |
| $g_2$ | $Z_1 Z_2 Z_5 Z_6$ | X-error detection |
| $g_3$ | $Z_0 Z_2 Z_4 Z_6$ | X-error detection |
| $g_4$ | $X_3 X_4 X_5 X_6$ | Z-error detection |
| $g_5$ | $X_1 X_2 X_5 X_6$ | Z-error detection |
| $g_6$ | $X_0 X_2 X_4 X_6$ | Z-error detection |

**Correction:** The 6-bit syndrome has two independent 3-bit halves — Z-stabilizer syndrome (bits 0–2) identifies X errors, X-stabilizer syndrome (bits 3–5) identifies Z errors, using the Hamming decoding map: syndrome $s \to$ qubit $j$ where column $j$ of the parity-check matrix equals the binary expansion of $s$.

**Logical operators:** $\bar{X} = X_0 X_1 \ldots X_6$, $\bar{Z} = Z_0 Z_1 \ldots Z_6$

```csharp
var steane = new SteaneCode7();
var stabs  = steane.GetStabilizers();   // 6 generators (3 Z-type + 3 X-type)
var map    = steane.GetCorrectionMap();  // 64 syndrome entries
```

---

## ✈️ Aerodynamics

The `Physics.FluidDynamics.Aerodynamics` namespace provides models for atmospheric flight — ISA atmosphere, airfoil lift/drag, control surfaces, and propulsion. These classes are used by the Game Engine's flight dynamics system but are pure physics with no engine dependency.

### Atmosphere Model (ISA)

International Standard Atmosphere from sea level to 86 km with all seven standard layers.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Thermodynamic properties at any altitude
double T   = AtmosphereModel.Temperature(11000);      // 216.65 K (tropopause)
double P   = AtmosphereModel.Pressure(11000);          // ~22 632 Pa
double rho = AtmosphereModel.Density(5000);             // ~0.736 kg/m³
double a   = AtmosphereModel.SpeedOfSound(0);           // 340.29 m/s at sea level

// Transport properties (Sutherland's law)
double mu  = AtmosphereModel.DynamicViscosity(10000);   // Pa·s
double nu  = AtmosphereModel.KinematicViscosity(10000);  // m²/s

// Sea-level reference
double rho0 = AtmosphereModel.SeaLevelDensity;          // 1.225 kg/m³
```

| Layer | Altitude | Lapse Rate |
|-------|----------|------------|
| Troposphere | 0–11 km | −6.5 K/km |
| Tropopause | 11–20 km | Isothermal |
| Stratosphere 1 | 20–32 km | +1.0 K/km |
| Stratosphere 2 | 32–47 km | +2.8 K/km |
| Stratopause | 47–51 km | Isothermal |
| Mesosphere 1 | 51–71 km | −2.8 K/km |
| Mesosphere 2 | 71–86 km | −2.0 K/km |

### Airfoil Model

Lift and drag coefficient lookup as a function of angle of attack (AoA). Supports built-in profiles and custom lookup tables with linear interpolation.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Built-in NACA symmetric airfoil (stall modelled)
var airfoil = AirfoilModel.NACASymmetric();
double cl = airfoil.Cl(5 * Math.PI / 180);    // Cl at 5° AoA
double cd = airfoil.Cd(5 * Math.PI / 180);    // Cd at 5° AoA
double ld = airfoil.LiftToDrag(0.05);          // L/D ratio

// Built-in flat plate model
var plate = AirfoilModel.FlatPlate(cd0: 0.02);

// Custom airfoil from wind-tunnel data
var custom = new AirfoilModel("MyWing",
    alphas: new[] { -0.3, -0.1, 0.0, 0.1, 0.2, 0.3 },
    cls:    new[] { -0.8, -0.3, 0.0, 0.3, 0.6, 0.4 },
    cds:    new[] { 0.05, 0.02, 0.01, 0.02, 0.04, 0.08 });
```

**Stall behaviour:** The NACA symmetric model includes post-stall Cl reduction — Cl rises linearly up to ~15° AoA then drops as flow separates.

### Control Surface

Models the aerodynamic effect of elevator, aileron, rudder, and flap deflections.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Built-in presets
var elevator = ControlSurface.Elevator();
var aileron  = ControlSurface.Aileron();
var rudder   = ControlSurface.Rudder();

// Incremental coefficients from deflection
double dCl = elevator.DeltaCl(0.1);    // ΔCl = τ · Cl_δ · δ
double dCd = elevator.DeltaCd(0.1);    // ΔCd = Cd_δ · δ²

// Deflection is clamped to mechanical limits
double maxDefl = elevator.MaxDeflection;  // ~25° (0.4363 rad)
```

| Property | Description |
|----------|-------------|
| `ClDelta` | Lift sensitivity dCl/dδ (1/rad) |
| `CdDelta` | Drag sensitivity dCd/dδ² (1/rad²) |
| `Tau` | Flap effectiveness factor (0–1) |
| `Axis` | Primary control axis (Pitch/Roll/Yaw) |

### Propulsion Model

Engine thrust as a function of throttle, altitude, and airspeed. Supports jet and propeller types.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Jet engine: T = T_max · throttle · (ρ/ρ₀)^n
var jet = PropulsionModel.Jet(maxThrustSeaLevel: 50000, engineCount: 2);
double thrust = jet.Thrust(throttle: 1.0, altitude: 0, airspeed: 100);  // ~100 kN

// Propeller engine: T = P·η/V, capped at static thrust
var prop = PropulsionModel.Propeller(maxPower: 200000, staticThrust: 5000);
double tProp = prop.Thrust(throttle: 0.8, altitude: 1000, airspeed: 50);

// Thrust decreases with altitude (density lapse)
double tHigh = jet.Thrust(1.0, 10000, 200);  // < thrust at sea level
```

### NACA Airfoil Geometry

Generate 2D surface coordinates for NACA 4-digit series airfoils. Uses cosine spacing for leading-edge resolution.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Generate NACA 0012 with 100 panels
var (x, y) = NACAGeometry.Generate("0012", numPanels: 100, chord: 1.0);

// Cambered airfoil (2% camber at 40% chord, 12% thickness)
var (xc, yc) = NACAGeometry.Generate("2412", numPanels: 80);

// Symmetric shorthand
var (xs, ys) = NACAGeometry.GenerateSymmetric(thicknessPercent: 15);

// Rotate for angle of attack (about quarter-chord)
double alpha = 5 * Math.PI / 180;
var (xr, yr) = NACAGeometry.Rotate(x, y, alpha);
```

| Parameter | Description |
|-----------|-------------|
| 1st digit | Max camber (% of chord) |
| 2nd digit | Location of max camber (tenths of chord) |
| 3rd–4th digits | Max thickness (% of chord) |

### Panel Method (Hess-Smith)

Solves 2D incompressible **potential flow** around arbitrary closed bodies using the Hess-Smith panel method with constant-strength source panels and uniform vortex distribution.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Generate airfoil and solve
var (x, y) = NACAGeometry.Generate("0012", 100);
double alpha = 5 * Math.PI / 180;
var (xr, yr) = NACAGeometry.Rotate(x, y, alpha);

var result = PanelMethod.Solve(xr, yr, alpha, freestream: 1.0);

// Surface pressure coefficient
double[] cp = result.Cp;         // Cp on each panel
double cl   = result.Cl;         // Integrated lift coefficient
double cm   = result.Cm;         // Quarter-chord moment coefficient
double gamma = result.Gamma;     // Circulation strength

// Velocity field at arbitrary points
var queryX = new double[] { -2.0, 0.5, 2.0 };
var queryY = new double[] { 0.5, 0.5, 0.5 };
var (u, v) = PanelMethod.VelocityField(result, xr, yr, queryX, queryY);
```

**Algorithm:**
1. Discretize airfoil into N flat panels with source strength $\sigma_i$ and uniform vortex $\gamma$
2. Apply N no-penetration boundary conditions on panel midpoints
3. Enforce Kutta condition at trailing edge (smooth flow departure)
4. Solve $(N+1) \times (N+1)$ linear system for $\sigma_i$ and $\gamma$
5. Compute tangential velocity $V_t$ and pressure coefficient $C_p = 1 - (V_t / U_\infty)^2$

| Output | Description |
|--------|-------------|
| `Cp` | Surface pressure coefficient distribution |
| `Vt` | Surface tangential velocity |
| `Cl` | Lift coefficient (integrated) |
| `Cm` | Moment coefficient about quarter-chord |
| `Gamma` | Vortex strength (circulation) |
| `Sigma` | Source strengths per panel |

---

## 🌪️ Turbulence Model (Smagorinsky SGS)

Subgrid-scale turbulence model for Large Eddy Simulation. Computes eddy viscosity from the local strain rate: $\nu_t = (C_s \cdot \Delta)^2 \cdot |S|$

```csharp
using CSharpNumerics.Physics.FluidDynamics.Turbulence;

var turb = new TurbulenceModel(gridSpacing: 0.01, cs: 0.17);

// Eddy viscosity from strain rate
double nuT = turb.EddyViscosity(strainRateMagnitude: 50.0);

// Total effective viscosity: ν + ν_t
double nuEff = turb.EffectiveViscosity(molecularViscosity: 1.5e-5, strainRateMagnitude: 50.0);

// Compute strain rate from velocity gradients
double s2d = TurbulenceModel.StrainRate2D(dudx: 10, dudy: 2, dvdx: 3, dvdy: -10);
double s3d = TurbulenceModel.StrainRate3D(
    dudx: 10, dudy: 2, dudz: 0,
    dvdx: 3, dvdy: -10, dvdz: 1,
    dwdx: 0, dwdy: -1, dwdz: 0);
```

---

## 🔥 Buoyancy Force (Boussinesq)

Thermal buoyancy for smoke, fire, and hot gas plume simulations. Uses the Boussinesq approximation where density variations only appear in the gravity term.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Buoyancy;

// Buoyancy acceleration: a = g · β · (T − T_ref)
double a = BuoyancyForce.BoussinesqAcceleration(
    temperature: 400, referenceTemperature: 300);  // hot gas → positive (upward)

// Apply as velocity impulse
double vNew = BuoyancyForce.ApplyBuoyancy(
    velocityZ: 0, temperature: 500, referenceTemperature: 300, dt: 0.01);

// Density-based buoyancy: a = g · (ρ_ref − ρ) / ρ_ref
double aDensity = BuoyancyForce.DensityBuoyancy(density: 0.5, referenceDensity: 1.0);
```

---

## 🧵 Soft Body Physics

The `Physics.Mechanics.SoftBody` namespace provides mass-spring deformable meshes and cloth simulation using Verlet integration with constraint projection.

### DeformableMesh

Mass-spring network on a triangulated mesh. Each vertex is a point mass connected to neighbours by springs.

```csharp
using CSharpNumerics.Physics.Mechanics.SoftBody;
using CSharpNumerics.Numerics.Objects;

// Create a rectangular grid mesh
var mesh = DeformableMesh.CreateGrid(
    width: 2.0, height: 2.0, resX: 20, resY: 20,
    mass: 0.1, stiffness: 0.9, origin: new Vector(0, 0, 5));

// Pin corners
mesh.Pin(0);
mesh.Pin(19);

// Step simulation
for (int i = 0; i < 1000; i++)
    mesh.Step(dt: 0.001);

// Collision with sphere
mesh.CollideWithSphere(center: new Vector(1, 1, 3), radius: 0.5);

// Collision with ground plane
mesh.CollideWithGround(height: 0);
```

| Property | Description |
|----------|-------------|
| `Gravity` | External acceleration (default (0,0,−9.81)) |
| `Damping` | Velocity damping 0–1 (default 0.99) |
| `Iterations` | Constraint iterations per step (default 5) |

**Custom mesh from triangle soup:**

```csharp
var verts = new Vector[] { /* positions */ };
var mesh = new DeformableMesh(verts, mass: 0.1);
mesh.GenerateSpringsFromTriangles(triangleIndices, stiffness: 0.8);
```

### ClothSimulation

Built on `DeformableMesh` with wind force and optional self-collision.

```csharp
using CSharpNumerics.Physics.Mechanics.SoftBody;

var cloth = new ClothSimulation(
    width: 3, height: 2, resX: 30, resY: 20,
    mass: 0.05, stiffness: 0.95);

cloth.PinTopEdge();                  // fixed along top row
cloth.Wind = new Vector(2, 0, 0.5);  // wind blowing in +X

for (int i = 0; i < 500; i++)
    cloth.Step(0.002);

Vector pos = cloth.GetPosition(15, 10);  // query vertex position
```

| Property | Description |
|----------|-------------|
| `Wind` | Wind force applied to each particle |
| `EnableSelfCollision` | Toggle self-collision (default false) |
| `SelfCollisionRadius` | Distance threshold for self-collision |

---

## 💦 SPH Solver (Smoothed Particle Hydrodynamics)

Lagrangian fluid simulation for water splashes, liquid in containers, and free-surface flows. Based on Müller et al. 2003 with Poly6/Spiky/Viscosity kernels.

```csharp
using CSharpNumerics.Physics.FluidDynamics.SPH;
using CSharpNumerics.Numerics.Objects;

// Create particles in a block
var positions = new List<Vector>();
for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++)
        for (int k = 0; k < 10; k++)
            positions.Add(new Vector(i * 0.02, j * 0.02, k * 0.02 + 0.5));

var sph = new SPHSolver(positions);
sph.SmoothingRadius = 0.04;
sph.ParticleMass = 0.02;
sph.RestDensity = 1000;
sph.GasConstant = 2000;
sph.Viscosity = 1.0;
sph.Gravity = new Vector(0, 0, -9.81);
sph.BoundsMin = new Vector(-0.5, -0.5, 0);
sph.BoundsMax = new Vector(0.5, 0.5, 1);

// Step simulation
sph.Step(numSteps: 100);

// Query particles
foreach (var p in sph.Particles)
    Console.WriteLine($"pos={p.Position} vel={p.Velocity} density={p.Density:F1}");
```

| Property | Description |
|----------|-------------|
| `SmoothingRadius` | Kernel support radius h |
| `ParticleMass` | Mass per particle |
| `RestDensity` | Reference fluid density (kg/m³) |
| `GasConstant` | Pressure stiffness (higher = less compressible) |
| `Viscosity` | Dynamic viscosity coefficient |
| `BoundaryRestitution` | Bounce coefficient at domain walls |

---

## 🌊 VOF Tracker (Volume of Fluid)

Eulerian free-surface tracking using the Volume-of-Fluid method. Each cell stores a volume fraction $F \in [0,1]$: $F=0$ is empty, $F=1$ is full liquid, $0 < F < 1$ is an interface cell.

```csharp
using CSharpNumerics.Physics.FluidDynamics.FreeSurface;

var vof = new VOFTracker(nx: 64, ny: 64, cellSize: 0.1);

// Fill a rectangular region with liquid
vof.FillRect(x0: 10, y0: 0, x1: 30, y1: 20);

// Fill a circular blob
vof.FillCircle(cx: 3.2, cy: 4.0, radius: 0.8);

// Advect with velocity field
vof.Advect(u: velocityX, v: velocityY, dt: 0.01);

// Query volume fraction
double f = vof.GetFraction(15, 10);

// Estimate surface normal at interface cells
var (nx, ny) = vof.EstimateNormal(15, 10);

// Total liquid volume (conserved quantity)
double vol = vof.TotalVolume();
```

| Property | Description |
|----------|-------------|
| `CellSize` | Size of each grid cell |
| `LiquidThreshold` | Minimum F to count as "liquid present" |

---

## 🚀 Rocket Propulsion

The `Physics.Mechanics.Propulsion` namespace provides a complete rocket engine and guidance model — thrust, propellant management, variable-mass dynamics, center of mass tracking, and launch guidance algorithms.

### RocketEngine

Models a bipropellant rocket engine with throttle and gimbal control. ISP varies linearly between sea-level and vacuum values based on ambient pressure.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

// Merlin-1D-class engine
var engine = new RocketEngine(
    ispSeaLevel: 282,
    ispVacuum: 311,
    maxThrustVacuum: 845000,
    minThrottle: 0.4,
    maxGimbalAngle: 5.0 * Math.PI / 180);

engine.SetThrottle(0.9);
engine.SetGimbal(pitch: 0.02, yaw: -0.01);

double thrust = engine.Thrust(pressureRatio: 0.5);    // N at 50% sea-level pressure
double isp    = engine.EffectiveIsp(pressureRatio: 0); // Vacuum ISP
double mdot   = engine.MassFlowRate(pressureRatio: 1); // kg/s at sea level
```

| Property | Description |
|----------|-------------|
| `IspSeaLevel` | Specific impulse at sea level (s) |
| `IspVacuum` | Specific impulse in vacuum (s) |
| `MaxThrustVacuum` | Maximum thrust at vacuum (N) |
| `MinThrottle` | Minimum throttle fraction |
| `MaxGimbalAngle` | Maximum gimbal deflection (rad) |
| `Throttle` | Current throttle setting |
| `IsActive` | Engine active state |

### ThrustCurve

Defines time-varying thrust profiles for solid boosters, throttle buckets, and custom profiles.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

// Constant thrust for 150s
var constant = ThrustCurve.Constant(burnDuration: 150);

// Solid booster with ramp-up and tail-off
var srb = ThrustCurve.SolidBooster(burnDuration: 120, rampFraction: 0.05, tailFraction: 0.1);

// MaxQ throttle bucket (reduce thrust between 60–80s)
var bucket = ThrustCurve.ThrottleBucket(burnDuration: 180,
    bucketStart: 60, bucketEnd: 80, bucketLevel: 0.7);

// Evaluate at any time
double fraction = srb.Evaluate(time: 10.0);  // 0.0–1.0
```

### PropellantTank

Tracks fuel and oxidizer consumption with mixture-ratio-preserving depletion.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

var tank = new PropellantTank(fuelMass: 25000, oxidizerMass: 60000);

double consumed = tank.Consume(massFlowRate: 250, dt: 1.0);
double remaining = tank.PropellantMass;
double fraction = tank.ConsumedFraction;
bool empty = tank.IsEmpty;

// Remaining ΔV via Tsiolkovsky equation
double dv = tank.RemainingDeltaV(dryMass: 5000, isp: 311);

tank.Reset(); // Refill for re-simulation
```

### VariableMassDynamics

Static methods for rocket rigid-body dynamics with variable mass, including thrust torque.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

Vector accel = VariableMassDynamics.ComputeAcceleration(
    totalMass: 500000, thrustForce: thrustVec, externalForces: gravity);

Vector angAccel = VariableMassDynamics.ComputeAngularAcceleration(
    inertiaTensor: inertia, torque: totalTorque, angularVelocity: omega);

Vector thrustTorque = VariableMassDynamics.ComputeThrustTorque(
    cgToEngine: armVector, thrustBody: thrustInBodyFrame);
```

### CenterOfMassTracker

Computes center of mass position as propellant drains, enabling correct moment-arm calculations.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

var tracker = new CenterOfMassTracker(
    dryMass: 5000, maxPropellantMass: 85000,
    cgFull: 20.0, cgEmpty: 15.0);

double cg = tracker.ComputeCg(currentPropellantMass: 42000);
double arm = tracker.CgToEngineDistance(currentPropellantMass: 42000, enginePosition: 0);
```

### GravityTurnGuidance

Open-loop guidance for the atmospheric ascent phase. Implements vertical ascent → pitch kick → gravity turn with azimuth targeting for desired orbital inclination.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

var guidance = new GravityTurnGuidance
{
    KickAltitude = 1000,
    KickAngle = 5.0 * Math.PI / 180,
    KickDuration = 10.0,
    TargetInclination = 28.5 * Math.PI / 180,
    LaunchLatitude = 28.6 * Math.PI / 180
};

double azimuth = guidance.LaunchAzimuth(); // Required heading

Quaternion cmd = guidance.ComputeAttitude(altitude: 5000, velocity: vel, time: 60);
// Phase progresses: VerticalAscent → PitchOver → GravityTurn
```

| Phase | Description |
|-------|-------------|
| `VerticalAscent` | Straight up until kick altitude |
| `PitchOver` | Gradual pitch kick over kick duration |
| `GravityTurn` | Velocity-aligned (prograde) |

### PEGGuidance (Powered Explicit Guidance)

Closed-loop terminal guidance that iteratively solves for the optimal thrust direction to achieve target orbital parameters. Based on the Space Shuttle's PEG algorithm.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

var peg = new PEGGuidance
{
    TargetSemiMajorAxis = 6571000,    // 200 km LEO
    TargetEccentricity = 0.0,
    TargetInclination = 28.5 * Math.PI / 180,
    Mu = 3.986004418e14
};

Vector thrustDir = peg.ComputeThrustDirection(
    position: posECI, velocity: velECI,
    thrustAccel: 30.0, exhaustVelocity: 3000);

bool converged = peg.HasConverged;
double tgo = peg.TimeToCutoff;
double dvRemaining = peg.VelocityToGain;
```

### AttitudeController

Quaternion-feedback PID controller with rate limiting and integral windup protection.

```csharp
using CSharpNumerics.Physics.Mechanics.Propulsion;

var controller = new AttitudeController
{
    Kp = 2.0, Kd = 1.0, Ki = 0.01,
    MaxRate = 5.0 * Math.PI / 180,
    MaxTorque = 50000
};

Vector torque = controller.ComputeTorque(
    currentAttitude: currentQuat,
    commandedAttitude: targetQuat,
    angularRate: omega,
    dt: 0.02);

Vector desiredRate = controller.ComputeDesiredRate(currentQuat, targetQuat);
```

---

## 🌍 Earth & Gravity Models

The `Physics.Mechanics` namespace provides WGS-84 Earth geometry and a J2-perturbed gravity model for accurate orbital and launch simulations.

### EarthModel (WGS-84)

```csharp
using CSharpNumerics.Physics.Mechanics;

double R = EarthModel.SemiMajorAxis;        // 6 378 137 m
double b = EarthModel.SemiMinorAxis;        // 6 356 752.3 m
double omega = EarthModel.RotationRate;     // 7.2921150e-5 rad/s
double mu = EarthModel.GM;                  // 3.986004418e14 m³/s²
double j2 = EarthModel.J2;                  // 1.08263e-3
double g0 = EarthModel.G0;                  // 9.80665 m/s²

// Geometry
double N = EarthModel.PrimeVerticalRadius(latitude: 0.5);  // Radius of curvature
double v = EarthModel.SurfaceVelocity(latitude: 0.5);      // m/s at surface
```

### GravityModel (J2 Zonal Harmonic)

```csharp
using CSharpNumerics.Physics.Mechanics;

// Full J2 gravity acceleration (includes oblateness)
Vector gJ2 = GravityModel.Acceleration(positionECI);

// Point-mass only (for comparison)
Vector gPM = GravityModel.PointMassAcceleration(positionECI);

// Scalar magnitude at radius and latitude
double g = GravityModel.Magnitude(radius: 6571000, latitude: 0.5);
```

### CoriolisForce

Non-inertial frame corrections for ECEF-referenced simulations.

```csharp
using CSharpNumerics.Physics.Mechanics;

Vector coriolis = CoriolisForce.Acceleration(velocityECEF);         // -2(ω × v)
Vector centrifugal = CoriolisForce.CentrifugalAcceleration(posECEF); // ω²(x, y, 0)
Vector total = CoriolisForce.TotalFictitiousAcceleration(posECEF, velECEF);
```

---

## 🛰️ Orbital Mechanics

The `Physics.OrbitalMechanics` namespace provides Keplerian elements, state-vector conversions, orbital propagation with J2 and drag perturbations, and impulsive maneuver planning.

### OrbitalElements

```csharp
using CSharpNumerics.Physics.OrbitalMechanics;

var elements = new OrbitalElements(
    a: 6771000,           // Semi-major axis (m)
    e: 0.001,             // Eccentricity
    i: 0.9,              // Inclination (rad)
    raan: 1.2,           // Right ascension of ascending node (rad)
    omega: 0.5,          // Argument of periapsis (rad)
    nu: 0.0,             // True anomaly (rad)
    mu: 3.986004418e14);  // Gravitational parameter

double period = elements.Period;              // Orbital period (s)
double rp = elements.Periapsis;              // Periapsis radius (m)
double ra = elements.Apoapsis;               // Apoapsis radius (m)
double v = elements.VelocityAtTrueAnomaly;   // Speed at current ν
bool bound = elements.IsBound;               // Negative energy?
```

| Property | Description |
|----------|-------------|
| `Period` | Orbital period ($2\pi\sqrt{a^3/\mu}$) |
| `MeanMotion` | Mean angular rate (rad/s) |
| `Periapsis` / `Apoapsis` | Apsidal radii (m) |
| `SpecificEnergy` | $-\mu / 2a$ |
| `SpecificAngularMomentum` | $h = \sqrt{\mu a(1-e^2)}$ |
| `IsCircular` | $e < 10^{-6}$ |

### State ↔ Elements Conversion

```csharp
using CSharpNumerics.Physics.OrbitalMechanics;

// State → Keplerian elements
OrbitalElements el = StateToElements.FromStateVector(position, velocity, mu);

// Elements → Cartesian state
var (pos, vel) = ElementsToState.ToStateVector(el);
```

### OrbitalPropagator

Velocity-Verlet integrator with optional J2 oblateness and atmospheric drag perturbations.

```csharp
using CSharpNumerics.Physics.OrbitalMechanics;

var prop = new OrbitalPropagator
{
    Mu = 3.986004418e14,
    IncludeJ2 = true,
    IncludeAtmosphericDrag = true,
    DragAltitudeThreshold = 600000,
    BallisticCoefficient = 50.0
};

prop.SetState(position, velocity);
prop.Propagate(duration: 5400, dt: 1.0);

Vector finalPos = prop.Position;
Vector finalVel = prop.Velocity;
OrbitalElements finalElements = prop.Elements;
```

### HohmannTransfer

Minimum-energy two-impulse transfer between circular orbits.

```csharp
using CSharpNumerics.Physics.OrbitalMechanics;

// Transfer from 200 km to 400 km orbit
var result = HohmannTransfer.FromAltitudes(alt1: 200000, alt2: 400000);

double dv1 = result.DeltaV1;          // First burn (m/s)
double dv2 = result.DeltaV2;          // Second burn (m/s)
double total = result.TotalDeltaV;     // Total ΔV
double tof = result.TransferTime;      // Coast duration (s)
```

### OrbitalManeuver & Circularization

```csharp
using CSharpNumerics.Physics.OrbitalMechanics;

// Prograde burn at specific time
var burn = OrbitalManeuver.Prograde(deltaV: 50, currentVelocity: vel, time: 1000);
Vector newVel = burn.Apply(vel);

// Circularization ΔV
double dvCirc = CircularizationBurn.Compute(position, velocity, mu);
double dvApo = CircularizationBurn.AtApoapsis(elements);
double dvPeri = CircularizationBurn.AtPeriapsis(elements);
```

---

## 🌡️ Compressible Aerodynamics

The `Physics.FluidDynamics.Aerodynamics` namespace includes compressible-flow models for high-speed vehicles — Mach-dependent drag, regime classification, and aerothermal heating.

### CompressibleDragModel

Mach-dependent drag coefficient with transonic rise and supersonic decline.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

var drag = new CompressibleDragModel(cd0: 0.3, cdPeak: 0.8, machPeak: 1.1);
double cd = drag.Evaluate(mach: 0.9);  // Transonic regime

// Presets
var rocket = CompressibleDragModel.SlenderRocket();   // Low Cd0
var capsule = CompressibleDragModel.BluntCapsule();   // High Cd0
```

### MachNumber

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

double mach = MachNumber.Compute(velocity: 340, altitude: 0);  // ~1.0
bool sub = MachNumber.IsSubsonic(mach);      // M < 0.8
bool trans = MachNumber.IsTransonic(mach);   // 0.8–1.2
bool sup = MachNumber.IsSupersonic(mach);    // 1.2–5.0
bool hyp = MachNumber.IsHypersonic(mach);    // M > 5.0
```

### HeatFlux (Sutton-Graves)

Stagnation-point convective heating for atmospheric re-entry and high-speed ascent.

```csharp
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

// Direct computation from density and velocity
double q = HeatFlux.StagnationPoint(
    density: 0.001, velocity: 7000, noseRadius: 1.0);  // W/m²

// From altitude and velocity (uses ISA atmosphere)
double qAlt = HeatFlux.FromAltitudeAndVelocity(
    altitude: 60000, velocity: 6000, noseRadius: 0.5);

double h0 = HeatFlux.StagnationEnthalpy(velocity: 7000);  // J/kg
```