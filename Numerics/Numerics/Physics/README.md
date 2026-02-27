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

This enables polymorphic use — the same analysis code works for `SimpleHarmonicOscillator`, `DampedOscillator`, `DrivenOscillator`, and `CoupledOscillators` (coming in future phases).

---

