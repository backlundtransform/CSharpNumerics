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

---

## 🌊 Fluid Extensions

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