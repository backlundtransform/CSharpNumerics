## 🛞 Kinematics Extensions

The `KinematicsExtensions` class provides a set of extension methods for performing common kinematic calculations in both **scalar** and **vector** form. It covers free fall, constant velocity, constant acceleration, time-independent SUVAT equations, and circular motion.

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

---

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


