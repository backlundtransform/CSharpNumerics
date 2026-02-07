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

Compute centripetal acceleration:

```csharp
double a = 3.0.CentripetalAcceleration(radius: 2); // scalar

Vector velocity = new Vector(2, 0, 0);
Vector radius = new Vector(0, 3, 0);
Vector ac = velocity.CentripetalAcceleration(radius); // vector, towards center
```

---


