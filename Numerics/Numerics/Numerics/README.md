
## 📘 Numeric Extensions

**Factorial**

```csharp
int result = 5.Factorial(); // 120
```

**Root Finding (Newton–Raphson)**

```csharp
Func<double, double> func = x => Math.Pow(x, 2) - 4;
double root = func.NewtonRaphson(); // 2
```

**Number Theory**

```csharp
bool prime = 79.IsPrime();                  // true
int[] factors = 78.GetPrimeFactors();       // [2, 3, 13]
bool happy = 19.IsHappy();                  // true
bool perfect = 6.IsPerfectNumber();         // true
int decimals = 0.01.GetDecimalPlaces();     // 2
```

**Trigonometry**

```csharp
double rad = 180.0.DegreeToRadians(); // π
```

---

## lim Limits

```csharp
Func<double, double> f = x => Math.Sin(x) / x;

double left  = f.LeftLimit(0);   // f(0 - ε)
double right = f.RightLimit(0);  // f(0 + ε)
double limit = f.Limit(0);       // two-sided: returns value if left == right, else NaN
```

---

## Δ Derivative

```csharp
Func<double, double> f = x => Math.Pow(x, 2);
Func<double, double> g = x => 4 * x - 3;
var result = f.Derivate(g, 1);
```

Supports **Chain**, **Product**, and **Quotient** rules via:

```csharp
var result = f.Derivate(g, Numerics.Enums.DerivateOperator.Product);
```

Multiple variables:

```csharp
Func<double[], double> func = vars => vars[0] * vars[1];
var dfdx = func.Derivate(new double[] { 2, 3 }, index: 0);
```

Or with vectors:

```csharp
Func<Vector, double> func = v => v.x * v.y;
var dfdx = func.Derivate(new Vector(2, 3, 0), Cartesian.X);
```

Derivate series:

```csharp
Func<double, double> displacement = t => 9.81 * Math.Pow(t, 2) / 2;
var velocity = displacement.GetSeries(0, 10, 1000).Derivate();
```

---

## ∫ Integrals

**Trapezoidal Rule**

```csharp
Func<double, double> f = x => Math.Sin(x);
double integral = f.Integrate(0, Math.PI);
```

**Simpson's 1/3 Rule** — O(h⁴) accuracy

```csharp
double result = f.IntegrateSimpson(0, Math.PI, subintervals: 200);
```

**Simpson's 3/8 Rule** — cubic interpolation, O(h⁴)

```csharp
double result = f.IntegrateSimpson38(0, Math.PI, subintervals: 999);
```

**Gauss-Legendre Quadrature** — 5-point per subinterval, very high accuracy

```csharp
double result = f.IntegrateGaussLegendre(0, Math.PI, subintervals: 10);
```

**Romberg Integration** — Richardson extrapolation for rapid convergence

```csharp
double result = f.IntegrateRomberg(0, Math.PI, maxLevel: 10);
```

**Adaptive Simpson** — recursive subdivision with error control

```csharp
double result = f.IntegrateAdaptive(0, Math.PI, tolerance: 1e-12);
```

**Monte Carlo Integration** — double and triple integrals

```csharp
Func<(double x, double y), double> func2d = p => p.x * p.y;
double result2d = func2d.Integrate((0, 1), (0, 1));

Func<Vector, double> func3d = v => v.x + v.y + v.z;
double result3d = func3d.Integrate(new Vector(-2, -2, -2), new Vector(2, 2, 2));
```

**Series / Time Series Integration**

```csharp
List<TimeSerie> ts = ...;
double total = ts.Integrate();
double bounded = ts.Integrate(startDate, endDate);
```

---

## 🧩 Complex Numbers

```csharp
var a = new ComplexNumber(3, 2);
var b = new ComplexNumber(5, 3);

var sum = a + b;
var product = a * b;
var power = a.Pow(2); // 5 + 12i
```

Exponential:

```csharp
new ComplexNumber(0, Math.PI).Exponential(); // -1
```

---

## 🧭 Vector

```csharp
var a = new Vector(5, 3, 0);
var b = new Vector(2, 6, 0);

var dot = a.Dot(b);
var cross = a.Cross(b);
```

From spherical coordinates:

```csharp
var v = Vector.FromSphericalCoordinates(radius, inclination, azimuth);
```

Vector of any length:

```csharp
double[] ydata =
{
    1,3,5,7,9,11,13,15,17,19
};
var y = new VectorN(ydata);
```

---

## 🔄 Quaternion

Quaternions generalize complex numbers to 4 dimensions: `q = w + xi + yj + zk`. They are the standard representation for 3D rotations — compact (4 doubles vs 9 for a matrix), numerically stable, and free of gimbal lock.

**Algebraic hierarchy:** `ℝ (double) ⊂ ℂ (ComplexNumber) ⊂ ℍ (Quaternion)`

```csharp
// Create from axis + angle
var q = Quaternion.FromAxisAngle(new Vector(0, 0, 1), Math.PI / 2); // 90° about Z

// Rotate a vector: v' = q·v·q*
var rotated = q.Rotate(new Vector(1, 0, 0)); // → (0, 1, 0)

// Compose rotations via multiplication (non-commutative)
var q2 = Quaternion.FromAxisAngle(new Vector(0, 1, 0), Math.PI / 3);
var combined = q2 * q; // q first, then q2

// Convert to/from 3×3 rotation matrix
Matrix m = q.ToMatrix();
Quaternion recovered = Quaternion.FromMatrix(m);

// Euler angles (ZYX convention)
var q3 = Quaternion.FromEulerAngles(roll: 0.3, pitch: 0.5, yaw: 0.7);
var (roll, pitch, yaw) = q3.ToEulerAngles();

// Axis-angle round-trip
var (axis, angle) = q.ToAxisAngle();
```

**Interpolation:**

```csharp
// Spherical linear interpolation (constant angular velocity)
var halfway = Quaternion.Slerp(qStart, qEnd, 0.5);

// Normalized linear interpolation (cheaper, good for small angles)
var approx = Quaternion.Lerp(qStart, qEnd, 0.5);
```

**Integration (for physics simulation):**

```csharp
// Integrate orientation with angular velocity ω over time step dt
var q = Quaternion.Identity;
var omega = new Vector(0, 0, 1); // 1 rad/s about Z
q = Quaternion.IntegrateOrientation(q, omega, dt: 0.001);
```

**Bridge from ComplexNumber:**

```csharp
// Embed ℂ into ℍ — preserves complex multiplication
var c = new ComplexNumber(3, 2);
var q = Quaternion.FromComplexNumber(c); // (3, 2, 0, 0)
```

---

## 🧮 Matrix

```csharp
var A = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });
var transpose = A.Transpose();
var det = A.Determinant();
var inv = A.Inverse();

Arithmetic:

var B = new Matrix(new double[,] { { 2, 5, 1 }, { 4, 3, 7 } });
var sum = A + B;
var product = A * B;

With vector:

var x = new Vector(2, 1, 3);
var y = A * x;
```

---

## 🔢 Complex Linear Algebra

`ComplexVector`, `ComplexVectorN`, and `ComplexMatrix` mirror the real-valued types with full complex number support. Existing real types convert implicitly — no API breakage.

### ComplexVector (3D)

```csharp
var a = new ComplexVector(
    new ComplexNumber(1, 2),
    new ComplexNumber(3, 0),
    new ComplexNumber(0, -1));

var b = new ComplexVector(
    new ComplexNumber(2, 1),
    new ComplexNumber(0, 3),
    new ComplexNumber(1, 1));

var sum = a + b;
var dot = a.Dot(b);                  // standard complex dot product
var hermitian = a.HermitianDot(b);   // ⟨a,b⟩ = Σ conj(aᵢ)·bᵢ
var cross = a.Cross(b);

double mag = a.GetMagnitude();       // hermitian norm: √(Σ|xᵢ|²)
var conj = a.GetConjugate();
var unit = a.GetUnitVector();

// Implicit from real Vector
Vector v = new Vector(1, 2, 3);
ComplexVector cv = v;                // imaginary parts are zero
```

### ComplexVectorN (N-dimensional)

```csharp
var v = new ComplexVectorN(new ComplexNumber[]
{
    new ComplexNumber(1, 2),
    new ComplexNumber(3, -1),
    new ComplexNumber(0, 4)
});

var dot = v.Dot(v);                  // Σ vᵢ·vᵢ
var hermitian = v.HermitianDot(v);   // Σ conj(vᵢ)·vᵢ  (always real for self)
var hadamard = v.Hadamard(v);
double norm = v.Norm();
var unit = v.Normalize();
var conj = v.GetConjugate();

// Implicit from real VectorN
VectorN real = new VectorN(new double[] { 1, 2, 3 });
ComplexVectorN complex = real;
```

### ComplexMatrix (NxM)

```csharp
var A = new ComplexMatrix(new ComplexNumber[,]
{
    { new ComplexNumber(1, 0), new ComplexNumber(0, 1) },
    { new ComplexNumber(0, -1), new ComplexNumber(1, 0) }
});

var transpose = A.Transpose();
var dagger = A.ConjugateTranspose();  // hermitian adjoint (A†)
var det = A.Determinant();            // returns ComplexNumber
var inv = A.Inverse();

// Arithmetic
var B = new ComplexNumber(2, 0) * A;
var C = A * A;

// Implicit from real Matrix
Matrix real = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });
ComplexMatrix complex = real;
```

---

## 📦 Tensor (multi-dimensionell)

```csharp

var tensor = new Tensor(2, 3);

tensor[0, 0] = 1;
tensor[0, 1] = 2;
tensor[0, 2] = 3;
tensor[1, 0] = 4;
tensor[1, 1] = 5;
tensor[1, 2] = 6;

tensor.Fill(10);

var tensorB = new Tensor(2, 3);
tensorB.Fill(5);

var sum = tensor + tensorB;
var diff = tensor - tensorB;
var prod = tensor * tensorB;
var div = tensor / tensorB;


var tensor1D = new Tensor(3);
tensor1D.Values[0] = 1;
tensor1D.Values[1] = 2;
tensor1D.Values[2] = 3;

var tensor1D2 = new Tensor(3);
tensor1D2.Values[0] = 4;
tensor1D2.Values[1] = 5;
tensor1D2.Values[2] = 6;

double dot = tensor1D.Dot(tensor1D2); // 1*4 + 2*5 + 3*6 = 32
```

---

## 📐 Scalar Field

A `ScalarField` wraps `Func<Vector, double>` and provides instance methods for the standard vector-calculus operations that are otherwise spread across extension classes. It converts implicitly to `Func<Vector, double>`, so it can be passed to any existing API.

```csharp
var V = new ScalarField(r => Math.Pow(r.x, 2) * Math.Pow(r.y, 3));

// Evaluate
double val = V.Evaluate(new Vector(1, -2, 0));
double val2 = V.Evaluate((1, -2, 0));

// Differential operators
Vector grad = V.Gradient((1, -2, 0));          // ∇f
double lap = V.Laplacian((1, -2, 0));          // ∇²f
double dfdx = V.Derivate(new Vector(1, -2, 0), Cartesian.x);

// Gradient as a full VectorField — ∇f(r) lazily evaluated
VectorField gradField = V.GradientField();
Vector curlOfGrad = gradField.Curl((1, -2, 0)); // ≈ 0 (vector identity)

// Integration (3D Monte Carlo)
double integral = V.Integrate(new Vector(-1, -1, -1), new Vector(1, 1, 1));
```

**Arithmetic** — build composite fields from simpler ones:

```csharp
var f = new ScalarField(r => r.x * r.x);
var g = new ScalarField(r => r.y * r.y);

ScalarField sum = f + g;       // f(r) + g(r)
ScalarField diff = f - g;      // f(r) - g(r)
ScalarField prod = f * g;      // f(r) · g(r)
ScalarField scaled = 3.0 * f;  // 3 · f(r)
ScalarField neg = -f;          // −f(r)
```

**Interop** — implicit conversion to `Func<Vector, double>`:

```csharp
var V = new ScalarField(r => r.x + r.y + r.z);

// Pass directly to any method that takes Func<Vector, double>
Func<Vector, double> func = V;
```

---

## 🌐 Vector Field

**Gradient**

```csharp
Func<Vector, double> f = p => Math.Pow(p.x, 2) * Math.Pow(p.y, 3);
var grad = f.Gradient((1, -2, 0));
```

**Laplacian**

```csharp
Func<Vector, double> f = p => p.x * p.x + p.y * p.y + p.z * p.z;
double laplacian = f.Laplacian((1, 2, 3)); // ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²
```

**Divergence**

```csharp
var field = new VectorField(p => Math.Sin(p.x * p.y),
                            p => Math.Cos(p.x * p.y),
                            p => Math.Exp(p.z));

double div = field.Divergence((1, 2, 2));
```

**Curl**

```csharp
var field = new VectorField(p => p.y, p => -p.x, p => 0);
var curl = field.Curl((1, 4, 2));
```

**Grid Evaluation**

`EvaluateRange` walks a diagonal line. For proper 2D arrow-plot rendering, use `EvaluateGrid2D`:

```csharp
var field = new VectorField(p => -p.y, p => p.x, p => 0);

// Structured nx × ny grid of vectors
var grid = field.EvaluateGrid2D(
    xmin: -2, xmax: 2,
    ymin: -2, ymax: 2,
    nx: 20, ny: 20);

// Returns IDictionary<Vector, Vector> — position → field value
foreach (var (pos, val) in grid)
    Console.WriteLine($"({pos.x:F1}, {pos.y:F1}) → ({val.x:F2}, {val.y:F2})");
```

---

## 🧊 Tensor Field

A `TensorField` maps each spatial point to a 3×3 `Matrix` — the rank-2 completion of the field hierarchy:

| Rank | Type | Value at point |
|---|---|---|
| 0 | `ScalarField` | `double` |
| 1 | `VectorField` | `Vector` |
| 2 | `TensorField` | `Matrix` (3×3) |

```csharp
// Stress tensor that varies with position
var T = new TensorField(r => new Matrix(new double[,]
{
    { r.x, 0,   0   },
    { 0,   r.y, 0   },
    { 0,   0,   r.z }
}));

Matrix val = T.Evaluate((1, 2, 3));
```

**Component access:**

```csharp
ScalarField Txy = T.Component(0, 1);          // Tᵢⱼ as ScalarField
VectorField row0 = T.Row(0);                   // (T₀₁, T₀₂, T₀₃)
```

**Divergence** — ∇·T → VectorField, where (∇·T)ᵢ = Σⱼ ∂Tᵢⱼ/∂xⱼ:

```csharp
VectorField forceDensity = T.Divergence();
Vector f = forceDensity.Curl((1, 2, 3));       // chain with existing operators
```

**Trace and contractions:**

```csharp
ScalarField tr = T.Trace();                    // Σᵢ Tᵢᵢ
VectorField Tv = T.Contract(velocityField);    // T·v
ScalarField energy = T.DoubleContract(strainField); // T : S = Σᵢⱼ TᵢⱼSᵢⱼ
```

**Jacobian** — gradient of a vector field → TensorField:

```csharp
var F = new VectorField(r => r.x * r.y, r => r.z, r => r.x);
TensorField J = TensorField.FromJacobian(F);  // (∇F)ᵢⱼ = ∂Fᵢ/∂xⱼ
```

**Arithmetic:**

```csharp
TensorField sum = T + T;
TensorField scaled = 2.0 * T;
TensorField transposed = T.Transpose();
TensorField weighted = scalarField * T;       // ScalarField × TensorField
```

---

## ⚙️ Transform

**FFT / DFT**

```csharp
Func<double, double> f = t => Math.Exp(-t * t / 0.02);
var freq = f.FastFourierTransform(-0.5, 0.5, 100)
             .ToFrequencyResolution(100);

// DFT from real or complex samples
var dft = samples.DiscreteFourierTransform();
```

**Inverse FFT / DFT**

```csharp
var timeDomain = freqDomain.InverseFastFourierTransform();
var timeDomain2 = freqDomain.InverseDiscreteFourierTransform();
```

**Laplace Transform**

```csharp
double result = f.LaplaceTransform(2.0);
double inverse = F.InverseLaplaceTransform(1.0);
```

**Low-Pass Filter**

```csharp
var filtered = input.LowPassFilter(output, alpha: 0.25);
```

---

## 📐 Differential Equations

**Scalar ODE — Runge–Kutta (RK4)**

```csharp
Func<(double t, double y), double> f = v => Math.Tan(v.y) + 1;
var result = f.RungeKutta(1, 1.1, 0.025, 1);
```

**Scalar ODE — Euler Method**

```csharp
var result = f.EulerMetod(min: 0, max: 1, stepSize: 0.01, yInitial: 1);
```

**Scalar ODE — Trapezoidal Rule**

```csharp
var result = f.TrapezoidalRule(min: 0, max: 1, stepSize: 0.01, yInitial: 1);
```

**Custom Butcher Tableau**

```csharp
var result = f.RungeKutta(min, max, stepSize, yInitial,
    rungeKuttaMatrix, weights, nodes);
```

**Vector ODE (3D)** — solve y' = f(t, y) where y is a Vector

```csharp
// Exponential decay in 3D
Func<(double t, Vector y), Vector> decay = v => -1.0 * v.y;
Vector result = decay.RungeKutta(0, 1, 0.001, new Vector(1, 2, 3));

// Euler method
Vector result2 = decay.EulerMethod(0, 1, 0.001, new Vector(1, 2, 3));

// Full trajectory
var trajectory = decay.RungeKuttaTrajectory(0, 1, 0.01, new Vector(1, 2, 3));
foreach (var (t, y) in trajectory) { /* ... */ }
```

**System ODE (VectorN)** — arbitrary-dimensional systems for dynamics

```csharp
// Free fall: state = [x, y, z, vx, vy, vz]
double g = 9.80665;
Func<(double t, VectorN y), VectorN> dynamics = v =>
    new VectorN([v.y[3], v.y[4], v.y[5], 0, 0, -g]);

var y0 = new VectorN([0, 0, 0, 10, 0, 10]);
VectorN result = dynamics.RungeKutta(0, 2, 0.001, y0);

// Euler method
VectorN result2 = dynamics.EulerMethod(0, 2, 0.001, y0);

// Full trajectory
var orbit = dynamics.RungeKuttaTrajectory(0, period, 1.0, y0);
double energy = orbit.Last().y.Dot(orbit.Last().y); // use VectorN operations
```

`double[]` convenience overloads are also available — they delegate to `VectorN` internally.

**Semi-Implicit (Symplectic) Euler** — stable for oscillatory systems

State is split as `[positions | velocities]`; velocities are updated first, then positions use the new velocities.

```csharp
// Simple harmonic oscillator: state = [x, v], dy/dt = [v, -x]
Func<(double t, VectorN y), VectorN> sho = v =>
    new VectorN([v.y[1], -v.y[0]]);

var y0 = new VectorN([1, 0]);
VectorN result = sho.SemiImplicitEuler(0, 100, 0.01, y0);

// Full trajectory
var traj = sho.SemiImplicitEulerTrajectory(0, 100, 0.01, y0);
```

**Velocity Verlet** — O(dt²) symplectic, excellent energy conservation

For second-order ODEs. State is `[positions | velocities]`; the acceleration function returns only accelerations (half the state length).

```csharp
// Spring: acceleration = -x (depends on position only)
Func<(double t, VectorN y), VectorN> accel = v =>
    new VectorN([-v.y[0]]);

var y0 = new VectorN([1, 0]); // [position, velocity]
VectorN result = accel.VelocityVerlet(0, 100, 0.01, y0);
```

---

## ⏱️ Time Stepping (Class-Based ODE Integration)

The extension methods above (`RungeKutta`, `EulerMethod`, `VelocityVerlet`, etc.) are functional and convenient for simple problems. The `ITimeStepper` interface provides a **class-based alternative** for advanced scenarios: adaptive step control, trajectory recording, per-step callbacks, and diagnostics.

### Interface

```csharp
public interface ITimeStepper
{
    string Name { get; }
    TimeStepResult Solve(
        Func<double, VectorN, VectorN> rhs,
        double t0, double tEnd, VectorN y0, double dt,
        bool recordTrajectory = false,
        Action<double, VectorN> onStep = null);
}
```

### Available Steppers

| Stepper | Class | Order | Step Size | Best For |
|---------|-------|-------|-----------|----------|
| Euler | `EulerStepper` | O(h) | Fixed | Quick prototyping, reference solutions |
| RK4 | `RK4Stepper` | O(h⁴) | Fixed | General-purpose, non-stiff problems |
| Dormand-Prince 4(5) | `AdaptiveRK45Stepper` | O(h⁴)/O(h⁵) | **Adaptive** | Automatic accuracy, varying timescales |
| Velocity Verlet | `VelocityVerletStepper` | O(h²) | Fixed | Symplectic, N-body, energy conservation |

### Result Object

Every stepper returns a `TimeStepResult`:

| Property | Description |
|----------|-------------|
| `T` | Final time reached |
| `Y` | Final state vector |
| `Trajectory` | Full `(t, y)` history (if `recordTrajectory = true`) |
| `Steps` | Total accepted steps |
| `RejectedSteps` | Rejected steps (adaptive only) |
| `FunctionEvaluations` | Total rhs evaluations |
| `LastError` | Estimated local error at final step (adaptive only) |
| `LastStepSize` | Actual step size used at final step |

### ➡️ Heat Equation (Method of Lines + RK4)

```csharp
var grid = new Grid2D(50, 50, 0.1);
VectorN u = grid.Initialize((x, y) => Math.Exp(-(x * x + y * y)));

Func<double, VectorN, VectorN> rhs = (t, state) =>
    GridOperators.Laplacian2D(state, grid, BoundaryCondition.Dirichlet);

var stepper = new RK4Stepper();
var result = stepper.Solve(rhs, 0, 1.0, u, dt: 0.0001, recordTrajectory: true);

double[,] finalField = grid.ToArray(result.Y);
Console.WriteLine($"Steps: {result.Steps}, Evaluations: {result.FunctionEvaluations}");
```

### ➡️ Adaptive Step Control (Dormand-Prince)

No manual `dt` tuning — the stepper finds the right step size automatically:

```csharp
Func<double, VectorN, VectorN> rhs = (t, state) =>
    GridOperators.Laplacian2D(state, grid, BoundaryCondition.Neumann);

var stepper = new AdaptiveRK45Stepper
{
    AbsoluteTolerance = 1e-8,
    RelativeTolerance = 1e-8
};

var result = stepper.Solve(rhs, 0, 10.0, u, dt: 0.01);

Console.WriteLine($"Steps: {result.Steps}, Rejected: {result.RejectedSteps}");
Console.WriteLine($"Final dt: {result.LastStepSize:E3}, Error: {result.LastError:E3}");
```

### ➡️ Spring System (Velocity Verlet)

Symplectic integration for oscillatory systems — state = `[positions | velocities]`:

```csharp
// Simple harmonic oscillator: x'' = -x
Func<double, VectorN, VectorN> sho = (t, y) =>
    new VectorN([y[1], -y[0]]);  // [velocity, acceleration]

var stepper = new VelocityVerletStepper();
var result = stepper.Solve(sho, 0, 100, new VectorN([1, 0]), dt: 0.01,
    recordTrajectory: true);

// Energy is conserved over long integrations
foreach (var (t, y) in result.Trajectory)
{
    double energy = 0.5 * (y[0] * y[0] + y[1] * y[1]);
    // energy ≈ 0.5 throughout
}
```

### ➡️ Per-Step Callback

Monitor or log intermediate states without storing the full trajectory:

```csharp
double maxTemp = 0;
var stepper = new RK4Stepper();
var result = stepper.Solve(rhs, 0, 10.0, u, dt: 0.001,
    onStep: (t, y) =>
    {
        double peak = y.Values.Max();
        if (peak > maxTemp) maxTemp = peak;
    });

Console.WriteLine($"Peak temperature: {maxTemp:F4}");
```

### ➡️ Swap Stepper at Runtime

All steppers share `ITimeStepper` — swap methods without changing application code:

```csharp
ITimeStepper stepper = needsAdaptive
    ? new AdaptiveRK45Stepper { AbsoluteTolerance = 1e-6 }
    : new RK4Stepper();

var result = stepper.Solve(rhs, 0, tEnd, y0, dt: 0.01);
```

**Key points:**

* Complements (does not replace) the functional extension methods
* `AdaptiveRK45Stepper` uses the Dormand-Prince 4(5) embedded pair — same method as MATLAB's `ode45`
* FSAL (first-same-as-last) optimization: only 6 new evaluations per accepted step
* `VelocityVerletStepper` is symplectic — ideal for Hamiltonian systems
* All steppers report diagnostics (step count, evaluations, error)
* `onStep` callback enables monitoring without memory overhead of full trajectory

// Full trajectory
var traj = accel.VelocityVerletTrajectory(0, 20, 0.01, y0);
foreach (var (t, y) in traj)
{
    double energy = 0.5 * (y[0] * y[0] + y[1] * y[1]); // bounded
}
```

`double[]` convenience overloads are available for both — they delegate to `VectorN` internally.

**Linear ODE System Solver** — eigenvalue-based analytical solution

```csharp
var solutions = matrix.OdeSolver(new List<double> { 1, 0, 0 });
double x_t = solutions[0](t);
double y_t = solutions[1](t);
```

---

## � Finite Difference (PDE via Method of Lines)

Discretize spatial domains and solve time-dependent PDEs by converting them to ODE systems. The spatial state lives in a `VectorN`, discrete operators produce the right-hand side, and the existing ODE solvers (RK4, Euler, Verlet) handle time integration.

### Grid2D

A uniform 2D grid that packs/unpacks between `double[,]` and flat `VectorN` (row-major).

```csharp
var grid = new Grid2D(nx: 50, ny: 50, dx: 0.1, dy: 0.1);

// Initialize from a function (cell centres)
var u0 = grid.Initialize((x, y) => Math.Exp(-(x * x + y * y)));

// Pack/unpack
double[,] field = grid.ToArray(u0);
VectorN v = grid.ToVector(field);

// Index helpers
int flat = grid.Index(ix: 10, iy: 20);
var (ix, iy) = grid.Index2D(flat);
```

### Boundary Conditions

All operators accept a `BoundaryCondition` parameter:

| Type | Description |
|------|-------------|
| `Dirichlet` | u = 0 outside domain |
| `Neumann` | ∂u/∂n = 0 (zero-flux, mirrors boundary cell) |
| `Periodic` | Domain wraps around |

### Discrete Operators

**Laplacian** — ∇²u (3-point stencil in 1D, 5-point in 2D)

```csharp
var lap1d = GridOperators.Laplacian1D(u, dx, BoundaryCondition.Dirichlet);
var lap2d = GridOperators.Laplacian2D(u, grid, BoundaryCondition.Neumann);
```

**Gradient** — ∂u/∂x (central differences)

```csharp
var grad1d = GridOperators.Gradient1D(u, dx, BoundaryCondition.Periodic);
var (dux, duy) = GridOperators.Gradient2D(u, grid);
```

**Divergence** — ∇·F = ∂Fx/∂x + ∂Fy/∂y

```csharp
var div = GridOperators.Divergence2D(fx, fy, grid);
```

**Advection** — v·∇u (first-order upwind for stability)

```csharp
var adv = GridOperators.Advection2D(u, vx, vy, grid);
```

**Biharmonic (4th derivative)** — d⁴u/dx⁴ (5-point central stencil, for Euler–Bernoulli beam equation)

```csharp
var d4u = GridOperators.Biharmonic1D(u, dx, BoundaryCondition.Neumann);
```

### Example: 2D Heat Equation

$$\frac{\partial u}{\partial t} = \alpha \nabla^2 u$$

```csharp
double alpha = 0.01;
var grid = new Grid2D(50, 50, 0.1);

// Hot spot in the centre
var u0 = grid.Initialize((x, y) =>
    (x > 2.0 && x < 3.0 && y > 2.0 && y < 3.0) ? 100.0 : 0.0);

// RHS: spatial operator → feeds into standard ODE solver
Func<(double t, VectorN y), VectorN> rhs = args =>
    alpha * GridOperators.Laplacian2D(args.y, grid, BoundaryCondition.Dirichlet);

// Solve with existing RK4 — no new solver needed
var trajectory = rhs.RungeKuttaTrajectory(0, 10.0, 0.005, u0);

// Each frame: trajectory[i].y → grid.ToArray() → 2D heatmap
double[,] finalState = grid.ToArray(trajectory.Last().y);
```

### Example: Advection-Diffusion

$$\frac{\partial u}{\partial t} = \alpha \nabla^2 u - \mathbf{v} \cdot \nabla u$$

```csharp
Func<(double t, VectorN y), VectorN> rhs = args =>
    alpha * GridOperators.Laplacian2D(args.y, grid)
    - GridOperators.Advection2D(args.y, vx, vy, grid);

var result = rhs.RungeKutta(0, 5.0, 0.001, u0);
```

### Poisson Solver (Gauss-Seidel)

Iterative solver for $\nabla^2 u = f$ on a `Grid2D` with arbitrary Dirichlet boundary masks:

```csharp
// Solve ∇²φ = rhs on a 50×50 grid with Dirichlet BCs
var grid = new Grid2D(50, 50, 0.02);
var rhs = grid.Zeros();                      // source term

// Mark boundary cells
bool[] mask = new bool[grid.Length];
double[] bcValues = new double[grid.Length];
for (int iy = 0; iy < 50; iy++)
{
    mask[grid.Index(0, iy)] = true;           // left = 0 V
    mask[grid.Index(49, iy)] = true;          // right = 100 V
    bcValues[grid.Index(49, iy)] = 100.0;
}

var (solution, iterations) = GridOperators.SolvePoisson2D(
    rhs, grid, mask, bcValues,
    tolerance: 1e-8,
    maxIterations: 20000);

// solution is a VectorN, unpack to 2D:
double[,] phi = grid.ToArray(solution);
```

### Grid3D — 3D Finite Difference

Structured 3D grid with row-major indexing (ix fastest, then iy, then iz):

```csharp
var grid = new Grid3D(nx: 20, ny: 20, nz: 10, dx: 0.5, dy: 0.5, dz: 0.5);

int flat = grid.Index(5, 10, 3);        // (ix, iy, iz) → flat index
var (ix, iy, iz) = grid.Index3D(flat);   // flat → (ix, iy, iz)

// Initialize from function
VectorN u = grid.Initialize((x, y, z) => x * x + y * y + z * z);

// Pack / unpack between flat VectorN and 3D array
double[,,] arr = grid.ToArray(u);
VectorN v = grid.ToVector(arr);
```

### 3D Discrete Operators

| Operator | Method | Description |
|----------|--------|-------------|
| Laplacian | `GridOperators3D.Laplacian3D(u, grid, bc)` | 7-point stencil |
| Gradient | `GridOperators3D.Gradient3D(u, grid, bc)` | Central differences |
| Divergence | `GridOperators3D.Divergence3D(fx, fy, fz, grid, bc)` | Divergence of vector field |
| Advection | `GridOperators3D.Advection3D(u, vx, vy, vz, grid, bc)` | First-order upwind |
| Poisson | `GridOperators3D.SolvePoisson3D(rhs, grid, mask, vals)` | Gauss-Seidel iterative |

```csharp
// Laplacian of a quadratic field (interior cells ≈ 6.0)
var u = grid.Initialize((x, y, z) => x * x + y * y + z * z);
var lap = GridOperators3D.Laplacian3D(u, grid, BoundaryCondition.Dirichlet);
```

---

## �📏 Linear Systems

```csharp
var result = A.LinearSystemSolver(b);
var eigenValues = A.EigenValues();
var eigenVector = A.EigenVector(eigenValue);
var dominant = A.DominantEigenVector();
```

**Gauss Elimination**

```csharp
var solution = matrix.GaussElimination(vector);
```

---
##  ✨ Interpolation

CSharpNumerics provides a rich interpolation toolkit — from simple piecewise methods to polynomial, spline, rational, trigonometric, and multivariate interpolation.

### Piecewise (two-point) interpolation

All piecewise methods are accessible via a unified dispatcher:

```csharp
double y = data.Interpolate(p => (p.Index, p.Value), 3.5, InterpolationType.Linear);
```

| Type | Enum | Description |
|------|------|-------------|
| Linear | `Linear` | $y = y_1 + (y_2 - y_1)\frac{x - x_1}{x_2 - x_1}$ |
| Log–Log | `Logarithmic` | Linear in log-space for both x and y |
| Lin–Log | `LinLog` | Linear in x, logarithmic in y |
| Log–Lin | `LogLin` | Logarithmic in x, linear in y |

### Polynomial interpolation

Passes a single polynomial of degree N−1 through all N data points.

```csharp
double[] x = { 0, 1, 2, 3 };
double[] y = { 1, 2, 0, 5 };
var poly = new PolynomialInterpolation(x, y);

double val  = poly.Evaluate(1.5);             // Lagrange basis form
double val2 = poly.EvaluateNewton(1.5);       // Newton divided-difference
var (val3, err) = poly.EvaluateNeville(1.5);  // Neville with error estimate

// Or via the extension method:
double val4 = data.Interpolate(p => (p.Index, p.Value), 1.5, InterpolationType.Polynomial);
```

### Cubic Spline interpolation

Piecewise cubic polynomials with C² continuity. Three boundary conditions:

| Boundary | Description |
|----------|-------------|
| `Natural` | $S''(x_0) = S''(x_n) = 0$ (free ends) |
| `Clamped` | First derivative specified at endpoints |
| `NotAKnot` | Third derivative continuous at second & second-to-last knot |

```csharp
double[] x = { 0, 1, 2, 3, 4 };
double[] y = { 0, 1, 0, 1, 0 };

var spline = new CubicSplineInterpolation(x, y);                       // Natural
var clamped = new CubicSplineInterpolation(x, y, SplineBoundary.Clamped, 1.0, -1.0);
var nak     = new CubicSplineInterpolation(x, y, SplineBoundary.NotAKnot);

double val   = spline.Evaluate(2.5);
double dydx  = spline.Derivative(2.5);         // first derivative
double d2y   = spline.SecondDerivative(2.5);    // curvature

// Or via extension method:
double val2 = data.Interpolate(p => (p.Index, p.Value), 2.5, InterpolationType.CubicSpline);
```

### Rational interpolation

Ratio of two polynomials — handles poles and near-singularities better than polynomials.

```csharp
double[] x = { 0, 1, 2, 3, 4 };
double[] y = { 1.0, 0.5, 0.333, 0.25, 0.2 };   // ≈ 1/(1+x)
var rat = new RationalInterpolation(x, y);

var (val, err) = rat.Evaluate(1.5);              // Bulirsch–Stoer + error estimate
double val2 = rat.EvaluateFloaterHormann(1.5);   // barycentric, guaranteed pole-free
```

### Trigonometric interpolation

Best for periodic functions. Builds a trigonometric polynomial from the data.

```csharp
int N = 16;
double[] x = new double[N], y = new double[N];
for (int i = 0; i < N; i++)
{
    x[i] = 2 * Math.PI * i / N;
    y[i] = Math.Sin(x[i]) + 0.5 * Math.Cos(2 * x[i]);
}

var trig = new TrigonometricInterpolation(x, y, period: 2 * Math.PI);
double val  = trig.Evaluate(Math.PI / 3);
double dydx = trig.Derivative(Math.PI / 3);

// Fourier coefficients
double[] a = trig.CosineCoefficients;   // a_0, a_1, ..., a_M
double[] b = trig.SineCoefficients;     // b_0, b_1, ..., b_M

// Or via extension method:
double val2 = data.Interpolate(p => (p.Index, p.Value), 1.5, InterpolationType.Trigonometric);
```

### Multivariate interpolation

For scattered data in multiple dimensions.

**Inverse Distance Weighting (IDW / Shepard)**

```csharp
double[][] points = {
    new[] { 0.0, 0.0 }, new[] { 1.0, 0.0 },
    new[] { 0.0, 1.0 }, new[] { 1.0, 1.0 }, new[] { 0.5, 0.5 }
};
double[] values = points.Select(p => p[0] * p[0] + p[1] * p[1]).ToArray();

var interp = new MultivariateInterpolation(points, values);
double val = interp.EvaluateIDW(new[] { 0.25, 0.75 }, power: 2);
```

**Radial Basis Functions (RBF)**

```csharp
double val = interp.EvaluateRBF(new[] { 0.25, 0.75 }, RbfKernel.Gaussian);
```

| Kernel | Formula |
|--------|---------|
| `Gaussian` | $\phi(r) = e^{-(r/\varepsilon)^2}$ |
| `Multiquadric` | $\phi(r) = \sqrt{1 + (r/\varepsilon)^2}$ |
| `InverseMultiquadric` | $\phi(r) = 1/\sqrt{1 + (r/\varepsilon)^2}$ |
| `ThinPlateSpline` | $\phi(r) = r^2 \ln(r)$ |
| `Cubic` | $\phi(r) = r^3$ |

**Bilinear / Trilinear (regular grids)**

```csharp
double val2d = MultivariateInterpolation.Bilinear(xGrid, yGrid, gridValues, xi, yi);
double val3d = MultivariateInterpolation.Trilinear(xGrid, yGrid, zGrid, gridValues, xi, yi, zi);
```

---

## 📊 Signal Processing

The `Numerics.SignalProcessing` namespace provides tools for spectral analysis and synthesis of periodic signals.

### FourierSeries

Real Fourier series analysis and synthesis for periodic functions:

$$f(t) \approx \frac{a_0}{2} + \sum_{n=1}^{N} \left[ a_n \cos(n\omega t) + b_n \sin(n\omega t) \right]$$

```csharp
using CSharpNumerics.Numerics.SignalProcessing;

var fs = new FourierSeries();

// Analyse a square wave
double period = 2 * Math.PI;
Func<double, double> square = t => (t % period) < period / 2 ? 1.0 : -1.0;
fs.Analyze(square, period, nTerms: 20, nSamples: 4096);

// Coefficients
double a0 = fs.A0;        // ≈ 0 (zero mean)
double b1 = fs.Bn[0];     // ≈ 4/π (first sine coefficient)
double b3 = fs.Bn[2];     // ≈ 4/(3π)

// Reconstruct signal
double val = fs.Synthesize(t: 0.5);

// Sweep across time — vary term count for Gibbs phenomenon demo
double[] tVals = Enumerable.Range(0, 500).Select(i => i * period / 500).ToArray();
List<Serie> synth5  = fs.SynthesizeRange(tVals, maxTerms: 5);
List<Serie> synth20 = fs.SynthesizeRange(tVals, maxTerms: 20);

// Power spectrum: |cₙ|² per harmonic
List<Serie> power = fs.PowerSpectrum();

// Parseval's theorem: frequency-domain energy = time-domain energy
double eParseval = fs.ParsevalEnergy();
double eTime     = fs.TimeDomainEnergy(square, nSamples: 4096);
// eParseval ≈ eTime
```

---

## 🎯 Optimization

The `Numerics.Optimization` sub-section provides single- and multi-objective optimisers, convergence strategies, and training utilities. ML models consume these — they never reimplement optimisation logic.

### Architecture

```
Optimization/
├── Interfaces/
│   ├── IOptimizer.cs           — Step(parameters, gradient) → updated parameters
│   ├── IObjectiveFunction.cs   — Evaluate(x), Gradient(x), Dimension
│   └── IConvergenceCriterion.cs— HasConverged(iteration, loss, gradientNorm)
├── SingleObjective/
│   ├── GradientDescent.cs      — Vanilla / Momentum / Nesterov / L2
│   ├── Adam.cs                 — Adam + AdamW (decoupled weight decay)
│   ├── CoordinateDescent.cs    — Cyclic CD with L1/L2 soft-thresholding
│   └── Minimizer.cs            — Ties IOptimizer + IConvergenceCriterion + IObjectiveFunction
├── MultiObjective/
│   ├── ParetoFront.cs          — Non-dominated sorting + crowding distance
│   └── NSGA2.cs                — NSGA-II evolutionary multi-objective optimiser
└── Strategies/
    ├── EarlyStopping.cs        — Patience-based training halt
    ├── LearningRateSchedule.cs — Constant / StepDecay / Exponential / InverseTime / Cosine
    └── MaxIterationsOrTolerance.cs — Simple IConvergenceCriterion
```

### Interfaces

| Interface | Key method | Purpose |
|-----------|-----------|---------|
| `IOptimizer` | `VectorN Step(VectorN parameters, VectorN gradient)` | Single gradient-based update step |
| `IObjectiveFunction` | `double Evaluate(double[] x)`, `double[] Gradient(double[] x)` | Objective + gradient for `Minimizer` |
| `IConvergenceCriterion` | `bool HasConverged(int iteration, double loss, double gradNorm)` | Termination check |

### Single-Objective Optimisers

**GradientDescent** — vanilla, momentum, Nesterov

```csharp
var gd = new GradientDescent(learningRate: 0.01, momentum: 0.9, nesterov: true);

VectorN w = new VectorN(new double[] { 0, 0, 0 });
VectorN grad = ComputeGradient(w);
w = gd.Step(w, grad);   // one update step
```

**Adam / AdamW** — adaptive moment estimation

```csharp
var adam = new Adam(learningRate: 0.001, beta1: 0.9, beta2: 0.999);
w = adam.Step(w, grad);

// AdamW (decoupled weight decay)
var adamw = new Adam(learningRate: 0.001, weightDecay: 0.01, decoupledWeightDecay: true);
```

**CoordinateDescent** — L1/L2 regularised least-squares

Minimises $\frac{1}{2}\|Xw - y\|^2 + \lambda_1\|w\|_1 + \frac{1}{2}\lambda_2\|w\|^2$ via cyclic coordinate updates with soft-thresholding. Used by `Lasso` and `ElasticNet`.

```csharp
var cd = new CoordinateDescent(maxIterations: 1000, tolerance: 1e-7);
double[] weights = cd.Solve(X, y, l1: 0.1, l2: 0.01, skipBiasRegularisation: true);
```

**Minimizer** — convenience runner

Ties an `IOptimizer`, `IConvergenceCriterion`, and `IObjectiveFunction` together:

```csharp
var minimizer = new Minimizer(
    new Adam(learningRate: 0.01),
    new MaxIterationsOrTolerance(maxIterations: 5000, tolerance: 1e-8));

VectorN x0 = new VectorN(new double[] { 5, -3 });
VectorN xMin = minimizer.Minimize(myObjective, x0);

Console.WriteLine($"Converged in {minimizer.IterationsUsed} iterations, loss = {minimizer.FinalLoss:E3}");
```

### Strategies

**EarlyStopping** — patience-based training halt

```csharp
var es = new EarlyStopping(patience: 10, minDelta: 1e-4);

for (int epoch = 0; epoch < maxEpochs; epoch++)
{
    double valLoss = Evaluate(model, validationData);
    if (es.Check(valLoss))
    {
        Console.WriteLine($"Early stop at epoch {epoch}");
        break;
    }
}
```

**LearningRateSchedule** — lr adjustment over training

```csharp
var schedule = new LearningRateSchedule(
    initialLr: 0.01,
    type: LearningRateSchedule.ScheduleType.CosineAnnealing,
    decaySteps: 1000);

for (int step = 0; step < 1000; step++)
{
    double lr = schedule.GetLearningRate(step);
    // use lr in optimizer
}
```

| Schedule | Formula |
|----------|---------|
| `Constant` | $\eta_0$ |
| `StepDecay` | $\eta_0 \cdot \gamma^{\lfloor t / s \rfloor}$ |
| `ExponentialDecay` | $\eta_0 \cdot \gamma^{t/s}$ |
| `InverseTimeDecay` | $\eta_0 / (1 + \gamma t)$ |
| `CosineAnnealing` | $\eta_{\min} + \frac{1}{2}(\eta_0 - \eta_{\min})(1 + \cos(\pi t / T))$ |

### Multi-Objective: Pareto Front

`ParetoSolution` holds decision variables and objective values. `ParetoFront` provides non-dominated sorting and crowding distance computation.

```csharp
// Create solutions (all objectives are minimised)
var solutions = new List<ParetoSolution>
{
    new ParetoSolution(new[] { 1.0 }, new[] { 0.2, 0.8 }),
    new ParetoSolution(new[] { 2.0 }, new[] { 0.5, 0.5 }),
    new ParetoSolution(new[] { 3.0 }, new[] { 0.9, 0.1 }),
    new ParetoSolution(new[] { 4.0 }, new[] { 0.6, 0.7 }),  // dominated
};

// Non-dominated sorting: partition into ranked fronts
List<List<ParetoSolution>> fronts = ParetoFront.NonDominatedSort(solutions);
// fronts[0] = Pareto-optimal, fronts[1] = next rank, ...

// Or just get the Pareto-optimal front directly
List<ParetoSolution> optimal = ParetoFront.GetFront(solutions);

// Crowding distance (higher = more isolated = more diverse)
ParetoFront.ComputeCrowdingDistance(optimal);
double cd = optimal[0].CrowdingDistance;
```

**Dominance check:**

```csharp
bool dominates = solutionA.Dominates(solutionB);
// true if A ≤ B on all objectives and A < B on at least one
```

### Multi-Objective: NSGA-II

NSGA-II (Deb et al. 2002) is a population-based evolutionary algorithm that converges towards the Pareto front while maintaining solution diversity via crowding distance. Uses SBX crossover and polynomial mutation.

```csharp
// Bi-objective: minimise f₁(x) = x² and f₂(x) = (x-2)²
var nsga2 = new NSGA2(
    evaluate: x => new[] { x[0] * x[0], (x[0] - 2) * (x[0] - 2) },
    numVariables: 1,
    lowerBounds: new[] { -5.0 },
    upperBounds: new[] { 5.0 },
    populationSize: 100,
    generations: 200,
    crossoverRate: 0.9,
    mutationRate: 0.1,
    seed: 42);

NSGA2Result result = nsga2.Run();

Console.WriteLine($"Pareto front: {result.FrontSize} solutions");
foreach (var sol in result.ParetoFront)
    Console.WriteLine($"  x = {sol.Variables[0]:F3}, f₁ = {sol.Objectives[0]:F3}, f₂ = {sol.Objectives[1]:F3}");
```

**Multi-variable, multi-objective (ZDT1 benchmark):**

```csharp
int n = 30;
var nsga2 = new NSGA2(
    evaluate: x =>
    {
        double f1 = x[0];
        double g = 1 + 9 * Enumerable.Range(1, n - 1).Sum(i => x[i]) / (n - 1);
        double f2 = g * (1 - Math.Sqrt(f1 / g));
        return new[] { f1, f2 };
    },
    numVariables: n,
    lowerBounds: Enumerable.Repeat(0.0, n).ToArray(),
    upperBounds: Enumerable.Repeat(1.0, n).ToArray(),
    populationSize: 200,
    generations: 500,
    seed: 123);

NSGA2Result result = nsga2.Run();

// result.ParetoFront contains the approximated Pareto front
// result.FinalPopulation contains the full final population
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `populationSize` | 100 | Number of individuals (rounded up to even) |
| `generations` | 200 | Evolution cycles |
| `crossoverRate` | 0.9 | SBX crossover probability per pair |
| `mutationRate` | 0.1 | Polynomial mutation probability per variable |
| `mutationScale` | 0.1 | Mutation range relative to variable bounds |
| `seed` | 42 | Random seed for reproducibility |

### ML Integration

ML models consume `Numerics.Optimization` — they never reimplement optimisation logic:

| ML Model | Optimiser Used |
|----------|---------------|
| `Logistic` | `GradientDescent` + `MaxIterationsOrTolerance` |
| `Lasso` | `CoordinateDescent` |
| `ElasticNet` | `CoordinateDescent` |
| `MLPRegressor` | `EarlyStopping` |
| `MLPClassifier` | `EarlyStopping` |
| `LinearSVC` / `LinearSVR` | `GradientDescent` |
| `KernelSVC` / `KernelSVR` | `GradientDescent` |
| `NeuralNetwork` | `IOptimizer` (any: GD, Adam, etc.) |

---

## 🔩 Finite Element

Finite-element primitives for 1D beams/bars and 2D plane elasticity.

### 1D Elements

### Elements

| Element | DOFs/Node | Total DOFs | Stiffness |
|---------|-----------|------------|-----------|
| `BarElement` | 1 (u) | 2 | $(EA/L)\begin{bmatrix}1&-1\\-1&1\end{bmatrix}$ |
| `BeamElement` | 2 (w, θ) | 4 | Standard 4×4 Hermite cubic |

### 1D Mesh + Assemble + Solve

```csharp
using CSharpNumerics.Numerics.FiniteElement;

// Cantilever beam: fixed at x=0, tip load P at x=L
double EI = 1e4, L = 1.0, P = 100.0;
int nElem = 10;

var mesh = new Mesh1D(0, L, nElem);
var elements = mesh.CreateElements((_, len) => new BeamElement(EI, len));

var asm = new Assembler1D(mesh, elements);
asm.Assemble();
asm.ApplyNodalLoad(nElem, 0, P); // tip force

var bc = new Dictionary<int, double>
{
    { 0, 0.0 }, // w = 0
    { 1, 0.0 }  // θ = 0
};

VectorN u = asm.Solve(bc);
double tipDeflection = u[nElem * 2]; // ≈ PL³/(3EI)
```

### 1D Bar element example

```csharp
// Fixed-free bar under axial load
double EA = 1000.0, L = 2.0, P = 100.0;
var mesh = new Mesh1D(0, L, 4);
var elements = mesh.CreateElements((_, len) => new BarElement(EA, len));

var asm = new Assembler1D(mesh, elements);
asm.Assemble();
asm.ApplyNodalLoad(4, 0, P);

var u = asm.Solve(new Dictionary<int, double> { { 0, 0.0 } });
// u[i] ≈ P * x[i] / EA
```

### 2D Plane Elasticity

| Type | Nodes | DOFs/Node | Notes |
|------|-------|-----------|-------|
| `TriElement` | 3 | 2 (`ux`, `uy`) | Constant strain triangle (CST) |
| `QuadElement` | 4 | 2 (`ux`, `uy`) | Bilinear Q4 with 2×2 Gauss integration |

Use `Mesh2D` to generate a structured rectangular mesh and `Assembler2D` to build and solve the global system with sparse storage and preconditioned conjugate gradient.

```csharp
using CSharpNumerics.Numerics.FiniteElement;
using CSharpNumerics.Numerics.FiniteElement.Enums;

double width = 2.0, height = 0.5;
int nx = 20, ny = 6;

var mesh = new Mesh2D(width, height, nx, ny, ElementType.Tri);
var asm = new Assembler2D(mesh, new TriElement(), E: 200e9, nu: 0.3, thickness: 0.01, planeStress: true);
asm.Assemble();

var bc = new Dictionary<int, double>();
for (int iy = 0; iy <= ny; iy++)
{
    int node = mesh.GetNodeIndex(0, iy);
    bc[node * 2] = 0.0;
    bc[node * 2 + 1] = 0.0;
}

int tipNode = mesh.GetNodeIndex(nx, ny / 2);
asm.ApplyNodalLoad(tipNode, direction: 1, value: -1000.0);

VectorN u = asm.Solve(bc);
```

### Sparse Matrix Backend

`SparseMatrix` stores global systems in CSR format and provides:

| Feature | API |
|---------|-----|
| Triplet assembly | `SparseMatrix.FromTriplets(...)` |
| Sparse matrix-vector product | `Multiply(VectorN x)` |
| Diagonal extraction | `Diagonal()` |
| Dirichlet elimination | `ApplyDirichlet(...)` |
| Iterative solve | `SolvePCG(...)` |



