
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

## 📏 Linear Systems

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







