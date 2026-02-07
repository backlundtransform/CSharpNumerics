
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

Trapezoidal rule:

```csharp
Func<double, double> f = x => Math.Sin(x);
double integral = f.Integrate(0, Math.PI);
```

Integrate a series or timeseries:

```csharp
List<TimeSerie> ts = ...;
double total = ts.Integrate();
```

**Monte Carlo Integration**

```csharp
Func<(double x, double y), double> func = p => p.x * p.y;
double result = func.Integrate((0, 1), (0, 1));
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

**FFT**

```csharp
Func<double, double> f = t => Math.Exp(-t * t / 0.02);
var freq = f.FastFouriertransform(-0.5, 0.5, 100)
             .ToFrequencyResolution(100);
```

**Laplace Transform**

```csharp
double result = f.LaplaceTransform(2.0);
```

---

## 📐 Differential Equations

**Runge–Kutta (RK4)**

```csharp
Func<(double t, double y), double> f = v => Math.Tan(v.y) + 1;
var result = f.RungeKutta(1, 1.1, 0.025, 1);
```

## 📏 Linear Systems

```csharp
var result = A.LinearSystemSolver(b);
var eigenValues = A.EigenValues();
```
```csharp
var matrix = new Matrix(new double[,] { { 1, -2, 3 }, { -1, 1, -2 }, { 2, -1, -1 } });

var vector = new VectorN(new double[] { 7, -5, 4 });

var result = matrix.LinearSystemSolver(vector);

```



