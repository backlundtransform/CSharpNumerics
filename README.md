# 🧮 CSharpNumerics

A comprehensive numerical library for **scientific computing**, **mathematical analysis**, and **iterative processes** in C#.
[NuGet Package](https://www.nuget.org/packages/CSharpNumerics/)

---

## ✨ Features

* 🔢 Numerical extensions (Factorial, derivatives, integrals, root finding, etc.)
* 📈 Vectors, matrices, and complex numbers
* 🌊 Vector fields (Gradient, Divergence, Curl, Laplacian)
* 🧠 Complex and real function analysis
* 🔬 Fourier, Laplace, and Monte Carlo transforms
* 📉 Differential equation solvers (Runge–Kutta, Trapezoidal, etc.)
* 📊 Statistics and regression tools
* ✨ Interpolation methods
* 🤖 Machine Learning
* 🔗 Full integration with LINQ and extension methods

---

## 📘 Numeric Extensions

### Factorial

```csharp
int result = 5.Factorial(); // 120
```

### Root Finding (Newton–Raphson)

```csharp
Func<double, double> func = x => Math.Pow(x, 2) - 4;
double root = func.NewtonRaphson(); // 2
```

---

### Derivative

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

### Monte Carlo Integration

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

---

## 🧮 Matrix

```csharp
var A = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });
var transpose = A.Transpose();
var det = A.Determinant();
var inv = A.Inverse();
```

Arithmetic:

```csharp
var B = new Matrix(new double[,] { { 2, 5, 1 }, { 4, 3, 7 } });
var sum = A + B;
var product = A * B;
```

With vector:

```csharp
var x = new Vector(2, 1, 3);
var y = A * x;
```

---

## 🌐 Vector Field

### Gradient

```csharp
Func<Vector, double> f = p => Math.Pow(p.x, 2) * Math.Pow(p.y, 3);
var grad = f.Gradient((1, -2, 0));
```

### Divergence

```csharp
var field = new VectorField(p => Math.Sin(p.x * p.y),
                            p => Math.Cos(p.x * p.y),
                            p => Math.Exp(p.z));

double div = field.Divergence((1, 2, 2));
```

### Curl

```csharp
var field = new VectorField(p => p.y, p => -p.x, p => 0);
var curl = field.Curl((1, 4, 2));
```

---

## ⚙️ Transform

### FFT

```csharp
Func<double, double> f = t => Math.Exp(-t * t / 0.02);
var freq = f.FastFouriertransform(-0.5, 0.5, 100)
             .ToFrequencyResolution(100);
```

### Laplace Transform

```csharp
double result = f.LaplaceTransform(2.0);
```

---

## 📐 Differential Equations

### Runge–Kutta (RK4)

```csharp
Func<(double t, double y), double> f = v => Math.Tan(v.y) + 1;
var result = f.RungeKutta(1, 1.1, 0.025, 1);
```

### Linear Systems

```csharp
var result = A.LinearSystemSolver(b);
var eigenValues = A.EigenValues();
```

---

## 📊 Statistics

```csharp
var noise = new Random().GenerateNoise(4);
double median = ts.Median(p => p.Value);
double std = ts.StandardDeviation(p => p.Value);
```
Coefficient of determination:
```csharp
var data = new[] { (1.0, 5.0), (2.0, 1.0), (3.0, 4.0), (4.0, 6.0) };
double r2 = data.CoefficientOfDetermination(p => (p.Item1, p.Item2));
```

Regression:

```csharp
var (slope, intercept, corr) = serie.LinearRegression(p => (p.Index, p.Value));
var expFunc = serie.ExponentialRegression(p => (p.Index, p.Value));
```


K-nearest neighbors:

```csharp
var data = new List<(double x, double y, int c)> { (7,7,0), (7,4,0), (3,4,1), (1,4,1) };
int classification = data.KnearestNeighbors(p => (p.x, p.y, p.c), (3,7), 3);
```

---
##  ✨ Interpolation

CSharpNumerics provides a unified interpolation API supporting linear and logarithmic scales:

* **Linear**
* **Log–Log** (log x, log y)
* **Lin–Log** (lin x, log y)
* **Log–Lin** (log x, lin y)

All methods are routed through one central function.


```csharp
public enum InterpolationType
{
    Linear,
    Logarithmic,   // log–log
    LogLin,
    LinLog
}
```

```csharp
double Interpolate<T>(
    this IEnumerable<T> source,
    Func<T, (double x, double y)> selector,
    double index,
    InterpolationType type);
```

Example:

```csharp
var data = new List<Serie>
{
    new Serie { Index = 1, Value = 10 },
    new Serie { Index = 10, Value = 100 }
};

double y = data.Interpolate(
    p => (p.Index, p.Value),
    3.5,
    InterpolationType.Linear
);
```

---

## 🤖 Machine Learning

CSharpNumerics includes a **lightweight, fully numerical machine learning framework** designed for **research, experimentation, and educational use**.
The focus is on **transparency**, **mathematical clarity**, and **pipeline-based model evaluation** — not black-box automation.

All models are implemented directly on top of the library’s **Matrix** and **Vector** primitives.


Models can be combined with:

* **Scalers** (e.g. StandardScaler)
* **Feature selectors** (e.g. SelectKBest)
* **Cross-validation strategies**
* **Hyperparameter search grids**

```csharp
    var pipelineGrid = new PipelineGrid()
    .AddModel<RandomForest>(g => g
        .Add("NumTrees", 50, 100, 200)
        .Add("MaxDepth", 5, 8, 10))
    .AddModel<Logistic>(g => g
        .Add("LearningRate", 0.05, 0.1)
        .Add("MaxIterations", 1000, 2000)
        .AddScaler<StandardScaler>(s => {})
        .AddSelector<SelectKBest>(s => s
           .Add("K", 1, 2)))
    .AddModel<DecisionTree>(g => g
        .Add("MaxDepth", 3, 5, 8))
    .AddModel<KNearestNeighbors>(g => g
        .Add("K", 3, 5, 7));
```

Rolling (time-aware) cross validation is supported for both **classification** and **regression** tasks.

```csharp
var cv = new RollingCrossValidator(pipelineGrid, folds: 5);
var result = cv.Run(X, y);

var bestModel = result.BestPipeline;
var score = result.BestScore;
```

####  📊 Classification Models

Supported classifiers include:

* Logistic Regression
* Decision Tree
* Random Forest
* K-Nearest Neighbors
* Naive Bayes
* Linear SVC


####  📈 Regression Models

Supported regression models:

* Linear Regression
* Ridge Regression (L2)
* Lasso Regression (L1)
* Elastic Net (L1 + L2)
* Linear SVR 

---

## 📎 Tips

* All methods are available as **extension methods** — just `using Numerics.Extensions`.
* You can **export data** with `.Save(path)` for CSV visualization.
* Works with LINQ pipelines for composable scientific workflows.





## 🧠 Example: Full Workflow

```csharp
Func<double, double> func = x => Math.Sin(x);
var integral = func.Integrate(0, Math.PI);
var derivative = func.Derivate(Math.PI / 4);
var fft = func.FastFouriertransform(-1, 1, 100);
```

---

## 🧾 License

MIT License © 2025 — [CSharpNumerics](https://www.nuget.org/packages/CSharpNumerics)


