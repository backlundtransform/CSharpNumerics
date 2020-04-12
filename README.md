  # CSharpNumerics
Numeric package


## Numeric Extensions

Get the Factorial of int 

`Factorial(this int number)`

E.g 

`5.Factorial()`

Outputs 120

### Derivative

To derivate a function use:

`Derivate(this Func<double, double> func, double variablevalue, int order=1)`

Calculate higher order derivative by setting the order parameter 

E.g with Chain rule

`double funcG(double x) => 4 * x - 3`
 
`Func<double, double> funcF=(double x) => Math.Pow(x, 2)`
 
`var result = funcF.Derivate(funcG,1)`

By using Numerics.Enums.DerivateOperator the Chain rule, Product rule, or Quotient rule can be used

`var result = funcF.Derivate(funcG,Numerics.Enums.DerivateOperator.Product)`

`var result = funcF.Derivate(funcG,Numerics.Enums.DerivateOperator.Quotient)`

If several variables use:

`Derivate(this Func<double[], double> func, double[] variables, int index, int order=1)`

Or use the vector (x,y,z)  

`Derivate(this Func<Vector, double> func, Vector variables, Cartesian cartesian, int order=1)`

### Integrals

To integrate a function with Trapezoidal rule use:

`Integrate(this Func<double, double> func, double lowerLimit, double upperLimit)`

To integrate a timeserie 

`Integrate(this List<Numerics.Models.TimeSerie> data)`

TimeSerie is a model with properties TimeStamp as DateTime and Value as double 

To integrate a serie  

`Integrate(this List<Numerics.Models.Serie> data)`

Serie model is a model with properties Index as souble and Value as double

### Monte Carlo Integration

To solve double integrals with Monte Carlo method use:

`Integrate(this Func<(double x, double y), double> func, (double lowerLimit, double upperLimit) xlimit, (double lowerLimit, double upperLimit) ylimit)`

or triple integral

`Integrate(this Func<Vector, double> func, Vector lowerLimit, Vector upperLimit)`



## The complex object

To work with Complex numbers use this struct: 

`ComplexNumber(double re, double im)`

E.g Arithmetics

  `var a = new ComplexNumber(3, 2);`
  
   `var b = new ComplexNumber(5, 3);`
   
   `var sum= a + b;`
   
   `var difference =a-b;`
   
   `var product =a*b;`
   
   `var quotient= a / b;`
   
Power of complex number
   
   `var i = new ComplexNumber(3, 2); i.Pow(2);`
        
output:  5+12*i

Calculate Imaginary exponents

  `var i = new ComplexNumber(0, Math.PI);
    i.Exponential()`

output:  -1

## The vector object

To work with vectors use this struct:

  `Vector(double x, double y, double z)`

or from two points

  `Vector((double,double, double) p1, (double, double, double) p2)`

Following methods could be used:

1. Scalar product

  `Dot(Vector b)`

2. Vector product

 `Cross(Vector b)`

3. Projection between two vectors

  `Projection(Vector b)`

4. Reflection between two vectors

  `Reflection(Vector b)`

 Using Sperical Coordinates

   `var v=Vector.FromSphericalCoordinates(radius, inclination, azimuth)`

or covert to sperical from cartesian

   `v.ToSphericalCoordinates()`

Metods to get radius, inclination, azimuth

   `GetMagnitude(), GetInclination(),GetAzimuth()`

E.g 

  `var a = new Vector(5, 3, 0);`

   `var b = new Vector(2, 6, 0);`

   `var skalar = a.Dot(b);`

   `var vector = a.Cross(b);`


E.g Arithmetics

  `var a = new Vector(2, 2, 0);`
  
   `var b = new Vector(2, 2, 0);`
   
   `var sum= a + b;`
   
   `var difference =a-b;`

   `var product =3*b;`

## The matrix object	

To work with matrix use this struct:

  `Matrix(double[,] values)`

E.g 

  `var matrix = new Matrix(new double[,] { { 1, 3, 7 }, { 5, 2, 9 } });`

Get Transpose

`var transposematrix = matrix.Transpose();`

Get Inverse

`var inv = matrix.Inverse();`


Get Pascal matrix

`var pascalMatrix = new Matrix(new double[6, 6]).Pascal()`

Get Adjugate

`var adj = matrix.Adjugate();`

Get Determinant

` var det = matrix.Determinant()`

output values 

`matrix.values`

or Identity matrix

 `matrix.identity`

E.g Arithmetics

  `var a = new Matrix(new double[,] { { 5, 7, 2 }, { -2, 9, 4 } });`
  
  `var b = new Matrix(new double[,] {{ 1, 3, 7 }, { 5, 2, 9} });`
   
   `var sum= a + b;`
   
   `var difference =a-b;`
   
   `var product =a*b;`
   
Or

  `var b = 3`

  `var product =b*a;`

  `var quotient= a / b;`


Or  

   `var a = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 }});`

   `var b = new Vector(2, 1, 3);`

   `var product =a*b;`

## The vectorfield object

### Gradient

Calculate for one point

  `Func<Vector, double> func = (Vector p) => Math.Pow(p.x, 2) * Math.Pow(p.y, 3);`

  `var v=func.Gradient((1, -2, 0))`

Calculate for range

 `var grad = func.Gradient(-4, -4, -4, 1, 8);`

Where the parameters are minimum value of x,y,z the step size, and the length. The range in this example is -4<=x<=4,-4<=y<=4,-4<=z<=4. 
This method will return new Dictionary<Vector, Vector> where the key is the point and value is the calculated function value for that point.
Save the data to csv 

  `grad.Save(@"${path}\${file}.csv");`

### Divergence

Calculate for one point

 `double fx(Vector p) => Math.Sin(p.x * p.y);`

 `double fy(Vector p) => Math.Cos(p.x * p.y);`

 `double fz(Vector p) => Math.Pow(Math.E, p.z);`

 `var field = new VectorField(fx, fy, fz);`

 `var div = field.Divergence((1, 2, 2))`

### Curl

Calculate for one point
           
 `double fx(Vector p) => 4*p.z;`

 `double fy(Vector p) => p.y *Math.Pow(p.x,3);`

 `double fz(Vector p) => p.z * Math.Pow(p.y,2);`

 `var field = new VectorField(fx, fy, fz);`

 `var v = field.Curl((1, 4, 2));`

Calculate for range. It is done in same way as for gradient E.g save both the vector field and curl to file:

  `double fx(Vector p) => p.y;`

  `double fy(Vector p) => -p.x;`

   `double fz(Vector p) => 0;`

   `var w = new VectorField(fx, fy, fz);`
    
   `var data= w.EvaluateRange(-4,-4,4,1,8);`
            
   `var curl = w.Curl(-4, -4, 4, 1, 8);`

  `data.Save(@"${path}\${file}.csv");`

  `curl.Save(@"${path}\${file}.csv");`

### Laplacian

Calculate for one point

  `Func<Vector, double> func = (Vector p) => Math.Pow(p.x, 2) * Math.Pow(p.y, 3);`

  `var v=func.Laplacian((1, -2, 0))`

## The complex function object

To work with Complex functions use this struct: 

 `ComplexFunction(Func<(double x, double y), double> re, Func<(double x, double y), double> im)`

E.g:

  `double fx((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Cos(p.y);`

  `double fy((double x, double y) p) => Math.Pow(Math.E, p.x) * Math.Sin(p.y);`

   `var fz = new ComplexFunction(fx, fy);`

or 

   `ComplexNumber fz(ComplexNumber z) => new ComplexNumber(Math.Pow(Math.E, z.realPart) * Math.Cos(z.imaginaryPart), Math.Pow(Math.E, z.realPart) * Math.Sin(z.imaginaryPart));` 

 To derivate a complex function use:
            
   `Derivate(this ComplexFunction func, ComplexNumber variables, int order = 1)`

### Cauchy–Riemann equations

To test if analytic fuction in a point using Cauchy–Riemann equations:
	
   `fz.IsAnalytical((x0,y0))`

### Jacobian

To get the Jacobian as a Matrix

   `fz.Jacobian((x0,y0))`


## Transform

E.g use a lowpass filter to remove noise from a signal

`var result =input.LowPassFilter(output).ToList()`

### Fast Fourier transform

To use a fast fourier transform use extentionsmethods
from a list of complexnumber 

`FastFouriertransform(this List<ComplexNumber> numbers)`

to calculate the inverse fast fourier transform

`List<ComplexNumber> InverseFastFouriertransform(this List<ComplexNumber> numbers)`

E.g Convert a Gaussian pulse from the time domain to the frequency domain and save result.

`Func<double, double> func = (double t) => 1 / (4 * Math.Sqrt(2 * Math.PI * 0.01)) * (Math.Exp(-t * t / (2 * 0.01)));`

`var timeseries = func.GetSeries(-0.5, 0.5, 100);`

`timeseries.Save(@"\timeserie.csv");`

GetSeries takes the interval and how many values to return

`var frequency = func.FastFouriertransform(-0.5, 0.5, 100).ToFrequencyResolution(100);`

`frequency.Save(@"\frequency.csv");`


ToFrequencyResolution takes the sample rate and will return the frequency as index and the magnitude of the complex number as value

### Discrete Fourier transform

Use in the same way as a fast fourier transform 

`DiscreteFourierTransform(this List<ComplexNumber> numbers)`

### LaplaceTransform

Calculate the laplace transform for s value

`LaplaceTransform(this Func<double, double> func, double s)`

or it's invers from t value

`InverseLaplaceTransform(this Func<double, double> func, double t)`

## Differential Equations

### Linear equation systems 

Extension methods to solve linear equation system  

`LinearSystemSolver(this Matrix matrix, Vector vector)`

or gauss elimination

`GaussElimination(this Matrix matrix, Vector vector)`

for N values

`GaussElimination(this Matrix matrix, List<double> vector)`

Find eigen values of matrix 

`var result = matrix.EigenValues()`

If knowing a eigenvalue of a matrix (in this example 1)

`var result =matrix.EigenVector(1)`

### Runge–Kutta

The Runge–Kutta (R4) method uses this extension method 

`RungeKutta(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)`

E.g yprim =tan(y)+1 with the initial-value problem  y0=1 and 1<= t <= 1.1 and step size 0.025

`Func<(double y, double t), double> func = ((double t, double y) v) => Math.Tan(v.y) +1`
         
`var result = func.RungeKutta(1,1.1,0.025,1)`

It is also possible using Explicit Runge–Kutta methods by defining the Runge–Kutta matrix, weights and nodes

`var result = func.RungeKutta(1,1.1,0.025,1,new Matrix(new double[,] { { 0, 0, 0 }, { 0.5, 0, 0 }, { 0, 0.5, 0 }, { 0, 0, 1 } }), new double[] { 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 }, new double[] { 0.0, 0.5, 0.5, 1 })`

Or solve ode by using the Trapezoidal rule

`var result = func.TrapezoidalRule(1, 1.1, 0.00025, 1);`


## Statistics

Generate zero-mean white noise with a variance of 4 using Random

`var rnd = new Random()`

`rnd.GenerateNoise(4)`


To calculate the median use linq in the same way as calculating avarerage, sum, max or min

`timeseries.Median(p => p.Value)`

To calculate the standard deviation

`timeseries.StandardDeviation(p => p.Value)`

To calculate the variance

`timeseries.Variance(p => p.Value)`

To calculate the covariance if model has X,Y properties

`series.Covariance(p => (p.X,p.Y))`

There is also a Statistics class containing static methods

E.g  get normal distribution curve 

 `Numerics.Methods.Statistics.NormalDistribution(variance, mean)`

Extension method for calculating cumulative sum of list

`CumulativeSum<T>(this IEnumerable<T> enumerable, Func<T, double> func)`

### Interpolation

Linear interpolation of a timeserie

`LinearInterpolationTimeSerie(this IEnumerable<TimeSerie> ts, DateTime timeStamp)`

Or serie

`LinearInterpolation<T>(this IEnumerable<T> ts,  Func<T, (double x, double y)> func, double index)`

E.g

`var value = serie.LinearInterpolation(p=>(p.Index, p.Value),1)` 


### Regression

Linear regression that will return intercept correlation and slope

`LinearRegression<T>(this IEnumerable<T> enumerable, Func<T, (double x, double y)> func)`

E.g

`var serie = new List<Serie>() 
{ new Serie() { Index = 3.0, Value = 0.62},
new Serie() { Index = 3.4, Value = 0.93 },
new Serie() { Index = 3.8, Value = 1.08 }};`

`var (slope, intercept, correlation) = serie.LinearRegression(p=>(p.Index, p.Value))`

Exponetial regression  that will return a exponetial function

`var func = serie.ExponentialRegression(p => (p.Index, p.Value));`






