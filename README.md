# CSharpNumerics
Simple numeric package


## Numeric Extensions

To derivate a function use:

`Derivate(this Func<double, double> func, double variablevalue, int order=1)`

Calculate higher order derivative by setting the order parameter 

If several variables use:

`Derivate(this Func<double[], double> func, double[] variables, int index, int order=1)`

Or use the vector (x,y,z)  

`Derivate(this Func<Vector, double> func, Vector variables, Cartesian cartesian, int order=1)`

To integrate a function use:

`Integrate(this Func<double, double> func, double lowerlimit, double upperlimit)`


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

 To derivate a complex function use:
            
   `Derivate(this ComplexFunction func, ComplexNumber variables, int order = 1)`

### Cauchy–Riemann equations

To test if analytic fuction in a point using Cauchy–Riemann equations:
	
   `fz.IsAnalytical((x0,y0))`

### Jacobian

To get the Jacobian as a Matrix

   `fz.Jacobian((x0,y0))`
