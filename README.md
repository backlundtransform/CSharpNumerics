# CSharpNumerics
Simple numeric package


## Numeric Extensions

To derivate a function use:

`Derivate(this Func<double, double> func, double variablevalue)`

If several variables use;

`PartialDerivate(this Func<double[], double> func, double[] variables, int index)`

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
   `var a = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });`

    `var b = new Vector(2, 1, 3);`

    `var product =a*b;`


