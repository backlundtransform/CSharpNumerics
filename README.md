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
      `var i = new ComplexNumber(3, 2);
        i.Pow(2);`
        
output:  5+12*i

Calculate Imaginary exponents

  `var i = new ComplexNumber(0, Math.PI);
    i.Exponential()`

output:  -1

## The vector object

To work with Vectors use this struct:

  `Vector(double x, double y, double z)`

or from two points

  `Vector((double,double, double) p1, (double, double, double) p2)`

Following extension method could be used:

1. Skalar product

  `Dot(this Vector a, Vector b)`

2. Vector product

 `Cross(this Vector a, Vector b)`

3. Projection between two vectors

  `Projection(this Vector a, Vector b)`

4. Reflection between two vectors

   `Reflection(this Vector a, Vector b)`

 

   
	





