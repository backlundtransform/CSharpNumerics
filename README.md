# CSharpNumerics
Simple numeric package


## Numeric Extensions

To derivate a function use:

`Derivate(this Func<double, double> func, double variablevalue)`

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
