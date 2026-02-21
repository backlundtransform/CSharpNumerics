using System;

namespace CSharpNumerics.Numerics.Objects
{
   
    public class ComplexNumber
    {
        public double realPart;

        public double imaginaryPart;

        public ComplexNumber(double re, double im)
        {
            realPart = re;
            imaginaryPart = im;
        }

        public override string ToString() => $"{realPart}+{imaginaryPart}*i";

        public double GetArgument() => Math.Atan2(imaginaryPart,realPart);

        public double GetMagnitude() => Math.Sqrt(Math.Pow(realPart, 2) + Math.Pow(imaginaryPart, 2));

        public ComplexNumber GetConjugate() => new ComplexNumber(realPart, -imaginaryPart);

     
        public static ComplexNumber FromPolarCoordinates(double r, double theta) => new ComplexNumber(r * Math.Cos(theta), r * Math.Sin(theta));

        public static ComplexNumber operator +(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.realPart + b.realPart, a.imaginaryPart + b.imaginaryPart);
        }
        public static ComplexNumber operator -(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.realPart - b.realPart, a.imaginaryPart - b.imaginaryPart);
        }
        public static ComplexNumber operator *(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber((a.realPart * b.realPart - a.imaginaryPart * b.imaginaryPart),
                (a.realPart * b.imaginaryPart + a.imaginaryPart * b.realPart));
        }
        public static ComplexNumber operator /(ComplexNumber a, ComplexNumber b)
        {
            var denominator = Math.Pow(b.realPart, 2) + Math.Pow(b.imaginaryPart, 2);
         
            return new ComplexNumber((a.realPart * b.realPart + a.imaginaryPart * b.imaginaryPart) / denominator,
                (a.imaginaryPart * b.realPart - a.realPart * b.imaginaryPart) / denominator);
        }

        public static ComplexNumber operator /(ComplexNumber a,double b)=>new ComplexNumber(a.realPart / b, a.imaginaryPart);

        public static ComplexNumber operator *(double a, ComplexNumber b) => new ComplexNumber(b.realPart * a, b.imaginaryPart);

        public static ComplexNumber operator *(ComplexNumber a, double b) => new ComplexNumber(a.realPart* b, a.imaginaryPart);

        public static ComplexNumber operator +(ComplexNumber a, double b) => new ComplexNumber(a.realPart + b, a.imaginaryPart);

        public static ComplexNumber operator -(double a, ComplexNumber b) => new ComplexNumber(a-b.realPart, b.imaginaryPart);

        public static ComplexNumber operator -(ComplexNumber a, double b) => new ComplexNumber(a.realPart - b, a.imaginaryPart);

        public static ComplexNumber operator +(double a, ComplexNumber b) => new ComplexNumber(a+b.realPart, b.imaginaryPart);


    }
}
