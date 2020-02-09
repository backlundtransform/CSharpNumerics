using System;

namespace Numerics.Objects
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

        public ComplexNumber GetConjugate() => new ComplexNumber(realPart, imaginaryPart);

        public void Exponential()
        {
            realPart = Math.Pow(Math.E, realPart) * Math.Cos(imaginaryPart);
            imaginaryPart = Math.Sin(imaginaryPart);
        }

        public void Pow(int power)
        {
            var mag = GetMagnitude();
            var arg = GetArgument();
            realPart = Math.Pow(mag, power) * Math.Cos(power * arg);
            imaginaryPart = Math.Pow(mag, power) * Math.Sin(power * arg);
        }

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
                
   
    }
}
