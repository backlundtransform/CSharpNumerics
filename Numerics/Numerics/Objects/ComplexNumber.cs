using System;

namespace Numerics
{
   
    public class ComplexNumber
    {
        public double RealPart;

        public double ImaginaryPart;

        public ComplexNumber(double re, double im)
        {
            RealPart = re;
            ImaginaryPart = im;
        }

        public override string ToString() => $"{RealPart}+{ImaginaryPart}*i";

        public double GetArgument() => Math.Atan2(ImaginaryPart, RealPart);

        public double GetMagnitude() => Math.Sqrt(Math.Pow(RealPart, 2) + Math.Pow(ImaginaryPart, 2));

        public ComplexNumber GetConjugate() => new ComplexNumber(RealPart, ImaginaryPart);

        public void Exponential()
        {
            var im = ImaginaryPart;
            RealPart = Math.Pow(Math.E, RealPart) * Math.Cos(ImaginaryPart);
            ImaginaryPart = Math.Sin(ImaginaryPart);
        }

        public void Pow(int power)
        {
            var mag = GetMagnitude();
            var arg = GetArgument();
            RealPart = Math.Pow(mag, power) * Math.Cos(power * arg);
            ImaginaryPart = Math.Pow(mag, power) * Math.Sin(power * arg);
        }

        public static ComplexNumber FromPolarCoordinates(double r, double theta) => new ComplexNumber(r * Math.Cos(theta), r * Math.Sin(theta));

        public static ComplexNumber operator +(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.RealPart + b.RealPart, a.ImaginaryPart + b.ImaginaryPart);
        }
        public static ComplexNumber operator -(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.RealPart - b.RealPart, a.ImaginaryPart - b.ImaginaryPart);
        }
        public static ComplexNumber operator *(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber((a.RealPart * b.RealPart - a.ImaginaryPart * b.ImaginaryPart),
                (a.RealPart * b.ImaginaryPart + a.ImaginaryPart * b.RealPart));
        }
        public static ComplexNumber operator /(ComplexNumber a, ComplexNumber b)
        {
            var denominator = Math.Pow(b.RealPart, 2) + Math.Pow(b.ImaginaryPart, 2);
            return new ComplexNumber((a.RealPart * b.RealPart + a.ImaginaryPart * b.ImaginaryPart) / denominator,
                (a.ImaginaryPart * b.RealPart - a.RealPart * b.ImaginaryPart) / denominator);
        }
    }
}
