using System;

namespace CSharpNumerics.Numerics.Objects
{
    public struct ComplexVector
    {
        public ComplexNumber x;

        public ComplexNumber y;

        public ComplexNumber z;


        public ComplexVector(ComplexNumber a, ComplexNumber b, ComplexNumber c)
        {
            x = a;
            y = b;
            z = c;
        }

        public ComplexVector(double a, double b, double c)
        {
            x = new ComplexNumber(a, 0);
            y = new ComplexNumber(b, 0);
            z = new ComplexNumber(c, 0);
        }

        public double GetMagnitude() => Math.Sqrt(
            x.GetMagnitude() * x.GetMagnitude() +
            y.GetMagnitude() * y.GetMagnitude() +
            z.GetMagnitude() * z.GetMagnitude());

        public ComplexVector GetConjugate() => new ComplexVector(x.GetConjugate(), y.GetConjugate(), z.GetConjugate());

        public ComplexVector GetUnitVector()
        {
            var mag = new ComplexNumber(GetMagnitude(), 0);
            return new ComplexVector(x / mag, y / mag, z / mag);
        }

        public static ComplexVector operator +(ComplexVector a, ComplexVector b) =>
            new ComplexVector(a.x + b.x, a.y + b.y, a.z + b.z);

        public static ComplexVector operator -(ComplexVector a, ComplexVector b) =>
            new ComplexVector(a.x - b.x, a.y - b.y, a.z - b.z);

        public static ComplexVector operator *(ComplexNumber a, ComplexVector b) =>
            new ComplexVector(a * b.x, a * b.y, a * b.z);

        public static ComplexVector operator *(double a, ComplexVector b) =>
            new ComplexNumber(a, 0) * b;

        public static ComplexVector operator /(ComplexVector a, ComplexNumber b) =>
            new ComplexVector(a.x / b, a.y / b, a.z / b);

        public static ComplexVector operator /(ComplexVector a, double b) =>
            a / new ComplexNumber(b, 0);

        public ComplexNumber Dot(ComplexVector b) => x * b.x + y * b.y + z * b.z;

        public ComplexNumber HermitianDot(ComplexVector b) =>
            x.GetConjugate() * b.x + y.GetConjugate() * b.y + z.GetConjugate() * b.z;

        public ComplexVector Cross(ComplexVector b) => new ComplexVector(
            y * b.z - b.y * z,
            new ComplexNumber(0, 0) - (x * b.z - b.x * z),
            x * b.y - b.x * y);

        public static implicit operator ComplexVector(Vector v) =>
            new ComplexVector(
                new ComplexNumber(v.x, 0),
                new ComplexNumber(v.y, 0),
                new ComplexNumber(v.z, 0));

        public override string ToString() => $"{x}*î̤+{y}*ĵ+{z}*k̂";
    }
}
