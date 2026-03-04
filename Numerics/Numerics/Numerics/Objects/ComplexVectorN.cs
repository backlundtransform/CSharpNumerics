using System;

namespace CSharpNumerics.Numerics.Objects
{
    public struct ComplexVectorN
    {
        public ComplexNumber[] Values;


        public ComplexVectorN(int length)
        {
            if (length <= 0)
                throw new ArgumentException("Length must be positive", nameof(length));
            Values = new ComplexNumber[length];
            for (int i = 0; i < length; i++)
                Values[i] = new ComplexNumber(0, 0);
        }


        public ComplexVectorN(ComplexNumber[] values)
        {
            if (values == null || values.Length == 0)
                throw new ArgumentException("Values cannot be null or empty", nameof(values));
            Values = (ComplexNumber[])values.Clone();
        }


        public int Length => Values.Length;


        public ComplexNumber this[int i]
        {
            get => Values[i];
            set => Values[i] = value;
        }


        public static ComplexVectorN operator +(ComplexVectorN a, ComplexVectorN b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");
            var result = new ComplexNumber[a.Length];
            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] + b[i];
            return new ComplexVectorN(result);
        }


        public static ComplexVectorN operator -(ComplexVectorN a, ComplexVectorN b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");
            var result = new ComplexNumber[a.Length];
            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] - b[i];
            return new ComplexVectorN(result);
        }

        public static ComplexVectorN operator *(ComplexNumber scalar, ComplexVectorN v)
        {
            var result = new ComplexNumber[v.Length];
            for (int i = 0; i < v.Length; i++)
                result[i] = scalar * v[i];
            return new ComplexVectorN(result);
        }

        public static ComplexVectorN operator *(double scalar, ComplexVectorN v) =>
            new ComplexNumber(scalar, 0) * v;

        public static ComplexVectorN operator /(ComplexVectorN v, ComplexNumber scalar)
        {
            var result = new ComplexNumber[v.Length];
            for (int i = 0; i < v.Length; i++)
                result[i] = v[i] / scalar;
            return new ComplexVectorN(result);
        }

        public static ComplexVectorN operator /(ComplexVectorN v, double scalar) =>
            v / new ComplexNumber(scalar, 0);


        public ComplexNumber Dot(ComplexVectorN other)
        {
            if (Length != other.Length)
                throw new ArgumentException("Vectors must be same length");
            var sum = new ComplexNumber(0, 0);
            for (int i = 0; i < Length; i++)
                sum = sum + Values[i] * other[i];
            return sum;
        }

        public ComplexNumber HermitianDot(ComplexVectorN other)
        {
            if (Length != other.Length)
                throw new ArgumentException("Vectors must be same length");
            var sum = new ComplexNumber(0, 0);
            for (int i = 0; i < Length; i++)
                sum = sum + Values[i].GetConjugate() * other[i];
            return sum;
        }

        public ComplexVectorN Hadamard(ComplexVectorN b)
        {
            if (Length != b.Length)
                throw new ArgumentException("Vectors must be of the same length for Hadamard product.");

            var result = new ComplexNumber[Length];
            for (int i = 0; i < Length; i++)
                result[i] = Values[i] * b.Values[i];
            return new ComplexVectorN(result);
        }

        public ComplexMatrix Outer(ComplexVectorN b)
        {
            var result = new ComplexNumber[Length, b.Length];
            for (int i = 0; i < Length; i++)
                for (int j = 0; j < b.Length; j++)
                    result[i, j] = Values[i] * b.Values[j];
            return new ComplexMatrix(result);
        }

        public double Norm()
        {
            double sum = 0;
            for (int i = 0; i < Length; i++)
                sum += Values[i].GetMagnitude() * Values[i].GetMagnitude();
            return Math.Sqrt(sum);
        }


        public ComplexVectorN Normalize()
        {
            double norm = Norm();
            if (norm == 0) return new ComplexVectorN(Length);
            return this / new ComplexNumber(norm, 0);
        }

        public ComplexVectorN GetConjugate()
        {
            var result = new ComplexNumber[Length];
            for (int i = 0; i < Length; i++)
                result[i] = Values[i].GetConjugate();
            return new ComplexVectorN(result);
        }

        public static implicit operator ComplexVectorN(VectorN v)
        {
            var values = new ComplexNumber[v.Length];
            for (int i = 0; i < v.Length; i++)
                values[i] = new ComplexNumber(v[i], 0);
            return new ComplexVectorN(values);
        }

        public override string ToString()
        {
            var parts = new string[Length];
            for (int i = 0; i < Length; i++)
                parts[i] = Values[i].ToString();
            return "[" + string.Join(", ", parts) + "]";
        }
    }
}
