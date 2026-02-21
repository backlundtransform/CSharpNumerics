using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;


namespace CSharpNumerics.Numerics.Objects
{
    public struct VectorN
    {
        public double[] Values;

 
        public VectorN(int length)
        {
            if (length <= 0)
                throw new ArgumentException("Length must be positive", nameof(length));
            Values = new double[length];
        }


        public VectorN(double[] values)
        {
            if (values == null || values.Length == 0)
                throw new ArgumentException("Values cannot be null or empty", nameof(values));
            Values = (double[])values.Clone();
        }

   
        public int Length => Values.Length;

       
        public double this[int i]
        {
            get => Values[i];
            set => Values[i] = value;
        }

   
        public static VectorN operator +(VectorN a, VectorN b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");
            var result = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] + b[i];
            return new VectorN(result);
        }

      
        public static VectorN operator -(VectorN a, VectorN b)
        {
            if (a.Length != b.Length)
                throw new ArgumentException("Vectors must be same length");
            var result = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
                result[i] = a[i] - b[i];
            return new VectorN(result);
        }

        public static VectorN operator *(double scalar, VectorN v)
        {
            var result = new double[v.Length];
            for (int i = 0; i < v.Length; i++)
                result[i] = scalar * v[i];
            return new VectorN(result);
        }

        public static VectorN operator /(VectorN v, double scalar)
        {
            if (scalar == 0)
                throw new DivideByZeroException();
            var result = new double[v.Length];
            for (int i = 0; i < v.Length; i++)
                result[i] = v[i] / scalar;
            return new VectorN(result);
        }

  
        public double Dot(VectorN other)
        {
            if (Length != other.Length)
                throw new ArgumentException("Vectors must be same length");
            double sum = 0;
            for (int i = 0; i < Length; i++)
                sum += Values[i] * other[i];
            return sum;
        }
        public VectorN Hadamard( VectorN b)
        {
            if (Length != b.Length)
                throw new ArgumentException("Vectors must be of the same length for Hadamard product.");

            var result = new double[Length];
            for (int i = 0; i < Length; i++)
                result[i] = Values[i] * b.Values[i];
            return new VectorN(result);
            
        }

        public Matrix Outer(VectorN b)
        {
            var result = new double[Length, b.Length];
            for (int i = 0; i < Length; i++)
            {
                for (int j = 0; j < b.Length; j++)
                {
                    result[i, j] = Values[i] * b.Values[j];
                }
            }
            return new Matrix { values = result, rowLength = Length, columnLength = b.Length };
        }

        public double Norm()
        {
            double sum = 0;
            for (int i = 0; i < Length; i++)
                sum += Values[i] * Values[i];
            return Math.Sqrt(sum);
        }

       
        public VectorN Normalize()
        {
            double norm = Norm();
            if (norm == 0) return new VectorN(Length);
            return this / norm;
        }

        public VectorN Slice(int[] indices)
        {
            var slice = new double[indices.Length];
            for (int i = 0; i < indices.Length; i++)
                slice[i] = Values[indices[i]];
            return new VectorN(slice);
        }

        public Matrix BuildConfusionMatrix( VectorN yPred)
        {
            double maxTrue =Values[0];
            double maxPred = yPred.Values[0];


            for (int i = 1; i < Length; i++)
            {
                if (Values[i] > maxTrue)
                    maxTrue = Values[i];
                if (yPred.Values[i] > maxPred)
                    maxPred = yPred.Values[i];
            }

            int numClasses = (int)Math.Max(maxTrue, maxPred) + 1;
            var cm = new Matrix(numClasses, numClasses);

            for (int i = 0; i < Length; i++)
            {
                int actual = (int)Values[i];
                int predicted = (int)Math.Round(yPred.Values[i]);

                if (predicted >= 0 && predicted < numClasses)
                    cm.values[actual, predicted]++;
            }

            return cm;
        }

        public VectorN SubVector( IReadOnlyList<int> indices)
        {
           
            if (indices is null) throw new ArgumentNullException(nameof(indices));

            var src = Values;
            var dst = new double[indices.Count];

            for (int i = 0; i < indices.Count; i++)
            {
                int srcIndex = indices[i];
                if ((uint)srcIndex >= (uint)Length)
                    throw new ArgumentOutOfRangeException(nameof(indices), $"Index {srcIndex} is out of range [0, {Length - 1}].");

                dst[i] = src[srcIndex];
            }

            return new VectorN(dst);

        }


        public override string ToString()
        {
            return "[" + string.Join(", ", Values) + "]";
        }
    }
}
