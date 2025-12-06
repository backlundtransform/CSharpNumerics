using System;


namespace CSharpNumerics.Objects
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

     
        public override string ToString()
        {
            return "[" + string.Join(", ", Values) + "]";
        }
    }
}
