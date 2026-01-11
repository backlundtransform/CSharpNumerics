using System;
using System.Linq;

namespace Numerics.Objects
{
    public struct Tensor
    {
        public double[] Values;
        public int[] Shape;
        public int Dimension;


        public Tensor(params int[] shape)
        {
            Dimension = shape.Length;
            Shape = shape.ToArray(); 
            int totalLength = 1;
            foreach (var s in shape) totalLength *= s;
            Values = new double[totalLength];
        }

        public Tensor(double[] values, params int[] shape)
        {
            int totalLength = 1;
            foreach (var s in shape) totalLength *= s;
            if (values.Length != totalLength)
                throw new ArgumentException("Values length does not match shape");

            Values = values;
            Shape = shape.ToArray();
            Dimension = shape.Length;
        }

 
        private int GetFlatIndex(int[] indices)
        {
            if (indices.Length != Dimension)
                throw new ArgumentException("Incorrect number of indices");

            int index = 0;
            int multiplier = 1;
            for (int i = Dimension - 1; i >= 0; i--)
            {
                if (indices[i] < 0 || indices[i] >= Shape[i])
                    throw new IndexOutOfRangeException();

                index += indices[i] * multiplier;
                multiplier *= Shape[i];
            }
            return index;
        }

     
        public double this[params int[] indices]
        {
            get => Values[GetFlatIndex(indices)];
            set => Values[GetFlatIndex(indices)] = value;
        }

        private static void CheckShape(Tensor a, Tensor b)
        {
            if (!a.Shape.SequenceEqual(b.Shape))
                throw new ArgumentException("Shape mismatch");
        }

        public static Tensor operator +(Tensor a, Tensor b)
        {
            CheckShape(a, b);
            var result = new Tensor(a.Shape);
            for (int i = 0; i < a.Values.Length; i++)
                result.Values[i] = a.Values[i] + b.Values[i];
            return result;
        }

        public static Tensor operator -(Tensor a, Tensor b)
        {
            CheckShape(a, b);
            var result = new Tensor(a.Shape);
            for (int i = 0; i < a.Values.Length; i++)
                result.Values[i] = a.Values[i] - b.Values[i];
            return result;
        }

        public static Tensor operator *(Tensor a, Tensor b)
        {
            CheckShape(a, b);
            var result = new Tensor(a.Shape);
            for (int i = 0; i < a.Values.Length; i++)
                result.Values[i] = a.Values[i] * b.Values[i];
            return result;
        }

        public static Tensor operator /(Tensor a, Tensor b)
        {
            CheckShape(a, b);
            var result = new Tensor(a.Shape);
            for (int i = 0; i < a.Values.Length; i++)
                result.Values[i] = a.Values[i] / b.Values[i];
            return result;
        }

        // Dot-produkt för lika långa tensorer
        public double Dot(Tensor b)
        {
            if (Values.Length != b.Values.Length)
                throw new ArgumentException("Tensor lengths must match for dot product");

            double sum = 0;
            for (int i = 0; i < Values.Length; i++)
                sum += Values[i] * b.Values[i];
            return sum;
        }

        // Fyll tensor med ett värde
        public void Fill(double value)
        {
            for (int i = 0; i < Values.Length; i++)
                Values[i] = value;
        }
    }
}
