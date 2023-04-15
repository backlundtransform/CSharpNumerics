
using System;
using System.Collections.Generic;
using System.Globalization;

namespace Numerics.Objects
{
  public struct Tensor
    {

        public Array values;

        public int dimension;

        public int[] shape;
     

        public Tensor(Array tensor)
        {
            dimension = tensor.Rank;
            shape = new int[dimension];

            for (var i=0;i < dimension; i++)
            {
                shape[i] =tensor.GetLength(i);

            }
          
            values = tensor;
         
        }

        public static Tensor operator +(Tensor a, Tensor b) => a.GetResult(b, (x, y) => x + y);
        public static Tensor operator -(Tensor a, Tensor b) => a.GetResult(b, (x, y) => x - y);
        public static Tensor operator *(Tensor a, Tensor b) => a.GetResult(b, (x, y) => x * y);
        public static Tensor operator /(Tensor a, Tensor b) => a.GetResult(b, (x, y) => x / y);

        public Tensor TensorDot(Tensor b) {

            var result = new List<Array>();


            foreach (var index in values)
            {
                if (!double.TryParse(index.ToString(), out double convertA))
                {
                    throw new ArgumentException("Is not numeric");

                }
                result.Add(GetResult(b, (x, y) => convertA * y).values);
            }

            return new Tensor(result.ToArray());

        }
        public Tensor GetResult(Tensor b, Func<double, double, double> operatoru)
        {
            
            var result = Array.CreateInstance(typeof(double), shape);
          
            var lastIndex = dimension- 1;

            var indices = new int[dimension];
            var lower = new int[dimension];
            var upper = new int[dimension];
            for (var i = 0; i < dimension; ++i)
            {
                indices[i] = lower[i] = values.GetLowerBound(i);
                upper[i] = values.GetUpperBound(i);
            }
            
            while (true)
            {

               if (!double.TryParse(values.GetValue(indices).ToString(), out double convertA) || !double.TryParse(b.values.GetValue(indices).ToString(), out double convertB))
                {

                    throw new ArgumentException("Is not numric");

                }

                var element = operatoru(convertA, convertB);
                result.SetValue(element, indices);
             
                
                if(!InnerLoop(lastIndex, indices, upper, lower)){

                    continue;
                }

                if (++indices[0] > upper[0])
                {
                    break;
                }

            }

            return new Tensor(result);

        }
       
        private static bool InnerLoop(int lastIndex, int[] indices, int[] upper, int[] lower)
        {
            for (var i = lastIndex; i > 0; i -= 1)
            {
                if (++indices[i] <= upper[i])
                {

                    return false;
                }

                indices[i] = lower[i];

            }

            return true;
        }

    }

}
