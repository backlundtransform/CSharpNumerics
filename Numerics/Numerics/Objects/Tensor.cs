
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

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

        public static Tensor operator +(Tensor a, Tensor b) => a.GetResult(b);


        public Tensor GetResult(Tensor b)
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
                var element = Convert.ToDouble(values.GetValue(indices), CultureInfo.InvariantCulture) + Convert.ToDouble(b.values.GetValue(indices), CultureInfo.InvariantCulture);
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
       
        private bool InnerLoop(int lastIndex, int[] indices, int[] upper, int[] lower)
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
