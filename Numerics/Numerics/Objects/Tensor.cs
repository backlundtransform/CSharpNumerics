
using System;
using System.Collections.Generic;
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

        public static Tensor operator +(Tensor a, Tensor b) => a.GetResult(b, 1);


        public Tensor GetResult(Tensor b, int sign, double multiplier = 1)
        {

            var result = Array.CreateInstance(typeof(double), shape);

            var indexArray= new int[dimension];
            var k = 0;
            foreach (var dim in shape.Reverse())
            {

              for (var i = dim-1; i > 0; i -= 1)
                 {
                    for (var j = k; j < dim; j++)
                    {
                        indexArray[i] = j;

                        var value = b.values.GetValue(indexArray);
                    }
                    k++;
                }

             
            }

            return new Tensor(result);

        }

    }

  
}
