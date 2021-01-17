
using System.Collections;

namespace Numerics.Objects
{
  public struct Tensor
    {

        public dynamic values;
       

        public Tensor(double[,,] tensor)
        {
            
            values = tensor;
        }

        public Tensor(double[,,,] tensor)
        {
          
            values = tensor;
        }
        public Tensor(double[,,,,] tensor)
        {
            values = tensor;
        }

    }

  
}
