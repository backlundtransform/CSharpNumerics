using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CSharpNumerics.ML.Models.Regression
{
    public class Linear : IRegressionModel
    {
     

        public void Fit(Matrix X, VectorN y)
        {
            throw new NotImplementedException();
        }

        public VectorN Predict(Matrix X)
        {
            throw new NotImplementedException();
        }
    }
}
 