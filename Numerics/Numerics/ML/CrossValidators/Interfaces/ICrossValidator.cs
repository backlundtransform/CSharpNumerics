using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.CrossValidators.Interfaces;

public interface ICrossValidator
{
    CrossValidationResult Run(Matrix X, VectorN y);
}
