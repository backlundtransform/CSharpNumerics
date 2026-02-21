using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Scalers.Interfaces;

public interface IScaler
{
    Matrix FitTransform(Matrix X);
    Matrix Transform(Matrix X);

    IScaler Clone();
}
