using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Models.Interfaces;

public interface IModel
{
    void Fit(Matrix X, double[] y);
    double[] Predict(Matrix X);
}
