
using CSharpNumerics.Numerics.Objects;


namespace CSharpNumerics.ML.Models.Interfaces;


public interface IModel
{
    void Fit(Matrix X, VectorN y);
    VectorN Predict(Matrix X);

    IModel Clone();
  
}
