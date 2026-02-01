using CSharpNumerics.Objects;
using Numerics.Objects;


namespace CSharpNumerics.ML.Selector.Interfaces;

public interface ISelector
{
    Matrix FitTransform(Matrix X, VectorN y);
    Matrix Transform(Matrix X);

    ISelector Clone();
}
