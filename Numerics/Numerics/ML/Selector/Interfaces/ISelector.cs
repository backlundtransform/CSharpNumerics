using Numerics.Objects;


namespace CSharpNumerics.ML.Selector.Interfaces;

public interface ISelector
{
    Matrix FitTransform(Matrix X, double[] y);
    Matrix Transform(Matrix X);
}
