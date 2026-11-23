using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.Models.Interfaces;

/// <summary>
/// A classifier that can say how sure it is: one probability per class for every row.
/// </summary>
public interface IProbabilisticClassifier : IClassificationModel
{
    /// <summary>
    /// Class probabilities, one row per row of <paramref name="X"/> and one column per class
    /// (<see cref="IClassificationModel.NumClasses"/>). Each row sums to 1. For a binary
    /// classifier, column 1 is the probability of class 1.
    /// </summary>
    Matrix PredictProbabilities(Matrix X);
}
