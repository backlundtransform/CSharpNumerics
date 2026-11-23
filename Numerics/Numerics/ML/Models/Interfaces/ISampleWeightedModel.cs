using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.ML.Models.Interfaces;

/// <summary>
/// A model that can be trained with a weight per sample, for example to give a rare class as
/// much influence as a common one (see <see cref="ClassWeights.Balanced"/>).
/// </summary>
public interface ISampleWeightedModel : IModel
{
    /// <summary>
    /// Trains with <paramref name="sampleWeights"/>: one non-negative weight per row. A weight of
    /// 2 counts the row as if it appeared twice. All weights equal gives the same model as
    /// <see cref="IModel.Fit"/>.
    /// </summary>
    void Fit(Matrix X, VectorN y, VectorN sampleWeights);
}
