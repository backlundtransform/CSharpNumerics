using CSharpNumerics.ML.Models.Interfaces;

namespace CSharpNumerics.ML.Sequence.Interfaces;

public interface ISequenceModel : IModel
{
    int TimeSteps { get; set; }

    int Features { get; set; }
}