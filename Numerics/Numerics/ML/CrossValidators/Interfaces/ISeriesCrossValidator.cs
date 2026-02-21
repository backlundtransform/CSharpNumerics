using CSharpNumerics.ML.CrossValidators.Grouping;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.Statistics.Data;


namespace CSharpNumerics.ML.CrossValidators.Interfaces;

public interface ISeriesCrossValidator
{
    CrossValidationResult Run(Series ts, string target);
}
