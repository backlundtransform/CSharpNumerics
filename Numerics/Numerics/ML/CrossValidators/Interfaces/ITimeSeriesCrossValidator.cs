using CSharpNumerics.ML.CrossValidators.Grouping;
using CSharpNumerics.ML.CrossValidators.Result;
using Numerics.Models;


namespace CSharpNumerics.ML.CrossValidators.Interfaces
{
    public interface ITimeSeriesCrossValidator
    {
        CrossValidationResult Run(TimeSeries ts, string target, ITimeGrouping grouping = null);
    }
}
