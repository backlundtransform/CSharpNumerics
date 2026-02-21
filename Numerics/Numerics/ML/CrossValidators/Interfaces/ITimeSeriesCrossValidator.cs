using CSharpNumerics.ML.CrossValidators.Grouping;
using CSharpNumerics.ML.CrossValidators.Result;
using CSharpNumerics.Statistics.Data;


namespace CSharpNumerics.ML.CrossValidators.Interfaces
{
    public interface ITimeSeriesCrossValidator
    {
        CrossValidationResult Run(TimeSeries ts, string target, ITimeGrouping grouping = null);
    }
}
