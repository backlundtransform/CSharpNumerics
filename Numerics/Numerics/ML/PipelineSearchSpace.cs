using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.ML.Selector.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML;

public class PipelineSearchSpace
{
    public List<Func<IScaler>> Scalers { get; } = new();
    public List<(Func<ISelector> factory, Dictionary<string, List<object>> grid)> Selectors { get; } = new();
    public List<(Func<IModel> factory, Dictionary<string, List<object>> grid)> Models { get; } = new();

    public List<Pipeline> Expand()
    {
        var pipelines = new List<Pipeline>();

        foreach (var scaler in Scalers.DefaultIfEmpty(null))
            foreach (var (selectorFactory, selectorGrid) in Selectors.DefaultIfEmpty((null, null)))
                foreach (var (modelFactory, modelGrid) in Models)
                {
                    var modelParamsList = CartesianProduct(modelGrid);

                    foreach (var modelParams in modelParamsList)
                    {
                        pipelines.Add(new Pipeline(
                            model: modelFactory(),
                            modelParams: modelParams,
                            scaler: scaler?.Invoke(),
                            selector: selectorFactory?.Invoke()
                        ));
                    }
                }

        return pipelines;
    }
    private static List<Dictionary<string, object>> CartesianProduct(
       Dictionary<string, List<object>> grid)
    {
        var result = new List<Dictionary<string, object>>
    {
        new Dictionary<string, object>()
    };

        foreach (var kv in grid)
        {
            var temp = new List<Dictionary<string, object>>();

            foreach (var value in kv.Value)
            {
                foreach (var dict in result)
                {
                    var newDict = new Dictionary<string, object>(dict)
                    {
                        [kv.Key] = value
                    };
                    temp.Add(newDict);
                }
            }

            result = temp;
        }

        return result;
    }
}