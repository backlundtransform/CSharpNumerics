using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.ML.Selector.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML;

public class PipelineGrid
{
    private readonly List<SearchGrid<IScaler>> _scalers = new();
    private readonly List<SearchGrid<ISelector>> _selectors = new();
    private readonly List<SearchGrid<IModel>> _models = new();

    public PipelineGrid AddModel<T>(Action<SearchGrid<IModel>> setup) where T : IModel, new()
    {
        var grid = new SearchGrid<IModel> { Factory = () => new T() };
        setup(grid);
        _models.Add(grid);
        return this;
    }

    public PipelineGrid AddScaler<T>(Action<SearchGrid<IScaler>> setup) where T : IScaler, new()
    {
        var grid = new SearchGrid<IScaler> { Factory = () => new T() };
        setup(grid);
        _scalers.Add(grid);
        return this;
    }
    public PipelineGrid AddSelector<T>(Action<SearchGrid<ISelector>> setup) where T : ISelector, new()
    {
        var grid = new SearchGrid<ISelector> { Factory = () => new T() };
        setup(grid);
        _selectors.Add(grid);
        return this;
    }

    public IEnumerable<Pipeline> Expand()
    {
       
        var scalerFactories = _scalers.Any() ? _scalers : new List<SearchGrid<IScaler>> {  };
        var selectorGrids = _selectors.Any() ? _selectors : new List<SearchGrid<ISelector>> { new SearchGrid<ISelector> { Factory = () => null } };

        foreach (var scalerFactory in scalerFactories)
        {
            foreach (var selectorGrid in selectorGrids)
            {
                var selectorParamsList = CartesianProduct(selectorGrid.Parameters);

                foreach (var selParams in selectorParamsList)
                {
                    foreach (var modelGrid in _models)
                    {
                        var modelParamsList = CartesianProduct(modelGrid.Parameters);

                        foreach (var modParams in modelParamsList)
                        {
                            yield return new Pipeline(
                                model: modelGrid.Factory(),
                                modelParams: modParams,
                                scaler: scalerFactory.Factory(),
                                selector: selectorGrid.Factory(),
                                selectorParams: selParams
                            );
                        }
                    }
                }
            }
        }
    }

    private static List<Dictionary<string, object>> CartesianProduct(
       Dictionary<string, object[]> grid)
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

   public class SearchGrid<T>
    {
        public Func<T> Factory { get; set; }
        public Dictionary<string, object[]> Parameters { get; set; } = new();

        public SearchGrid<T> Add(string name, params object[] values)
        {
            Parameters[name] = values;
            return this;
        }
    }
}