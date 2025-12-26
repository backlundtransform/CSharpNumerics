using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.ML.Selector.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Emit;
using System.Text;

namespace CSharpNumerics.ML;

public class PipelineGrid
{
    private readonly List<ModelGrid> _models = new();

    public PipelineGrid AddModel<T>(Action<ModelBuilder> setup)
        where T : IModel, new()
    {
        var mg = new ModelGrid
        {
            Model = new SearchGrid<IModel> { Factory = () => new T() }
        };

        var builder = new ModelBuilder(mg);
        setup(builder);

        _models.Add(mg);
        return this;
    }

    public IEnumerable<Pipeline> Expand()
    {
        foreach (var mg in _models)
        {
           

            var scalers = mg.Scalers.Any()
                ? mg.Scalers
                : new List<SearchGrid<IScaler>> { new() { Factory = () => null } };

            var selectors = mg.Selectors.Any()
                ? mg.Selectors
                : new List<SearchGrid<ISelector>> { new() { Factory = () => null } };

            foreach (var scalerGrid in scalers)
                foreach (var selectorGrid in selectors)
                {
                    

                    foreach (var modelParams in CartesianProduct(mg.Model.Parameters))
                        foreach (var scalerParams in CartesianProduct(scalerGrid.Parameters))
                            foreach (var selectorParams in CartesianProduct(selectorGrid.Parameters))
                            {
                                yield return new Pipeline(
                                    model: mg.Model.Factory(),
                                    modelParams: modelParams,
                                    scaler: scalerGrid.Factory(),
                                    scalerParams: scalerParams,
                                    selector: selectorGrid.Factory(),
                                    selectorParams: selectorParams
                                );
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
   public class ModelGrid
    {
        public SearchGrid<IModel> Model { get; set; }

        public List<SearchGrid<IScaler>> Scalers { get; } = new();
        public List<SearchGrid<ISelector>> Selectors { get; } = new();
    }
    public class ModelBuilder
    {
        private readonly ModelGrid _grid;

        public ModelBuilder(ModelGrid grid)
        {
            _grid = grid;
        }

        public ModelBuilder Add(string name, params object[] values)
        {
            _grid.Model.Add(name, values);
            return this;
        }

        public ModelBuilder AddScaler<T>(Action<SearchGrid<IScaler>> setup)
            where T : IScaler, new()
        {
            var g = new SearchGrid<IScaler> { Factory = () => new T() };
            setup(g);
            _grid.Scalers.Add(g);
            return this;
        }

        public ModelBuilder AddSelector<T>(Action<SearchGrid<ISelector>> setup)
            where T : ISelector, new()
        {
            var g = new SearchGrid<ISelector> { Factory = () => new T() };
            setup(g);
            _grid.Selectors.Add(g);
            return this;
        }
    }
}