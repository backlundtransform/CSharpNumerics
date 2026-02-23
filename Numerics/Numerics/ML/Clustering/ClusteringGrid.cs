using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Scalers.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Clustering;

/// <summary>
/// Grid search over clustering algorithms and hyperparameters.
/// Mirrors <see cref="PipelineGrid"/> but for unsupervised models.
/// </summary>
public class ClusteringGrid
{
    private readonly List<ClusteringModelGrid> _models = new();

    public ClusteringGrid AddModel<T>(Action<ClusteringModelBuilder> setup)
        where T : IClusteringModel, new()
    {
        var mg = new ClusteringModelGrid
        {
            Model = new SearchGrid<IClusteringModel> { Factory = () => new T() }
        };

        var builder = new ClusteringModelBuilder(mg);
        setup(builder);

        _models.Add(mg);
        return this;
    }

    /// <summary>
    /// Expands all model × hyperparameter × scaler combinations
    /// into individual <see cref="ClusteringPipeline"/> instances.
    /// </summary>
    public IEnumerable<ClusteringPipeline> Expand()
    {
        foreach (var mg in _models)
        {
            var scalers = mg.Scalers.Any()
                ? mg.Scalers
                : new List<SearchGrid<IScaler>> { new() { Factory = () => null } };

            foreach (var scalerGrid in scalers)
            {
                foreach (var modelParams in CartesianProduct(mg.Model.Parameters))
                    foreach (var scalerParams in CartesianProduct(scalerGrid.Parameters))
                    {
                        yield return new ClusteringPipeline(
                            model: mg.Model.Factory(),
                            modelParams: modelParams,
                            scaler: scalerGrid.Factory(),
                            scalerParams: scalerParams
                        );
                    }
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════
    //  Cartesian product (same logic as PipelineGrid)
    // ═══════════════════════════════════════════════════════════════

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

    // ═══════════════════════════════════════════════════════════════
    //  Nested types
    // ═══════════════════════════════════════════════════════════════

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

    public class ClusteringModelGrid
    {
        public SearchGrid<IClusteringModel> Model { get; set; }
        public List<SearchGrid<IScaler>> Scalers { get; } = new();
    }

    public class ClusteringModelBuilder
    {
        private readonly ClusteringModelGrid _grid;

        public ClusteringModelBuilder(ClusteringModelGrid grid)
        {
            _grid = grid;
        }

        public ClusteringModelBuilder Add(string name, params object[] values)
        {
            _grid.Model.Add(name, values);
            return this;
        }

        public ClusteringModelBuilder AddScaler<T>(Action<SearchGrid<IScaler>> setup)
            where T : IScaler, new()
        {
            var g = new SearchGrid<IScaler> { Factory = () => new T() };
            setup(g);
            _grid.Scalers.Add(g);
            return this;
        }
    }
}
