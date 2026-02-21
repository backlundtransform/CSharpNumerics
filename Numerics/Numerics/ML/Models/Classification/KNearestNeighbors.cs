using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML.Models.Classification;

public class KNearestNeighbors :
    IClassificationModel,
    IHasHyperparameters
{
    private Matrix _X;
    private VectorN _y;
    
    public int K { get; private set; } = 5;


    public int NumClasses { get; private set; }

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters != null && parameters.ContainsKey("K"))
            K = Convert.ToInt32(parameters["K"]);
    }

    public void Fit(Matrix X, VectorN y)
    {
        _X = X;
        _y = y;
        NumClasses = (int)y.Values.Max() + 1;
    }

    public VectorN Predict(Matrix X)
    {
        VectorN preds = new VectorN(X.rowLength);

        for (int i = 0; i < X.rowLength; i++)
        {
            var distances = new List<(double d, int cls)>();

            for (int j = 0; j < _X.rowLength; j++)
            {
                double dist = 0;
                for (int k = 0; k < X.columnLength; k++)
                    dist += Math.Pow(X.values[i, k] - _X.values[j, k], 2);

                distances.Add((Math.Sqrt(dist), (int)_y[j]));
            }

            var vote = distances
                .OrderBy(t => t.d)
                .Take(K)
                .GroupBy(t => t.cls)
                .OrderByDescending(g => g.Count())
                .First().Key;

            preds[i] = vote;
        }

        return preds;
    }

    public IModel Clone()
    {
        var clone = new KNearestNeighbors();
        clone.SetHyperParameters(new Dictionary<string, object>
        {
            ["K"] = K
        });
        return clone;
    }
}
