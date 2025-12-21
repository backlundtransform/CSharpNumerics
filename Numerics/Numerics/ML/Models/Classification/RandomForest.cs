using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Models.Classification;

public class RandomForest : IClassificationModel, IHasHyperparameters
{
    private List<DecisionTree> _trees = new();
    private Random _rnd = new Random(123);

    public int NumClasses { get; private set; }

    public int NumTrees { get; private set; } = 50;
    public int MaxDepth { get; private set; } = 10;
    public int MinSamplesSplit { get; private set; } = 2;

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("NumTrees"))
            NumTrees = Convert.ToInt32(parameters["NumTrees"]);

        if (parameters.ContainsKey("MaxDepth"))
            MaxDepth = Convert.ToInt32(parameters["MaxDepth"]);

        if (parameters.ContainsKey("MinSamplesSplit"))
            MinSamplesSplit = Convert.ToInt32(parameters["MinSamplesSplit"]);
    }

    public void Fit(Matrix X, VectorN y)
    {
        NumClasses = (int)y.Values.Max() + 1;
        _trees.Clear();

        for (int i = 0; i < NumTrees; i++)
        {
            var (Xb, yb) = BootstrapSample(X, y);

            var tree = new DecisionTree();
            tree.SetHyperParameters(new()
            {
                ["MaxDepth"] = MaxDepth,
                ["MinSamplesSplit"] = MinSamplesSplit
            });

            tree.Fit(Xb, yb);
            _trees.Add(tree);
        }
    }

    public VectorN Predict(Matrix X)
    {
        VectorN preds = new VectorN(X.rowLength);

        for (int i = 0; i < X.rowLength; i++)
        {
            int[] votes = new int[NumClasses];

            foreach (var tree in _trees)
            {
                int cls = (int)tree.PredictRow(X, i);
                votes[cls]++;
            }

            preds[i] = ArgMax(votes);
        }

        return preds;
    }
    private (Matrix, VectorN) BootstrapSample(Matrix X, VectorN y)
    {
        int n = X.rowLength;
        int p = X.columnLength;

        double[,] Xb = new double[n, p];
        double[] yb = new double[n];

        for (int i = 0; i < n; i++)
        {
            int idx = _rnd.Next(n);

            for (int j = 0; j < p; j++)
                Xb[i, j] = X.values[idx, j];

            yb[i] = y[idx];
        }

        return (new Matrix(Xb), new VectorN(yb));
    }

    private int ArgMax(int[] votes)
    {
        int best = 0;
        for (int i = 1; i < votes.Length; i++)
            if (votes[i] > votes[best])
                best = i;
        return best;
    }
}
