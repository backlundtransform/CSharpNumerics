using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace CSharpNumerics.ML.Models.Classification;

/// <summary>
/// Random forest of <see cref="DecisionTree"/>s, each trained on a bootstrap sample.
/// </summary>
/// <remarks>
/// The bootstrap is drawn as row counts, so trees train on the original matrix without copies,
/// and the trees are trained in parallel. Probabilities are the mean of the trees' leaf class
/// distributions; <see cref="Predict"/> returns the most probable class. Setting
/// <see cref="MaxFeatures"/> (commonly √features) makes the trees more different from each other,
/// which usually helps. Training is repeatable for a given <see cref="Seed"/>.
/// </remarks>
public class RandomForest : IProbabilisticClassifier, ISampleWeightedModel, IHasHyperparameters
{
    private DecisionTree[] _trees = Array.Empty<DecisionTree>();

    public int NumClasses { get; private set; }

    public int NumTrees { get; private set; } = 50;
    public int MaxDepth { get; private set; } = 10;
    public int MinSamplesSplit { get; private set; } = 2;

    /// <summary>Fewest rows a leaf may hold. Default 1.</summary>
    public int MinSamplesLeaf { get; private set; } = 1;

    /// <summary>Features tried at each split, chosen at random per split; 0 (default) tries all.</summary>
    public int MaxFeatures { get; private set; } = 0;

    /// <summary>Seed for the bootstrap samples and feature choices. Default 123.</summary>
    public int Seed { get; private set; } = 123;

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("NumTrees"))
            NumTrees = Convert.ToInt32(parameters["NumTrees"]);

        if (parameters.ContainsKey("MaxDepth"))
            MaxDepth = Convert.ToInt32(parameters["MaxDepth"]);

        if (parameters.ContainsKey("MinSamplesSplit"))
            MinSamplesSplit = Convert.ToInt32(parameters["MinSamplesSplit"]);

        if (parameters.ContainsKey("MinSamplesLeaf"))
            MinSamplesLeaf = Math.Max(1, Convert.ToInt32(parameters["MinSamplesLeaf"]));

        if (parameters.ContainsKey("MaxFeatures"))
            MaxFeatures = Math.Max(0, Convert.ToInt32(parameters["MaxFeatures"]));

        if (parameters.ContainsKey("Seed"))
            Seed = Convert.ToInt32(parameters["Seed"]);
    }

    public void Fit(Matrix X, VectorN y) => FitWeighted(X, y, null);

    public void Fit(Matrix X, VectorN y, VectorN sampleWeights) => FitWeighted(X, y, sampleWeights);

    private void FitWeighted(Matrix X, VectorN y, VectorN? sampleWeights)
    {
        if (X.rowLength != y.Length)
            throw new ArgumentException("X.Rows must match y.Length.");
        if (y.Length == 0)
            throw new ArgumentException("At least one training row is required.");
        if (NumTrees < 1)
            throw new InvalidOperationException("NumTrees must be at least 1.");

        double[] weights = ClassWeights.Validate(sampleWeights, y.Length);
        int[] labels = DecisionTree.Labels(y);
        int n = y.Length;
        NumClasses = labels.Max() + 1;

        // Draw every bootstrap sample up front, in order, so the result does not depend on
        // how the parallel training below is scheduled.
        var random = new Random(Seed);
        var treeRows = new int[NumTrees][];
        var treeWeights = new double[NumTrees][];
        for (int t = 0; t < NumTrees; t++)
        {
            var counts = new int[n];
            for (int i = 0; i < n; i++)
                counts[random.Next(n)]++;

            var w = new double[n];
            var rows = new List<int>();
            for (int i = 0; i < n; i++)
            {
                w[i] = counts[i] * weights[i];
                if (w[i] > 0) rows.Add(i);
            }
            treeRows[t] = rows.ToArray();
            treeWeights[t] = w;
        }

        var trees = new DecisionTree[NumTrees];
        int numClasses = NumClasses;
        Parallel.For(0, NumTrees, t =>
        {
            var tree = new DecisionTree();
            tree.SetHyperParameters(new()
            {
                ["MaxDepth"] = MaxDepth,
                ["MinSamplesSplit"] = MinSamplesSplit,
                ["MinSamplesLeaf"] = MinSamplesLeaf,
                ["MaxFeatures"] = MaxFeatures
            });
            tree.FitRows(X, labels, treeWeights[t], treeRows[t], numClasses, new Random(unchecked(Seed * 7919 + t)));
            trees[t] = tree;
        });
        _trees = trees;
    }

    public VectorN Predict(Matrix X)
    {
        var probabilities = PredictProbabilities(X);
        var preds = new VectorN(X.rowLength);

        for (int i = 0; i < X.rowLength; i++)
        {
            int best = 0;
            for (int c = 1; c < NumClasses; c++)
                if (probabilities.values[i, c] > probabilities.values[i, best])
                    best = c;
            preds[i] = best;
        }

        return preds;
    }

    public Matrix PredictProbabilities(Matrix X)
    {
        if (_trees.Length == 0)
            throw new InvalidOperationException("Model has not been fitted.");

        var result = new Matrix(X.rowLength, NumClasses);
        for (int i = 0; i < X.rowLength; i++)
        {
            foreach (var tree in _trees)
            {
                double[] p = tree.ProbabilitiesRow(X, i);
                for (int c = 0; c < NumClasses; c++)
                    result.values[i, c] += p[c];
            }

            for (int c = 0; c < NumClasses; c++)
                result.values[i, c] /= _trees.Length;
        }

        return result;
    }

    public IModel Clone()
    {
        var clone = new RandomForest();
        clone.SetHyperParameters(new Dictionary<string, object>
        {
            ["NumTrees"] = NumTrees,
            ["MaxDepth"] = MaxDepth,
            ["MinSamplesSplit"] = MinSamplesSplit,
            ["MinSamplesLeaf"] = MinSamplesLeaf,
            ["MaxFeatures"] = MaxFeatures,
            ["Seed"] = Seed
        });
        return clone;
    }
}
