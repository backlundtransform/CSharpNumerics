using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;


namespace CSharpNumerics.ML.Models.Classification;

/// <summary>
/// CART classification tree with Gini impurity.
/// </summary>
/// <remarks>
/// Splits are searched by sorting each feature once per node, so a node with m rows costs
/// O(p · m log m). A split sends rows with <c>x[f] &lt;= threshold</c> left, where the threshold is
/// one of the feature's values. Leaves keep the (weighted) class distribution of their rows, which
/// <see cref="PredictProbabilities"/> returns.
/// </remarks>
public class DecisionTree : IProbabilisticClassifier, ISampleWeightedModel, IHasHyperparameters
{
    private TreeNode _root;
    public int NumClasses { get; private set; }
    public int MaxDepth { get; private set; } = 10;
    public int MinSamplesSplit { get; private set; } = 2;

    /// <summary>Fewest rows a leaf may hold. Larger values give smoother probabilities. Default 1.</summary>
    public int MinSamplesLeaf { get; private set; } = 1;

    /// <summary>Features tried at each split, chosen at random; 0 (default) tries all of them.</summary>
    public int MaxFeatures { get; private set; } = 0;

    /// <summary>Seed for the feature choice when <see cref="MaxFeatures"/> is set.</summary>
    public int Seed { get; private set; } = 123;

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

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

        double[] weights = ClassWeights.Validate(sampleWeights, y.Length);
        int[] labels = Labels(y);
        int[] rows = Enumerable.Range(0, y.Length).Where(i => weights[i] > 0).ToArray();
        FitRows(X, labels, weights, rows, labels.Max() + 1, new Random(Seed));
    }

    /// <summary>
    /// Trains on the rows <paramref name="rows"/> of <paramref name="X"/> with the given weights,
    /// without copying the data. Used by <see cref="RandomForest"/> for its bootstrap samples.
    /// </summary>
    internal void FitRows(Matrix X, int[] labels, double[] weights, int[] rows, int numClasses, Random random)
    {
        NumClasses = numClasses;
        _root = Build(X, labels, weights, rows, 0, random);
    }

    public VectorN Predict(Matrix X)
    {
        VectorN preds = new VectorN(X.rowLength);

        for (int i = 0; i < X.rowLength; i++)
            preds[i] = Leaf(X, i).PredictedClass;

        return preds;
    }

    public Matrix PredictProbabilities(Matrix X)
    {
        var result = new Matrix(X.rowLength, NumClasses);
        for (int i = 0; i < X.rowLength; i++)
        {
            double[] p = Leaf(X, i).Probabilities;
            for (int c = 0; c < NumClasses; c++)
                result.values[i, c] = p[c];
        }
        return result;
    }

    internal int PredictRow(Matrix X, int row) => Leaf(X, row).PredictedClass;

    internal double[] ProbabilitiesRow(Matrix X, int row) => Leaf(X, row).Probabilities;

    private TreeNode Leaf(Matrix X, int row)
    {
        if (_root == null)
            throw new InvalidOperationException("Model has not been fitted.");

        var node = _root;
        while (!node.IsLeaf)
            node = X.values[row, node.FeatureIndex] <= node.Threshold ? node.Left : node.Right;
        return node;
    }

    private TreeNode Build(Matrix X, int[] labels, double[] weights, int[] rows, int depth, Random random)
    {
        double[] counts = new double[NumClasses];
        double total = 0;
        foreach (int i in rows)
        {
            counts[labels[i]] += weights[i];
            total += weights[i];
        }

        bool pure = counts.Count(c => c > 0) <= 1;
        if (depth >= MaxDepth || rows.Length < MinSamplesSplit || rows.Length < 2 * MinSamplesLeaf || pure)
            return LeafNode(counts, total);

        var (feature, threshold) = FindBestSplit(X, labels, weights, rows, counts, total, random);
        if (feature < 0)
            return LeafNode(counts, total);

        var left = new List<int>(rows.Length);
        var right = new List<int>(rows.Length);
        foreach (int i in rows)
            (X.values[i, feature] <= threshold ? left : right).Add(i);

        return new TreeNode
        {
            FeatureIndex = feature,
            Threshold = threshold,
            Left = Build(X, labels, weights, left.ToArray(), depth + 1, random),
            Right = Build(X, labels, weights, right.ToArray(), depth + 1, random)
        };
    }

    private (int Feature, double Threshold) FindBestSplit(Matrix X, int[] labels, double[] weights, int[] rows,
        double[] counts, double total, Random random)
    {
        int m = rows.Length;
        var keys = new double[m];
        var order = new int[m];
        var leftCounts = new double[NumClasses];

        double bestScore = double.MaxValue;
        int bestFeature = -1;
        double bestThreshold = 0;

        foreach (int f in Features(X.columnLength, random))
        {
            for (int k = 0; k < m; k++)
            {
                keys[k] = X.values[rows[k], f];
                order[k] = rows[k];
            }
            Array.Sort(keys, order);
            Array.Clear(leftCounts, 0, leftCounts.Length);
            double leftTotal = 0;

            for (int k = 0; k < m - 1; k++)
            {
                int i = order[k];
                leftCounts[labels[i]] += weights[i];
                leftTotal += weights[i];

                // A threshold must separate distinct values: rows with equal values go the same way.
                // NaN sorts first and compares false, so it can be neither a threshold nor split off.
                if (keys[k] == keys[k + 1] || double.IsNaN(keys[k])) continue;
                int leftRows = k + 1;
                if (leftRows < MinSamplesLeaf || m - leftRows < MinSamplesLeaf) continue;

                double rightTotal = total - leftTotal;
                double giniLeft = 1.0, giniRight = 1.0;
                for (int c = 0; c < NumClasses; c++)
                {
                    double pl = leftCounts[c] / leftTotal;
                    double pr = (counts[c] - leftCounts[c]) / rightTotal;
                    giniLeft -= pl * pl;
                    giniRight -= pr * pr;
                }

                double score = (leftTotal * giniLeft + rightTotal * giniRight) / total;
                if (score < bestScore)
                {
                    bestScore = score;
                    bestFeature = f;
                    bestThreshold = keys[k];
                }
            }
        }

        return (bestFeature, bestThreshold);
    }

    /// <summary>All features in order, or <see cref="MaxFeatures"/> of them at random.</summary>
    private IEnumerable<int> Features(int count, Random random)
    {
        if (MaxFeatures <= 0 || MaxFeatures >= count)
            return Enumerable.Range(0, count);

        var all = Enumerable.Range(0, count).ToArray();
        for (int i = 0; i < MaxFeatures; i++)
        {
            int j = i + random.Next(count - i);
            (all[i], all[j]) = (all[j], all[i]);
        }
        return all.Take(MaxFeatures);
    }

    private static TreeNode LeafNode(double[] counts, double total)
    {
        var probabilities = new double[counts.Length];
        int best = 0;
        for (int c = 0; c < counts.Length; c++)
        {
            probabilities[c] = total > 0 ? counts[c] / total : 1.0 / counts.Length;
            if (counts[c] > counts[best]) best = c;
        }

        return new TreeNode { IsLeaf = true, PredictedClass = best, Probabilities = probabilities };
    }

    /// <summary>Class labels as integers 0, 1, 2, … (the values of y).</summary>
    internal static int[] Labels(VectorN y)
    {
        var labels = new int[y.Length];
        for (int i = 0; i < y.Length; i++)
        {
            double v = y[i];
            if (v < 0 || v != Math.Floor(v))
                throw new ArgumentException($"Class labels must be 0, 1, 2, …; got {v}.");
            labels[i] = (int)v;
        }
        return labels;
    }

    internal class TreeNode
    {
        public bool IsLeaf;
        public int PredictedClass;
        public double[] Probabilities;

        public int FeatureIndex;
        public double Threshold;

        public TreeNode Left;
        public TreeNode Right;
    }

    public IModel Clone()
    {
        var clone = new DecisionTree();
        clone.SetHyperParameters(new Dictionary<string, object>
        {
            ["MaxDepth"] = MaxDepth,
            ["MinSamplesSplit"] = MinSamplesSplit,
            ["MinSamplesLeaf"] = MinSamplesLeaf,
            ["MaxFeatures"] = MaxFeatures,
            ["Seed"] = Seed
        });
        return clone;
    }
}
