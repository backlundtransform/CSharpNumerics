using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;


namespace CSharpNumerics.ML.Models.Classification;

public class DecisionTree : IClassificationModel, IHasHyperparameters
{
    private TreeNode _root;
    public int NumClasses { get; private set; }
    public int MaxDepth { get; private set; } = 10;
    public int MinSamplesSplit { get; private set; } = 2;

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        if (parameters.ContainsKey("MaxDepth"))
            MaxDepth = Convert.ToInt32(parameters["MaxDepth"]);

        if (parameters.ContainsKey("MinSamplesSplit"))
            MinSamplesSplit = Convert.ToInt32(parameters["MinSamplesSplit"]);
    }

    public void Fit(Matrix X, VectorN y)
    {
        NumClasses = (int)y.Values.Max() + 1;
        _root = BuildTree(X, y, depth: 0);
    }

    public VectorN Predict(Matrix X)
    {
        VectorN preds = new VectorN(X.rowLength);

        for (int i = 0; i < X.rowLength; i++)
            preds[i] = PredictSample(X, i, _root);

        return preds;
    }

    private TreeNode BuildTree(Matrix X, VectorN y, int depth)
    {
        
        if (depth >= MaxDepth || X.rowLength < MinSamplesSplit || IsPure(y))
        {
            return new TreeNode
            {
                IsLeaf = true,
                PredictedClass = MajorityClass(y)
            };
        }

        var (bestFeature, bestThreshold) = FindBestSplit(X, y);

        if (bestFeature == -1)
        {
            return new TreeNode
            {
                IsLeaf = true,
                PredictedClass = MajorityClass(y)
            };
        }

        var (Xl, yl, Xr, yr) = Split(X, y, bestFeature, bestThreshold);

        // If either split is empty, return a leaf node
        if (yl.Length == 0 || yr.Length == 0)
        {
            return new TreeNode
            {
                IsLeaf = true,
                PredictedClass = MajorityClass(y)
            };
        }

        return new TreeNode
        {
            FeatureIndex = bestFeature,
            Threshold = bestThreshold,
            Left = BuildTree(Xl, new VectorN(yl), depth + 1),
            Right = BuildTree(Xr, new VectorN(yr), depth + 1)
        };
    }
    private double Gini(VectorN y)
    {
        double gini = 1.0;
        for (int c = 0; c < NumClasses; c++)
        {
            double p = y.Values.Count(v => (int)v == c) / (double)y.Length;
            gini -= p * p;
        }
        return gini;
    }

    private (int, double) FindBestSplit(Matrix X, VectorN y)
    {
        double bestGini = double.MaxValue;
        int bestFeature = -1;
        double bestThreshold = 0;

        for (int f = 0; f < X.columnLength; f++)
        {
            foreach (double t in X.ColumnSlice(f).Values.Distinct())
            {
                var (Xl, yl, Xr, yr) = Split(X, y, f, t);
                if (yl.Length == 0 || yr.Length == 0) continue;

                double g =
                    (yl.Length * Gini(new VectorN(yl)) + yr.Length * Gini(new VectorN(yr))) / y.Length;

                if (g < bestGini)
                {
                    bestGini = g;
                    bestFeature = f;
                    bestThreshold = t;
                }
            }
        }

        return (bestFeature, bestThreshold);
    }
    private int PredictSample(Matrix X, int row, TreeNode node)
    {
        if (node.IsLeaf)
            return node.PredictedClass;

        if (X.values[row, node.FeatureIndex] <= node.Threshold)
            return PredictSample(X, row, node.Left);
        else
            return PredictSample(X, row, node.Right);
    }
    private bool IsPure(VectorN y)
    {
        int first = (int)y.Values[0];

        for (int i = 1; i < y.Values.Length; i++)
            if ((int)y.Values[i] != first)
                return false;

        return true;
    }
    private int MajorityClass(VectorN y)
    {
        int[] counts = new int[NumClasses];

        for (int i = 0; i < y.Values.Length; i++)
            counts[(int)y.Values[i]]++;

        int maxClass = 0;
        for (int c = 1; c < NumClasses; c++)
            if (counts[c] > counts[maxClass])
                maxClass = c;

        return maxClass;
    }

    private (Matrix Xl, double[] yl, Matrix Xr, double[] yr)
    Split(Matrix X, VectorN y, int feature, double threshold)
    {
        int leftCount = 0;
        int rightCount = 0;

        for (int i = 0; i < X.rowLength; i++)
        {
            if (X.values[i, feature] <= threshold)
                leftCount++;
            else
                rightCount++;
        }

        double[,] Xl = new double[leftCount, X.columnLength];
        double[,] Xr = new double[rightCount, X.columnLength];
        double[] yl = new double[leftCount];
        double[] yr = new double[rightCount];

        int li = 0, ri = 0;

        for (int i = 0; i < X.rowLength; i++)
        {
            if (X.values[i, feature] <= threshold)
            {
                for (int j = 0; j < X.columnLength; j++)
                    Xl[li, j] = X.values[i, j];
                yl[li++] = y[i];
            }
            else
            {
                for (int j = 0; j < X.columnLength; j++)
                    Xr[ri, j] = X.values[i, j];
                yr[ri++] = y[i];
            }
        }

        return (new Matrix(Xl), yl, new Matrix(Xr), yr);
    }
    internal int PredictRow(Matrix X, int row)
    {
        return PredictSample(X, row, _root);
    }
    internal class TreeNode
    {
        public bool IsLeaf;
        public int PredictedClass;

        public int FeatureIndex;
        public double Threshold;

        public TreeNode Left;
        public TreeNode Right;
    }
}