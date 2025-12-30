
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Models.Regression;

public class LinearSVR : IRegressionModel, IHasHyperparameters
{
    public double C { get; set; } = 1.0;
    public double Epsilon { get; set; } = 0.1;
    public double LearningRate { get; set; } = 0.01;
    public int Epochs { get; set; } = 1000;

    private VectorN _weights;
    private double _bias;
    public void Fit(Matrix X, VectorN y)
    {
        int nSamples = X.rowLength;
        int nFeatures = X.columnLength;

        _weights = new VectorN(new double[nFeatures]);
        _bias = 0;

        for (int epoch = 0; epoch < Epochs; epoch++)
        {
            for (int i = 0; i < nSamples; i++)
            {
                var xi = X.RowSlice(i);
                double pred = _weights.Dot(xi) + _bias;
                double error = pred - y[i];

                if (Math.Abs(error) <= Epsilon)
                {
                    _weights -= LearningRate * _weights;
                }
                else
                {
                    double sign = Math.Sign(error);
                    _weights -= LearningRate * (_weights + C * sign * xi);
                    _bias -= LearningRate * C * sign;
                }
            }
        }
    }


    public VectorN Predict(Matrix X)
    {
        int n = X.rowLength;
        double[] preds = new double[n];

        for (int i = 0; i < n; i++)
            preds[i] = _weights.Dot(X.RowSlice(i)) + _bias;

        return new VectorN(preds);
    }


    public void SetHyperParameters(Dictionary<string, object> p)
    {
        if (p.TryGetValue("C", out var c)) C = (double)c;
        if (p.TryGetValue("Epsilon", out var e)) Epsilon = (double)e;
        if (p.TryGetValue("LearningRate", out var lr)) LearningRate = (double)lr;
        if (p.TryGetValue("Epochs", out var it)) Epochs = (int)it;
    }   
}
