using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Models.Regression;

public class KernelSVR : IModel,IHasHyperparameters, IRegressionModel
{
    public double C { get; set; } = 1.0;
    public double Epsilon { get; set; } = 0.1;

    public KernelType Kernel { get; set; } = KernelType.RBF;
    public double Gamma { get; set; } = 0.5;
    public int Degree { get; set; } = 3;


    private VectorN _alpha;
    private VectorN _supportY;
    private Matrix _supportX;

    public void Fit(Matrix X, VectorN y)
    {
        int n = X.rowLength;
        _alpha = new VectorN(n);

        double lr = 0.01;
        int epochs = 500;

        for (int epoch = 0; epoch < epochs; epoch++)
        {
            for (int i = 0; i < n; i++)
            {
                double prediction = 0;
                for (int j = 0; j < n; j++)
                    prediction += _alpha[j] * KernelFunc(
                        X.RowSlice(j), X.RowSlice(i));

                double error = prediction - y[i];

                if (Math.Abs(error) > Epsilon)
                {
                    _alpha[i] -= lr * (
                        Math.Sign(error) * C
                    );
                }
            }
        }

        _supportX = X;
        _supportY = y;
    }

    public VectorN Predict(Matrix X)
    {
        int n = X.rowLength;
        var preds = new double[n];

        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < _supportX.rowLength; j++)
            {
                sum += _alpha[j] *
                       KernelFunc(_supportX.RowSlice(j), X.RowSlice(i));
            }
            preds[i] = sum;
        }

        return new VectorN(preds);
    }

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters.TryGetValue("C", out var c)) C = (double)c;
        if (parameters.TryGetValue("Epsilon", out var e)) Epsilon = (double)e;
        if (parameters.TryGetValue("Gamma", out var g)) Gamma = (double)g;
        if (parameters.TryGetValue("Degree", out var d)) Degree = (int)d;
        if (parameters.TryGetValue("Kernel", out var k)) Kernel = (KernelType)k;
    }
    private double KernelFunc(VectorN x1, VectorN x2)
    {
        return Kernel switch
        {
            KernelType.Linear => x1.Dot(x2),

            KernelType.RBF =>
                Math.Exp(-Gamma * Math.Sqrt((x1 - x2).Norm())),


            KernelType.Polynomial =>
                Math.Pow(Gamma * x1.Dot(x2) + 1, Degree),

            _ => throw new NotSupportedException()
        };
    }
}

