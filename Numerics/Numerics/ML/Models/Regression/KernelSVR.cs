using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Objects;
using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Models.Regression;

public class KernelSVR : IHasHyperparameters, IRegressionModel
{
    public double C { get; set; } = 1.0;
    public double Epsilon { get; set; } = 0.1;
    public double LearningRate { get; set; } = 0.01;
    public int Epochs { get; set; } = 500;

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

        for (int epoch = 0; epoch < Epochs; epoch++)
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
                    _alpha[i] -= LearningRate * (
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
        if (parameters.TryGetValue("C", out var c)) C = Convert.ToDouble(c);
        if (parameters.TryGetValue("Epsilon", out var e)) Epsilon = Convert.ToDouble(e);
        if (parameters.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (parameters.TryGetValue("Degree", out var d)) Degree = Convert.ToInt32(d);
        if (parameters.TryGetValue("Kernel", out var k)) Kernel = (KernelType)k;
        if (parameters.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (parameters.TryGetValue("Epochs", out var ep)) Epochs = Convert.ToInt32(ep);
    }
    private double KernelFunc(VectorN x1, VectorN x2)
    {
        return Kernel switch
        {
            KernelType.Linear => x1.Dot(x2),

            KernelType.RBF =>
                Math.Exp(-Gamma * Math.Pow((x1 - x2).Norm(), 2)),


            KernelType.Polynomial =>
                Math.Pow(Gamma * x1.Dot(x2) + 1, Degree),

            _ => throw new NotSupportedException()
        };
    }

    public IModel Clone()
    {
        return new KernelSVR
        {
            C = C,
            Epsilon = Epsilon,
            Kernel = Kernel,
            Gamma = Gamma,
            Degree = Degree,
            LearningRate = LearningRate,
            Epochs = Epochs
        };
    }
}

