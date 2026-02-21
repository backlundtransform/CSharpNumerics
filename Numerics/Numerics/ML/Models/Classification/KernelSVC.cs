using CSharpNumerics.ML.Enums;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML.Models.Classification;

public class KernelSVC :
    IClassificationModel,
    IHasHyperparameters
{
    public double C { get; set; } = 1.0;
    public int Epochs { get; set; } = 1000;
    public double LearningRate { get; set; } = 0.01;

    public KernelType Kernel { get; set; } = KernelType.RBF;
    public double Gamma { get; set; } = 0.5; 
    public int Degree { get; set; } = 3;

    public int NumClasses { get; private set; }

    private Matrix Xtrain;
    private VectorN ytrain;
    private VectorN alpha;

    public void Fit(Matrix X, VectorN y)
    {
        NumClasses = (int)y.Values.Max() + 1;
        Xtrain = X;
        ytrain = new VectorN([.. y.Values.Select(v => v == 0 ? -1 : 1)]);
        int n = X.rowLength;

        alpha = new VectorN(n);

        for (int epoch = 0; epoch < Epochs; epoch++)
        {
            for (int i = 0; i < n; i++)
            {
                double decision = 0.0;

                for (int j = 0; j < n; j++)
                {
                    if (alpha[j] != 0)
                    {
                        decision += alpha[j] * ytrain[j] *
                            KernelFunc(Xtrain.RowSlice(j), Xtrain.RowSlice(i));
                    }

 
                }

                if (ytrain[i] * decision < 1)
                {
                    alpha[i] += LearningRate * C;
                }
            }
        }
    }
    public IModel Clone()
    {
        return new KernelSVC
        {
            C = C,
            Epochs = Epochs,
            LearningRate = LearningRate,
            Kernel = Kernel,
            Gamma = Gamma,
            Degree = Degree
        };
    }

    public VectorN Predict(Matrix X)
    {
        int m = X.rowLength;
        var preds = new double[m];

        for (int i = 0; i < m; i++)
        {
            double sum = 0.0;

            for (int j = 0; j < Xtrain.rowLength; j++)
            {
                if (alpha[j] != 0)
                {
                    sum += alpha[j] * ytrain[j] *
                        KernelFunc(Xtrain.RowSlice(j), X.RowSlice(i));
                }
            }

            preds[i] = sum >= 0 ? 1 : 0;
        }

        return new VectorN(preds);
    }


    double KernelFunc(VectorN xi, VectorN xj)
    {
        return Kernel switch
        {
            KernelType.Linear =>
                xi.Dot(xj),

            KernelType.Polynomial =>
                Math.Pow(Gamma * xi.Dot(xj) + 1, Degree),

            KernelType.RBF =>
                Math.Exp(-Gamma * Math.Pow((xi - xj).Norm(), 2)),

            _ => throw new NotSupportedException()
        };
    }

    public void SetHyperParameters(Dictionary<string, object> p)
    {

        if (p.TryGetValue("C", out var c)) C = Convert.ToDouble(c);
        if (p.TryGetValue("LearningRate", out var lr)) LearningRate = Convert.ToDouble(lr);
        if (p.TryGetValue("Epochs", out var it)) Epochs = Convert.ToInt32(it);
        if (p.TryGetValue("Kernel", out var k)) Kernel = (KernelType)k;
        if (p.TryGetValue("Gamma", out var g)) Gamma = Convert.ToDouble(g);
        if (p.TryGetValue("Degree", out var d)) Degree = Convert.ToInt32(d);


    }
}

