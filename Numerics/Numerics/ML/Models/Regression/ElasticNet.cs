using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;


namespace CSharpNumerics.ML.Models.Regression;

public class ElasticNet :IHasHyperparameters, IRegressionModel
{
    public double Lambda { get; set; } = 1.0;
    public double L1Ratio { get; set; } = 0.5;
    public int MaxIterations { get; set; } = 1000;

    private VectorN _weights;

    public void Fit(Matrix X, VectorN y)
    {
        var Xb = X.WithBiasColumn();
        int n = Xb.columnLength;

        _weights = new VectorN(new double[n]);

        for (int iter = 0; iter < MaxIterations; iter++)
        {
            for (int j = 0; j < n; j++)
            {
                double rho = 0.0;
                double norm = 0.0;

                for (int i = 0; i < Xb.rowLength; i++)
                {
                    double pred = 0.0;
                    for (int k = 0; k < n; k++)
                        if (k != j)
                            pred += Xb.values[i, k] * _weights[k];

                    rho += Xb.values[i, j] * (y[i] - pred);
                    norm += Xb.values[i, j] * Xb.values[i, j];
                }

                if (j == 0)
                {
                    _weights[j] = rho / norm;
                }
                else
                {
                    double l1 = Lambda * L1Ratio;
                    double l2 = Lambda * (1 - L1Ratio);

                    _weights[j] =
                        SoftThreshold(rho / norm, l1) /
                        (1 + l2);
                }
            }
        }
    }

    public VectorN Predict(Matrix X)
        => X.WithBiasColumn() * _weights;

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters.TryGetValue("Lambda", out var l))
            Lambda = Convert.ToDouble(l);
        if (parameters.TryGetValue("L1Ratio", out var r))
            L1Ratio = Convert.ToDouble(r);
    }

    private static double SoftThreshold(double z, double g)
        => Math.Sign(z) * Math.Max(Math.Abs(z) - g, 0.0);

    public IModel Clone()
    {
        return new ElasticNet
        {
            Lambda = Lambda,
            L1Ratio = L1Ratio,
            MaxIterations = MaxIterations
        };
    }

}
