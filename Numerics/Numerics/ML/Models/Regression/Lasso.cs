using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.ML.Models.Regression;



public class Lasso: IHasHyperparameters, IRegressionModel
{
    private VectorN _weights;
    private double _alpha = 0.1;
    private int _iterations = 1000;


    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters.TryGetValue("Alpha", out var a))
            _alpha = Convert.ToDouble(a);

        if (parameters.TryGetValue("MaxIterations", out var it))
            _iterations = Convert.ToInt32(it);
    }

    public void Fit(Matrix X, VectorN y)
    {
        var Xb = X.WithBiasColumn();
        int n = Xb.rowLength;
        int p = Xb.columnLength;

        _weights = new VectorN(new double[p]);

        for (int iter = 0; iter < _iterations; iter++)
        {
            for (int j = 0; j < p; j++)
            {
                double num = 0;
                double den = 0;

                for (int i = 0; i < n; i++)
                {
                    double pred = 0;
                    for (int k = 0; k < p; k++)
                        if (k != j)
                            pred += Xb.values[i, k] * _weights[k];

                    num += Xb.values[i, j] * (y[i] - pred);
                    den += Xb.values[i, j] * Xb.values[i, j];
                }

                _weights[j] = SoftThreshold(num / den, _alpha);
            }
        }
    }

    public VectorN Predict(Matrix X)
    {
        var Xb = X.WithBiasColumn();
        return Xb * _weights;
    }

    private double SoftThreshold(double z, double a)
    {
        if (z > a) return z - a;
        if (z < -a) return z + a;
        return 0;
    }

    public IModel Clone()
    {
        var clone = new Lasso();
        clone.SetHyperParameters(new Dictionary<string, object>
        {
            ["Alpha"] = _alpha,
            ["MaxIterations"] = _iterations
        });
        return clone;
    }

}
