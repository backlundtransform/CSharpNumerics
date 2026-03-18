using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.SingleObjective;
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

        double l1 = Lambda * L1Ratio;
        double l2 = Lambda * (1 - L1Ratio);

        var cd = new CoordinateDescent(maxIterations: MaxIterations);
        var result = cd.Solve(Xb.values, y.Values, l1, l2, skipBiasRegularisation: true);

        _weights = new VectorN(result);
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
