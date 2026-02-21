
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
namespace CSharpNumerics.ML.Models.Regression;

public class Ridge : IRegressionModel, IHasHyperparameters
{
    public double Alpha { get; set; } = 1.0; 
    public bool FitIntercept { get; set; } = true;

    private VectorN _weights;

    public bool SupportsScaling => true;
    public bool SupportsFeatureSelection => true;
    public bool IsProbabilistic => false;

    public void Fit(Matrix X, VectorN y)
    {

        var Xb = FitIntercept ? X.WithBiasColumn() : X;

        var Xt = Xb.Transpose();
        var XtX = Xt * Xb;


        var reg = new Matrix(XtX.rowLength, XtX.columnLength);
        int start = FitIntercept ? 1 : 0;
        for (int i = start; i < reg.rowLength; i++)
            reg.values[i, i] = Alpha; 


        var inv = (XtX + reg).Inverse();
        _weights = inv * (Xt * y);
    }

    public VectorN Predict(Matrix X)
    {
        var Xb = FitIntercept ? X.WithBiasColumn() : X;
        return Xb * _weights;
    }

    public void SetHyperParameters(Dictionary<string, object> parameters)
    {
        if (parameters.TryGetValue("Alpha", out var alpha))
            Alpha = Convert.ToDouble(alpha);
        if (parameters.TryGetValue("FitIntercept", out object value))
            FitIntercept = Convert.ToBoolean(value);
    }

    public IModel Clone()
    {
        return new Ridge
        {
            Alpha = Alpha,
            FitIntercept = FitIntercept
        };
    }
}
