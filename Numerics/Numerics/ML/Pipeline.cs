using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.ML.Selector.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CSharpNumerics.ML;

public class Pipeline
{
    public IModel Model { get; }
    public IScaler Scaler { get; }
    public ISelector Selector { get; }

    public Dictionary<string, object> ModelParams { get; }
    public Dictionary<string, object> ScalerParams { get; }
    public Dictionary<string, object> SelectorParams { get; }

    private bool _isFitted = false;

    public Pipeline(
        IModel model,
        Dictionary<string, object> modelParams,
        IScaler scaler = null,
         Dictionary<string, object> scalerParams = null,
        ISelector selector = null,
         Dictionary<string, object> selectorParams = null)
    {
        Model = model;
        ModelParams = modelParams;

        Scaler = scaler;
        ScalerParams = scalerParams ?? [];

        Selector = selector;
        SelectorParams = selectorParams ?? [];

        if (model is IHasHyperparameters hpModel)
            hpModel.SetHyperParameters(modelParams);

        if (selector is IHasHyperparameters hpSelector)
            hpSelector.SetHyperParameters(selectorParams);
    }




    public void Fit(Matrix X, VectorN y)
    {
        // Apply selector
        if (Selector != null)
            X = Selector.FitTransform(X, y);

        // Apply scaler
        if (Scaler != null)
            X = Scaler.FitTransform(X);

        // Apply model hyperparameters
        ApplyHyperparameters(Model, ModelParams);

        Model.Fit(X, y);
        _isFitted = true;
    }

    public VectorN Predict(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("Pipeline has not been fitted.");

        if (Selector != null)
            X = Selector.Transform(X);

        if (Scaler != null)
            X = Scaler.Transform(X);

        return Model.Predict(X);
    }
    public Pipeline Clone()
    {
        return new Pipeline(
            Model.Clone(),
            ModelParams,
            Scaler?.Clone(),
            ScalerParams,
            Selector?.Clone(),
            SelectorParams
        );
    }

    private void ApplyHyperparameters(object target, Dictionary<string, object> parameters)
    {
        if (parameters != null)
        {
            var type = target.GetType();
            foreach (var param in parameters)
            {
                var prop = type.GetProperty(param.Key);
                prop?.SetValue(target, param.Value);
            }
        }
    }

}