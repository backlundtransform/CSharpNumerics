using CSharpNumerics.ML.DimensionalityReduction.Interfaces;
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
    public IDimensionalityReducer Reducer { get; }

    public Dictionary<string, object> ModelParams { get; }
    public Dictionary<string, object> ScalerParams { get; }
    public Dictionary<string, object> SelectorParams { get; }
    public Dictionary<string, object> ReducerParams { get; }

    private bool _isFitted = false;

    public Pipeline(
        IModel model,
        Dictionary<string, object> modelParams,
        IScaler scaler = null,
         Dictionary<string, object> scalerParams = null,
        ISelector selector = null,
         Dictionary<string, object> selectorParams = null,
        IDimensionalityReducer reducer = null,
        Dictionary<string, object> reducerParams = null)
    {
        Model = model;
        ModelParams = modelParams;

        Scaler = scaler;
        ScalerParams = scalerParams ?? [];

        Selector = selector;
        SelectorParams = selectorParams ?? [];

        Reducer = reducer;
        ReducerParams = reducerParams ?? [];

        if (model is IHasHyperparameters hpModel)
            hpModel.SetHyperParameters(modelParams);

        if (selector is IHasHyperparameters hpSelector)
            hpSelector.SetHyperParameters(selectorParams);

        if (reducer is IHasHyperparameters hpReducer)
            hpReducer.SetHyperParameters(reducerParams);
    }




    public void Fit(Matrix X, VectorN y) => Fit(X, y, (VectorN?)null);

    /// <summary>
    /// Fits the pipeline, training the model with a weight per row. The selector, reducer and
    /// scaler are fitted without weights. Null weights is the same as <see cref="Fit(Matrix, VectorN)"/>;
    /// otherwise the model must implement <see cref="ISampleWeightedModel"/>.
    /// </summary>
    public void Fit(Matrix X, VectorN y, VectorN? sampleWeights)
    {
        if (sampleWeights.HasValue && Model is not ISampleWeightedModel)
            throw new InvalidOperationException(
                $"{Model.GetType().Name} cannot be trained with sample weights (it does not implement ISampleWeightedModel).");

        // Apply selector
        if (Selector != null)
            X = Selector.FitTransform(X, y);

        // Apply reducer
        if (Reducer != null)
            X = Reducer.FitTransform(X);

        // Apply scaler
        if (Scaler != null)
            X = Scaler.FitTransform(X);

        // Apply model hyperparameters
        ApplyHyperparameters(Model, ModelParams);

        if (sampleWeights.HasValue)
            ((ISampleWeightedModel)Model).Fit(X, y, sampleWeights.Value);
        else
            Model.Fit(X, y);
        _isFitted = true;
    }

    public VectorN Predict(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("Pipeline has not been fitted.");

        if (Selector != null)
            X = Selector.Transform(X);

        if (Reducer != null)
            X = Reducer.Transform(X);

        if (Scaler != null)
            X = Scaler.Transform(X);

        return Model.Predict(X);
    }

    /// <summary>
    /// Class probabilities from the model (one column per class), after the same selector, reducer
    /// and scaler as <see cref="Predict"/>. The model must implement <see cref="IProbabilisticClassifier"/>.
    /// </summary>
    public Matrix PredictProbabilities(Matrix X)
    {
        if (Model is not IProbabilisticClassifier probabilistic)
            throw new InvalidOperationException(
                $"{Model.GetType().Name} does not give probabilities (it does not implement IProbabilisticClassifier).");

        return probabilistic.PredictProbabilities(Transform(X));
    }

    private Matrix Transform(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("Pipeline has not been fitted.");

        if (Selector != null)
            X = Selector.Transform(X);

        if (Reducer != null)
            X = Reducer.Transform(X);

        if (Scaler != null)
            X = Scaler.Transform(X);

        return X;
    }
    public Pipeline Clone()
    {
        return new Pipeline(
            Model.Clone(),
            ModelParams,
            Scaler?.Clone(),
            ScalerParams,
            Selector?.Clone(),
            SelectorParams,
            Reducer?.Clone(),
            ReducerParams
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