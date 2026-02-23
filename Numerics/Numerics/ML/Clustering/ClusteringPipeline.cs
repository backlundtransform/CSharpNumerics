using CSharpNumerics.ML.Clustering.Interfaces;
using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.ML.Scalers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.ML.Clustering;

/// <summary>
/// Unsupervised clustering pipeline: optional scaler + clustering model.
/// Mirrors the supervised <see cref="Pipeline"/> but without selectors or target vectors.
/// </summary>
public class ClusteringPipeline
{
    public IClusteringModel Model { get; }
    public IScaler Scaler { get; }

    public Dictionary<string, object> ModelParams { get; }
    public Dictionary<string, object> ScalerParams { get; }

    private bool _isFitted;

    public ClusteringPipeline(
        IClusteringModel model,
        Dictionary<string, object> modelParams,
        IScaler scaler = null,
        Dictionary<string, object> scalerParams = null)
    {
        Model = model;
        ModelParams = modelParams ?? new Dictionary<string, object>();

        Scaler = scaler;
        ScalerParams = scalerParams ?? new Dictionary<string, object>();

        if (model is IHasHyperparameters hp)
            hp.SetHyperParameters(ModelParams);
    }

    // ── Fit ──────────────────────────────────────────────────────
    public void Fit(Matrix X)
    {
        if (Scaler != null)
            X = Scaler.FitTransform(X);

        ApplyHyperparameters(Model, ModelParams);
        Model.Fit(X);
        _isFitted = true;
    }

    // ── Predict ──────────────────────────────────────────────────
    public VectorN Predict(Matrix X)
    {
        if (!_isFitted)
            throw new InvalidOperationException("ClusteringPipeline has not been fitted.");

        if (Scaler != null)
            X = Scaler.Transform(X);

        return Model.Predict(X);
    }

    // ── FitPredict ───────────────────────────────────────────────
    public VectorN FitPredict(Matrix X)
    {
        if (Scaler != null)
            X = Scaler.FitTransform(X);

        ApplyHyperparameters(Model, ModelParams);
        var labels = Model.FitPredict(X);
        _isFitted = true;
        return labels;
    }

    // ── Clone ────────────────────────────────────────────────────
    public ClusteringPipeline Clone()
    {
        return new ClusteringPipeline(
            Model.Clone(),
            new Dictionary<string, object>(ModelParams),
            Scaler?.Clone(),
            ScalerParams != null ? new Dictionary<string, object>(ScalerParams) : null
        );
    }

    public override string ToString()
    {
        var name = Model.GetType().Name;
        if (Scaler != null)
            name = $"{Scaler.GetType().Name} → {name}";

        if (ModelParams.Count > 0)
        {
            var parts = new List<string>();
            foreach (var kvp in ModelParams)
                parts.Add($"{kvp.Key}={kvp.Value}");
            name += $" ({string.Join(", ", parts)})";
        }

        return name;
    }

    private void ApplyHyperparameters(object target, Dictionary<string, object> parameters)
    {
        if (parameters == null) return;

        var type = target.GetType();
        foreach (var param in parameters)
        {
            var prop = type.GetProperty(param.Key);
            prop?.SetValue(target, param.Value);
        }
    }
}
