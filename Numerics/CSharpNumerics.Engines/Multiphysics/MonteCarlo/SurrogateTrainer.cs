using CSharpNumerics.ML.Models.Interfaces;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.MonteCarlo;

/// <summary>
/// Trains a regression surrogate on Monte Carlo scenario data.
/// The surrogate can predict simulation output without re-running the solver,
/// enabling fast what-if analysis and sensitivity studies.
/// <para>
/// Input features are sampled parameter values (material props, BCs, loads).
/// Targets are output values (per-cell temperature, deflection, etc.) or a
/// scalar summary (max value).
/// </para>
/// </summary>
public class SurrogateTrainer
{
    private readonly MultiphysicsScenarioResult _scenarios;
    private readonly IModel _model;
    private IModel _trained;

    /// <summary>Whether the surrogate has been trained.</summary>
    public bool IsTrained => _trained != null;

    /// <summary>
    /// Creates a surrogate trainer for the given MC results.
    /// </summary>
    /// <param name="scenarios">Monte Carlo scenario results (provides the training data).</param>
    /// <param name="model">Regression model to train (e.g. <c>new Linear()</c>, <c>new MLPRegressor()</c>).</param>
    public SurrogateTrainer(MultiphysicsScenarioResult scenarios, IModel model)
    {
        _scenarios = scenarios ?? throw new ArgumentNullException(nameof(scenarios));
        _model = model ?? throw new ArgumentNullException(nameof(model));
    }

    /// <summary>
    /// Trains the surrogate on the scenario matrix.
    /// Uses the scenario matrix as features X and a target extraction function y.
    /// </summary>
    /// <param name="targetExtractor">
    /// Function that extracts a scalar target from each <see cref="SimulationResult"/>.
    /// Defaults to <c>r => r.MaxValue</c>.
    /// </param>
    public void Train(Func<SimulationResult, double> targetExtractor = null)
    {
        targetExtractor ??= r => r.MaxValue;

        int n = _scenarios.Iterations;
        int cols = _scenarios.FeatureCount;

        // X = scenario matrix (each row is the flattened output of one trial)
        // y = scalar target per trial
        var targets = new double[n];
        for (int i = 0; i < n; i++)
            targets[i] = targetExtractor(_scenarios.Results[i]);

        _trained = _model.Clone();
        _trained.Fit(_scenarios.ScenarioMatrix, new VectorN(targets));
    }

    /// <summary>
    /// Trains the surrogate with explicit feature/target matrices.
    /// Use when you want to provide custom inputs (e.g. sampled parameters)
    /// rather than the full scenario matrix.
    /// </summary>
    /// <param name="X">Feature matrix (rows = iterations, cols = input parameters).</param>
    /// <param name="y">Target vector (one scalar per iteration).</param>
    public void Train(Matrix X, VectorN y)
    {
        _trained = _model.Clone();
        _trained.Fit(X, y);
    }

    /// <summary>
    /// Predict targets for new input scenarios.
    /// </summary>
    /// <param name="X">Input matrix (same column schema as training data).</param>
    /// <returns>Predicted target vector.</returns>
    public VectorN Predict(Matrix X)
    {
        if (_trained == null)
            throw new InvalidOperationException("Surrogate not trained. Call Train() first.");
        return _trained.Predict(X);
    }

    /// <summary>
    /// Predict the target for a single scenario vector.
    /// </summary>
    /// <param name="scenario">Feature vector (same length as training columns).</param>
    /// <returns>Predicted scalar target.</returns>
    public double PredictOne(double[] scenario)
    {
        if (_trained == null)
            throw new InvalidOperationException("Surrogate not trained. Call Train() first.");

        var matrix = new double[1, scenario.Length];
        for (int j = 0; j < scenario.Length; j++)
            matrix[0, j] = scenario[j];

        var result = _trained.Predict(new Matrix(matrix));
        return result[0];
    }
}
