using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Numerics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Multiphysics.MonteCarlo;

/// <summary>
/// Result of a batch Monte Carlo multiphysics run.
/// Contains the scenario matrix (rows = iterations, columns = flattened output)
/// and per-iteration results for further analysis.
/// </summary>
public class MultiphysicsScenarioResult
{
    /// <summary>
    /// Scenario matrix: each row is one iteration, each column is one output value.
    /// For 2D simulations: columns = Nx × Ny (final field, flattened).
    /// For 1D simulations: columns = NodeCount (final values).
    /// </summary>
    public Matrix ScenarioMatrix { get; }

    /// <summary>Full result per iteration.</summary>
    public List<SimulationResult> Results { get; }

    /// <summary>Simulation type.</summary>
    public MultiphysicsType Type { get; }

    /// <summary>Number of Monte Carlo iterations.</summary>
    public int Iterations => Results.Count;

    /// <summary>Number of features (columns) per iteration.</summary>
    public int FeatureCount => ScenarioMatrix.columnLength;

    internal MultiphysicsScenarioResult(
        Matrix scenarioMatrix,
        List<SimulationResult> results,
        MultiphysicsType type)
    {
        ScenarioMatrix = scenarioMatrix;
        Results = results;
        Type = type;
    }

    /// <summary>
    /// Returns the distribution of a single feature (column) across all iterations.
    /// </summary>
    /// <param name="featureIndex">Column index in the scenario matrix.</param>
    /// <returns>Array of length <see cref="Iterations"/>.</returns>
    public double[] GetFeatureDistribution(int featureIndex)
    {
        var result = new double[Iterations];
        for (int i = 0; i < Iterations; i++)
            result[i] = ScenarioMatrix.values[i, featureIndex];
        return result;
    }

    /// <summary>
    /// Returns the flattened output vector for a single iteration.
    /// </summary>
    public double[] GetScenarioVector(int iterationIndex)
    {
        int cols = ScenarioMatrix.columnLength;
        var row = new double[cols];
        for (int j = 0; j < cols; j++)
            row[j] = ScenarioMatrix.values[iterationIndex, j];
        return row;
    }

    /// <summary>
    /// Computes percentile maps across iterations.
    /// For 2D: returns [ix, iy] array where each cell is the percentile value.
    /// For 1D: returns a 1D array where each node is the percentile value.
    /// </summary>
    /// <param name="percentile">Percentile in [0, 100] (e.g. 95).</param>
    public double[] ComputePercentile(double percentile)
    {
        int n = Iterations;
        int cols = FeatureCount;
        var result = new double[cols];
        var column = new double[n];

        for (int j = 0; j < cols; j++)
        {
            for (int i = 0; i < n; i++)
                column[i] = ScenarioMatrix.values[i, j];
            System.Array.Sort(column);

            double rank = percentile / 100.0 * (n - 1);
            int lo = (int)rank;
            double frac = rank - lo;
            if (lo >= n - 1)
                result[j] = column[n - 1];
            else
                result[j] = column[lo] * (1.0 - frac) + column[lo + 1] * frac;
        }

        return result;
    }
}
