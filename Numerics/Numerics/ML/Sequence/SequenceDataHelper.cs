using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;

namespace CSharpNumerics.ML.Sequence;

/// <summary>
/// Utilities for converting time series data into windowed samples suitable for
/// sequence models (CNN1D, LSTM, BiLSTM). Each window becomes a single row in
/// the output Matrix with shape (windowSize × features) flattened.
/// </summary>
public static class SequenceDataHelper
{
    /// <summary>
    /// Creates sliding-window samples from a <see cref="TimeSeries"/>.
    /// Each window of <paramref name="windowSize"/> consecutive timesteps becomes one flattened row.
    /// </summary>
    /// <param name="timeSeries">Source time series with one or more feature columns.</param>
    /// <param name="windowSize">Number of consecutive timesteps per sample.</param>
    /// <param name="labelColumnIndex">
    /// Column index in <see cref="TimeSeries.Data"/> to use as the label.
    /// The label is taken from the <em>last</em> timestep in each window.
    /// This column is excluded from the feature matrix.
    /// </param>
    /// <param name="stride">Step size between consecutive windows (default: 1).</param>
    /// <returns>
    /// A tuple (X, y) where X has shape [numWindows × (windowSize × featureCount)]
    /// and y has length numWindows.
    /// </returns>
    public static (Matrix X, VectorN y) CreateWindows(TimeSeries timeSeries, int windowSize, int labelColumnIndex, int stride = 1)
    {
        if (timeSeries is null) throw new ArgumentNullException(nameof(timeSeries));
        if (windowSize <= 0) throw new ArgumentOutOfRangeException(nameof(windowSize), "Window size must be positive.");
        if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride), "Stride must be positive.");
        if (labelColumnIndex < 0 || labelColumnIndex >= timeSeries.ColumnCount)
            throw new ArgumentOutOfRangeException(nameof(labelColumnIndex), "Label column index is out of range.");
        if (windowSize > timeSeries.RowCount)
            throw new ArgumentException("Window size cannot exceed the number of rows in the time series.");

        int totalRows = timeSeries.RowCount;
        int featureCount = timeSeries.ColumnCount - 1; // Exclude label column

        if (featureCount <= 0)
            throw new ArgumentException("Time series must have at least one feature column besides the label column.");

        int numWindows = ((totalRows - windowSize) / stride) + 1;
        int flattenedWidth = windowSize * featureCount;

        var X = new Matrix(numWindows, flattenedWidth);
        var y = new VectorN(numWindows);

        // Build feature column index map (skip labelColumnIndex)
        int[] featureCols = new int[featureCount];
        int fi = 0;
        for (int c = 0; c < timeSeries.ColumnCount; c++)
        {
            if (c != labelColumnIndex)
                featureCols[fi++] = c;
        }

        for (int w = 0; w < numWindows; w++)
        {
            int windowStart = w * stride;

            // Label: value at the last timestep of this window
            y[w] = timeSeries.Data[labelColumnIndex][windowStart + windowSize - 1];

            // Features: flatten (timestep, feature) in row-major order
            int col = 0;
            for (int t = 0; t < windowSize; t++)
            {
                int rowIndex = windowStart + t;
                for (int f = 0; f < featureCount; f++)
                {
                    X.values[w, col++] = timeSeries.Data[featureCols[f]][rowIndex];
                }
            }
        }

        return (X, y);
    }

    /// <summary>
    /// Creates sliding-window samples from raw arrays (no label column).
    /// Useful when labels are computed separately (e.g., binary transit/no-transit).
    /// </summary>
    /// <param name="featureColumns">
    /// Array of feature columns, each of length numTimesteps. Column-major like <see cref="TimeSeries.Data"/>.
    /// </param>
    /// <param name="labels">
    /// One label per timestep. The label for each window is taken from the <em>last</em> timestep.
    /// </param>
    /// <param name="windowSize">Number of consecutive timesteps per sample.</param>
    /// <param name="stride">Step size between consecutive windows (default: 1).</param>
    public static (Matrix X, VectorN y) CreateWindows(double[][] featureColumns, double[] labels, int windowSize, int stride = 1)
    {
        if (featureColumns is null) throw new ArgumentNullException(nameof(featureColumns));
        if (labels is null) throw new ArgumentNullException(nameof(labels));
        if (featureColumns.Length == 0) throw new ArgumentException("At least one feature column is required.", nameof(featureColumns));
        if (windowSize <= 0) throw new ArgumentOutOfRangeException(nameof(windowSize));
        if (stride <= 0) throw new ArgumentOutOfRangeException(nameof(stride));

        int totalRows = featureColumns[0].Length;
        if (labels.Length != totalRows)
            throw new ArgumentException("Labels length must match the number of timesteps.");

        int featureCount = featureColumns.Length;
        int numWindows = ((totalRows - windowSize) / stride) + 1;
        int flattenedWidth = windowSize * featureCount;

        var X = new Matrix(numWindows, flattenedWidth);
        var y = new VectorN(numWindows);

        for (int w = 0; w < numWindows; w++)
        {
            int windowStart = w * stride;
            y[w] = labels[windowStart + windowSize - 1];

            int col = 0;
            for (int t = 0; t < windowSize; t++)
            {
                int rowIndex = windowStart + t;
                for (int f = 0; f < featureCount; f++)
                {
                    X.values[w, col++] = featureColumns[f][rowIndex];
                }
            }
        }

        return (X, y);
    }
}
