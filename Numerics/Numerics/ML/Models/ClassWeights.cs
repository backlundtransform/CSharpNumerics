using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.ML.Models;

/// <summary>Sample weights that balance the classes of a classification target.</summary>
public static class ClassWeights
{
    /// <summary>
    /// A weight per sample of n / (classes × n_c), where n_c is the size of the sample's class:
    /// every class then carries the same total weight, and the weights average 1.
    /// </summary>
    public static VectorN Balanced(VectorN y)
    {
        if (y.Values == null || y.Length == 0) throw new ArgumentException("y is empty.", nameof(y));

        var counts = new Dictionary<double, int>();
        foreach (double label in y.Values)
            counts[label] = counts.TryGetValue(label, out int c) ? c + 1 : 1;

        double n = y.Length;
        int classes = counts.Count;
        return new VectorN(y.Values.Select(label => n / (classes * counts[label])).ToArray());
    }

    /// <summary>
    /// Validates sample weights for <paramref name="rows"/> rows and returns them as an array:
    /// non-negative, finite, not all zero. No weights (null or an empty vector) gives all ones.
    /// </summary>
    internal static double[] Validate(VectorN? sampleWeights, int rows)
    {
        if (sampleWeights is not { Values: not null } given)
            return Enumerable.Repeat(1.0, rows).ToArray();
        if (given.Length != rows)
            throw new ArgumentException($"Expected {rows} sample weights, got {given.Length}.", nameof(sampleWeights));

        double total = 0;
        foreach (double w in given.Values)
        {
            if (!(w >= 0) || double.IsInfinity(w))
                throw new ArgumentException("Sample weights must be finite and non-negative.", nameof(sampleWeights));
            total += w;
        }
        if (total <= 0)
            throw new ArgumentException("At least one sample weight must be positive.", nameof(sampleWeights));

        return (double[])given.Values.Clone();
    }
}
