using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace CSharpNumerics.Engines.Exoplanet.Features;

/// <summary>
/// Creates windowed training data for ML classification by combining raw phase-folded
/// light curve windows with extracted transit features as additional columns.
/// Compatible with <see cref="CSharpNumerics.ML.Sequence.SequenceDataHelper.CreateWindows"/>.
/// </summary>
public static class WindowedFeatureExtractor
{
    /// <summary>
    /// Creates ML training data from a set of light curves and their labels.
    /// Each light curve is phase-folded, windowed, and augmented with extracted features.
    /// </summary>
    /// <param name="curves">Input light curves.</param>
    /// <param name="candidates">Transit candidates (one per light curve) with extracted features.</param>
    /// <param name="labels">Binary labels: 1.0 = transit, 0.0 = no transit.</param>
    /// <param name="windowSize">Number of phase bins per sample window.</param>
    /// <param name="phaseBins">Number of bins for phase folding (default: 200).</param>
    /// <returns>
    /// (X, y) where X has shape [numSamples × (windowSize + numFeatures)] and y has length numSamples.
    /// Each row contains the phase-folded flux window followed by the extracted feature values.
    /// </returns>
    public static (Matrix X, VectorN y) CreateTrainingData(
        LightCurve[] curves, TransitCandidate[] candidates, double[] labels,
        int windowSize, int phaseBins = 200)
    {
        if (curves == null) throw new ArgumentNullException(nameof(curves));
        if (candidates == null) throw new ArgumentNullException(nameof(candidates));
        if (labels == null) throw new ArgumentNullException(nameof(labels));
        if (curves.Length != candidates.Length || curves.Length != labels.Length)
            throw new ArgumentException("Curves, candidates, and labels must have the same length.");
        if (windowSize <= 0) throw new ArgumentOutOfRangeException(nameof(windowSize));

        int numFeatures = TransitFeatureSet.FeatureNames.All.Length;
        int rowWidth = windowSize + numFeatures;
        int numSamples = curves.Length;

        var X = new Matrix(numSamples, rowWidth);
        var y = new VectorN(numSamples);

        for (int s = 0; s < numSamples; s++)
        {
            y[s] = labels[s];

            // Phase-fold and bin the light curve
            double[] phaseBinFlux = PhaseFoldToBins(
                curves[s], candidates[s].Parameters, windowSize, phaseBins);

            // Fill flux window columns
            for (int j = 0; j < windowSize; j++)
                X.values[s, j] = phaseBinFlux[j];

            // Append extracted features
            var features = candidates[s].Features;
            string[] featureNames = TransitFeatureSet.FeatureNames.All;
            for (int f = 0; f < featureNames.Length; f++)
            {
                double val = features != null && features.HasFeature(featureNames[f])
                    ? features[featureNames[f]]
                    : 0.0;
                X.values[s, windowSize + f] = val;
            }
        }

        return (X, y);
    }

    /// <summary>
    /// Creates windowed data from a single light curve for inference (no labels).
    /// </summary>
    public static Matrix CreateInferenceData(
        LightCurve curve, TransitCandidate candidate, int windowSize, int phaseBins = 200)
    {
        if (curve == null) throw new ArgumentNullException(nameof(curve));
        if (candidate == null) throw new ArgumentNullException(nameof(candidate));
        if (windowSize <= 0) throw new ArgumentOutOfRangeException(nameof(windowSize));

        int numFeatures = TransitFeatureSet.FeatureNames.All.Length;
        int rowWidth = windowSize + numFeatures;
        var X = new Matrix(1, rowWidth);

        double[] phaseBinFlux = PhaseFoldToBins(curve, candidate.Parameters, windowSize, phaseBins);

        for (int j = 0; j < windowSize; j++)
            X.values[0, j] = phaseBinFlux[j];

        var features = candidate.Features;
        string[] featureNames = TransitFeatureSet.FeatureNames.All;
        for (int f = 0; f < featureNames.Length; f++)
        {
            double val = features != null && features.HasFeature(featureNames[f])
                ? features[featureNames[f]]
                : 0.0;
            X.values[0, windowSize + f] = val;
        }

        return X;
    }

    /// <summary>
    /// Phase-folds a light curve and bins into a fixed-size array centered on the transit.
    /// </summary>
    private static double[] PhaseFoldToBins(LightCurve lc, TransitParameters p, int windowSize, int phaseBins)
    {
        if (lc.Length < 3 || p.Period <= 0)
            return new double[windowSize];

        var times = new VectorN(lc.Time);
        var values = new VectorN(lc.Flux);
        var (binCenters, binMeans) = PhaseFolding.FoldAndBin(
            times, values, p.Period, phaseBins, p.Epoch);

        // Center phases around 0 and sort
        double[] centeredPhases = new double[binCenters.Length];
        double[] sortedMeans = new double[binMeans.Length];
        int[] order = new int[binCenters.Length];

        for (int i = 0; i < binCenters.Length; i++)
        {
            centeredPhases[i] = binCenters[i] > 0.5 ? binCenters[i] - 1.0 : binCenters[i];
            order[i] = i;
        }

        Array.Sort((double[])centeredPhases.Clone(), order);
        double[] sortedPhases = new double[centeredPhases.Length];
        for (int i = 0; i < order.Length; i++)
        {
            sortedPhases[i] = centeredPhases[order[i]];
            sortedMeans[i] = binMeans[order[i]];
        }

        // Extract a centered window of windowSize bins
        int centerIdx = sortedPhases.Length / 2;
        int halfWindow = windowSize / 2;
        int startIdx = centerIdx - halfWindow;

        double[] result = new double[windowSize];
        for (int i = 0; i < windowSize; i++)
        {
            int idx = startIdx + i;
            if (idx >= 0 && idx < sortedMeans.Length)
                result[i] = sortedMeans[idx];
            else
                result[i] = 1.0; // baseline
        }

        return result;
    }
}
