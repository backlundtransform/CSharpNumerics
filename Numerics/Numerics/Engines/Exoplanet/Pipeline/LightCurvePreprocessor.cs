using System;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public static class LightCurvePreprocessor
{
    public static LightCurve Preprocess(LightCurve lc, TransitDetectionConfig config)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));
        if (config == null) throw new ArgumentNullException(nameof(config));

        var result = LightCurveSanitizer.RemoveBadQuality(lc);
        result = Detrend(result, config);
        result = RemoveUpwardOutliers(result, config.OutlierSigmaThreshold);
        result = LightCurveSanitizer.NormalizeFlux(result);

        return result;
    }

    public static LightCurve Detrend(LightCurve lc, TransitDetectionConfig config)
    {
        var times = new VectorN(lc.Time);
        var values = new VectorN(lc.Flux);

        VectorN detrended;
        VectorN trend;

        switch (config.DetrendingMethod)
        {
            case DetrendingMethod.Polynomial:
                (detrended, trend) = TimeSeriesDetrending.PolynomialDetrend(times, values, config.DetrendingPolyDegree);
                break;
            case DetrendingMethod.MovingAverage:
                (detrended, trend) = TimeSeriesDetrending.MovingAverageDetrend(values, config.DetrendingWindowSize);
                break;
            case DetrendingMethod.SavitzkyGolay:
                (detrended, trend) = TimeSeriesDetrending.SavitzkyGolayDetrend(values, config.DetrendingWindowSize, config.DetrendingPolyDegree);
                break;
            case DetrendingMethod.MedianFilter:
            default:
                (detrended, trend) = TimeSeriesDetrending.MedianFilterDetrend(values, config.DetrendingWindowSize);
                break;
        }

        // For transit analysis, divide by trend rather than subtract,
        // to preserve relative flux and transit depths.
        double[] correctedFlux = new double[lc.Length];
        double[] correctedErr = new double[lc.Length];
        for (int i = 0; i < lc.Length; i++)
        {
            double t = trend[i];
            if (Math.Abs(t) > 1e-10)
            {
                correctedFlux[i] = values[i] / t;
                correctedErr[i] = lc.FluxError[i] / Math.Abs(t);
            }
            else
            {
                correctedFlux[i] = values[i];
                correctedErr[i] = lc.FluxError[i];
            }
        }

        return new LightCurve(
            (double[])lc.Time.Clone(), correctedFlux, correctedErr,
            (int[])lc.QualityFlags.Clone(), lc.Metadata);
    }

    /// <summary>
    /// Removes only upward outliers (flares, cosmic rays) while preserving
    /// downward outliers (transits) that should not be clipped.
    /// </summary>
    private static LightCurve RemoveUpwardOutliers(LightCurve lc, double sigmaThreshold)
    {
        int count = lc.Length;
        if (count < 5) return lc;

        double median = MedianOf(lc.Flux);
        double mad = MedianAbsDeviation(lc.Flux, median);
        double sigma = mad * 1.4826; // scale MAD to approximate std dev
        if (sigma < 1e-15) return lc;

        double upperLimit = median + sigmaThreshold * sigma;

        var indices = new System.Collections.Generic.List<int>();
        for (int i = 0; i < count; i++)
        {
            if (lc.Flux[i] <= upperLimit)
                indices.Add(i);
        }

        if (indices.Count == count) return lc;

        int n = indices.Count;
        double[] time = new double[n];
        double[] flux = new double[n];
        double[] err = new double[n];
        int[] qual = new int[n];

        for (int i = 0; i < n; i++)
        {
            int idx = indices[i];
            time[i] = lc.Time[idx];
            flux[i] = lc.Flux[idx];
            err[i] = lc.FluxError[idx];
            qual[i] = lc.QualityFlags[idx];
        }

        return new LightCurve(time, flux, err, qual, lc.Metadata);
    }

    private static double MedianOf(double[] arr)
    {
        double[] sorted = (double[])arr.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 1 ? sorted[n / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }

    private static double MedianAbsDeviation(double[] arr, double median)
    {
        double[] absDevs = new double[arr.Length];
        for (int i = 0; i < arr.Length; i++)
            absDevs[i] = Math.Abs(arr[i] - median);
        Array.Sort(absDevs);
        int n = absDevs.Length;
        return n % 2 == 1 ? absDevs[n / 2] : (absDevs[n / 2 - 1] + absDevs[n / 2]) / 2.0;
    }
}
