using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Exoplanet.Data;

public static class LightCurveSanitizer
{
    public static LightCurve RemoveBadQuality(LightCurve lc, int qualityMask = 0)
    {
        var indices = new List<int>();
        for (int i = 0; i < lc.Length; i++)
        {
            if ((lc.QualityFlags[i] & ~qualityMask) == 0)
                indices.Add(i);
        }

        return SelectIndices(lc, indices);
    }

    public static LightCurve RemoveOutliers(LightCurve lc, double sigmaThreshold = 5.0)
    {
        double mean = 0;
        int count = lc.Length;
        for (int i = 0; i < count; i++)
            mean += lc.Flux[i];
        mean /= count;

        double variance = 0;
        for (int i = 0; i < count; i++)
        {
            double diff = lc.Flux[i] - mean;
            variance += diff * diff;
        }
        double std = Math.Sqrt(variance / count);

        if (std < 1e-15)
            return lc;

        var indices = new List<int>();
        for (int i = 0; i < count; i++)
        {
            double z = Math.Abs(lc.Flux[i] - mean) / std;
            if (z <= sigmaThreshold)
                indices.Add(i);
        }

        return SelectIndices(lc, indices);
    }

    public static LightCurve FillGaps(LightCurve lc, double maxGapSize)
    {
        if (lc.Length < 2) return lc;

        var time = new List<double>(lc.Time);
        var flux = new List<double>(lc.Flux);
        var err = new List<double>(lc.FluxError);
        var qual = new List<int>(lc.QualityFlags);

        double medianCadence = MedianCadence(lc.Time);
        if (medianCadence < 1e-15) return lc;

        int i = 0;
        while (i < time.Count - 1)
        {
            double gap = time[i + 1] - time[i];
            if (gap > medianCadence * 1.5 && gap <= maxGapSize)
            {
                int nFill = (int)Math.Round(gap / medianCadence) - 1;
                double dt = gap / (nFill + 1);

                double fluxLeft = flux[i];
                double fluxRight = flux[i + 1];

                for (int j = 1; j <= nFill; j++)
                {
                    double t = time[i] + j * dt;
                    double f = fluxLeft + (fluxRight - fluxLeft) * j / (nFill + 1);
                    double e = (err[i] + err[i + 1]) / 2.0;

                    time.Insert(i + j, t);
                    flux.Insert(i + j, f);
                    err.Insert(i + j, e);
                    qual.Insert(i + j, 0);
                }

                i += nFill + 1;
            }
            else
            {
                i++;
            }
        }

        return new LightCurve(
            time.ToArray(), flux.ToArray(), err.ToArray(), qual.ToArray(), lc.Metadata);
    }

    public static LightCurve NormalizeFlux(LightCurve lc)
    {
        if (lc.Length == 0) return lc;

        double median = Median(lc.Flux);
        if (Math.Abs(median) < 1e-15)
            return lc;

        double[] normalizedFlux = new double[lc.Length];
        double[] normalizedErr = new double[lc.Length];

        for (int i = 0; i < lc.Length; i++)
        {
            normalizedFlux[i] = lc.Flux[i] / median;
            normalizedErr[i] = lc.FluxError[i] / median;
        }

        return new LightCurve(
            (double[])lc.Time.Clone(), normalizedFlux, normalizedErr,
            (int[])lc.QualityFlags.Clone(), lc.Metadata);
    }

    private static LightCurve SelectIndices(LightCurve lc, List<int> indices)
    {
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

    private static double Median(double[] values)
    {
        double[] sorted = (double[])values.Clone();
        Array.Sort(sorted);
        int n = sorted.Length;
        if (n % 2 == 1)
            return sorted[n / 2];
        return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }

    private static double MedianCadence(double[] time)
    {
        if (time.Length < 2) return 0;
        double[] diffs = new double[time.Length - 1];
        for (int i = 0; i < diffs.Length; i++)
            diffs[i] = time[i + 1] - time[i];
        return Median(diffs);
    }
}
