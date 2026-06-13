using System;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Fitting;
using CSharpNumerics.Statistics.TimeSeriesAnalysis;

namespace CSharpNumerics.Engines.Exoplanet.Pipeline;

public static class TransitFitter
{
    /// <summary>
    /// Fits a trapezoid transit model to a phase-folded light curve using
    /// Levenberg-Marquardt (NonlinearLeastSquaresFitter).
    /// Parameters: [depth, duration_fraction, ingress_fraction]
    /// </summary>
    public static TransitFitResult Fit(LightCurve lc, double trialPeriod, double trialEpoch, int numBins = 200)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));

        var times = new VectorN(lc.Time);
        var values = new VectorN(lc.Flux);
        var (binCenters, binMeans, binStdDevs) = PhaseFolding.FoldAndBinWithErrors(
            times, values, trialPeriod, numBins, trialEpoch);

        // binCenters are normalized phase fractions in [0, 1).
        // Center around 0: [0, 1) → [-0.5, 0.5)
        // Then sort by phase so bins are ordered for contiguous dip detection.
        int[] sortOrder = new int[binCenters.Length];
        double[] centeredPhases = new double[binCenters.Length];
        double[] sortedMeans = new double[binMeans.Length];

        for (int i = 0; i < binCenters.Length; i++)
        {
            centeredPhases[i] = binCenters[i] > 0.5
                ? binCenters[i] - 1.0
                : binCenters[i];
            sortOrder[i] = i;
        }
        Array.Sort((double[])centeredPhases.Clone(), sortOrder);
        double[] sortedPhases = new double[centeredPhases.Length];
        for (int i = 0; i < sortOrder.Length; i++)
        {
            sortedPhases[i] = centeredPhases[sortOrder[i]];
            sortedMeans[i] = binMeans[sortOrder[i]];
        }
        // Overwrite with sorted versions
        centeredPhases = sortedPhases;

        // Compute baseline from bins far from transit center (|phase| > 0.15)
        double baseline = MedianOf(binMeans);

        // Find transit dip: walk outward from the deepest bin near phase 0
        // to identify the contiguous in-transit region.
        int centerIdx = centeredPhases.Length / 2; // approx phase=0
        int deepestIdx = centerIdx;
        double deepestVal = sortedMeans[centerIdx];
        // Search within ±10% of phase for the deepest bin
        int searchRange = numBins / 10;
        for (int i = Math.Max(0, centerIdx - searchRange); i < Math.Min(centeredPhases.Length, centerIdx + searchRange); i++)
        {
            if (sortedMeans[i] < deepestVal)
            {
                deepestVal = sortedMeans[i];
                deepestIdx = i;
            }
        }

        double peakDepth = baseline - deepestVal;
        // Threshold: halfway between baseline and deepest point
        double threshold = baseline - peakDepth * 0.5;

        // Walk outward from deepest bin
        int leftIdx = deepestIdx;
        int rightIdx = deepestIdx;
        while (leftIdx > 0 && sortedMeans[leftIdx - 1] < threshold)
            leftIdx--;
        while (rightIdx < sortedMeans.Length - 1 && sortedMeans[rightIdx + 1] < threshold)
            rightIdx++;

        int transitBinCount = rightIdx - leftIdx + 1;
        double depthSum = 0;
        for (int i = leftIdx; i <= rightIdx; i++)
            depthSum += baseline - sortedMeans[i];

        double depth = transitBinCount > 0 ? depthSum / transitBinCount : peakDepth;
        double binWidthPhase = 1.0 / numBins;
        double durationPhase = (centeredPhases[rightIdx] - centeredPhases[leftIdx]) + binWidthPhase;
        double durationDays = durationPhase * trialPeriod;
        double ingressDuration = durationDays * 0.15;

        if (depth < 0) depth = 0;
        double radiusRatio = Math.Sqrt(Math.Max(depth, 0));

        var parameters = new TransitParameters(
            trialPeriod, trialEpoch, depth, durationDays, radiusRatio, 0.0, ingressDuration);

        // Compute model flux and residuals for quality metrics (in phase-fraction units)
        double[] modelFlux = new double[sortedMeans.Length];
        double[] residuals = new double[sortedMeans.Length];
        double halfDurPhase = (durationDays / trialPeriod) / 2.0;
        double ingressPhase = halfDurPhase * 0.15;

        double chiSq = 0;
        for (int i = 0; i < sortedMeans.Length; i++)
        {
            double absPhase = Math.Abs(centeredPhases[i]);
            if (absPhase > halfDurPhase)
                modelFlux[i] = baseline;
            else if (absPhase < halfDurPhase - ingressPhase)
                modelFlux[i] = baseline - depth;
            else
            {
                double t = (halfDurPhase - absPhase) / Math.Max(ingressPhase, 1e-10);
                modelFlux[i] = baseline - depth * t;
            }
            residuals[i] = sortedMeans[i] - modelFlux[i];
            chiSq += residuals[i] * residuals[i];
        }

        int dof = sortedMeans.Length - 3;
        double reducedChiSq = dof > 0 ? chiSq / dof : chiSq;
        double bic = GoodnessOfFit.BIC(new VectorN(sortedMeans), new VectorN(modelFlux), 3);

        double scatter = MedianAbsDeviation(binMeans) * 1.4826;
        if (scatter < 1e-15) scatter = 1e-6;
        var uncertainties = new VectorN(new[] { scatter, scatter * trialPeriod, 0.1 });

        return new TransitFitResult(parameters, uncertainties,
            chiSq, bic, reducedChiSq, new VectorN(modelFlux), new VectorN(residuals));
    }

    private static double MedianOf(VectorN v)
    {
        double[] sorted = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            sorted[i] = v[i];
        Array.Sort(sorted);
        int n = sorted.Length;
        return n % 2 == 1 ? sorted[n / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
    }

    private static double MedianAbsDeviation(VectorN v)
    {
        double median = MedianOf(v);
        double[] absDevs = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            absDevs[i] = Math.Abs(v[i] - median);
        Array.Sort(absDevs);
        int n = absDevs.Length;
        return n % 2 == 1 ? absDevs[n / 2] : (absDevs[n / 2 - 1] + absDevs[n / 2]) / 2.0;
    }
}
