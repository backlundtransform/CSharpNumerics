using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Fitting;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Statistics.TimeSeriesAnalysis;

/// <summary>
/// Peak detection and fitting for time-series event analysis.
/// Supports Gaussian and Lorentzian profile fitting via Levenberg-Marquardt.
/// </summary>
public static class PeakFitting
{
    /// <summary>
    /// Find peaks (local maxima) in a 1-D signal with prominence and distance filtering.
    /// </summary>
    /// <param name="values">Signal values.</param>
    /// <param name="prominence">Minimum prominence (height above the higher of the two neighbouring troughs). 0 = no filtering.</param>
    /// <param name="minDistance">Minimum sample distance between peaks (default 1).</param>
    /// <returns>List of peak indices sorted by position.</returns>
    public static List<int> FindPeaks(VectorN values, double prominence = 0, int minDistance = 1)
    {
        if (values.Length < 3)
            return new List<int>();

        double[] y = values.Values;
        int n = y.Length;
        var candidates = new List<int>();

        // Step 1: find all local maxima
        for (int i = 1; i < n - 1; i++)
        {
            if (y[i] > y[i - 1] && y[i] >= y[i + 1])
                candidates.Add(i);
        }

        // Step 2: prominence filter
        if (prominence > 0)
        {
            var filtered = new List<int>();
            foreach (int idx in candidates)
            {
                double prom = ComputeProminence(y, idx, n);
                if (prom >= prominence)
                    filtered.Add(idx);
            }
            candidates = filtered;
        }

        // Step 3: distance filter — keep tallest peaks first
        if (minDistance > 1 && candidates.Count > 1)
        {
            // Sort by height desc
            candidates.Sort((a, b) => y[b].CompareTo(y[a]));
            var kept = new List<int>();
            foreach (int idx in candidates)
            {
                bool tooClose = false;
                foreach (int k in kept)
                {
                    if (Math.Abs(idx - k) < minDistance)
                    {
                        tooClose = true;
                        break;
                    }
                }
                if (!tooClose) kept.Add(idx);
            }
            kept.Sort();
            candidates = kept;
        }

        return candidates;
    }

    /// <summary>
    /// Find valleys (local minima) in a 1-D signal.
    /// </summary>
    public static List<int> FindValleys(VectorN values, double prominence = 0, int minDistance = 1)
    {
        double[] y = values.Values;
        double[] neg = new double[y.Length];
        for (int i = 0; i < y.Length; i++)
            neg[i] = -y[i];
        return FindPeaks(new VectorN(neg), prominence, minDistance);
    }

    /// <summary>
    /// Fit a Gaussian profile to a peak region: y = A · exp(−(x−μ)² / (2σ²)) + C.
    /// </summary>
    /// <param name="x">X-values of the region around the peak.</param>
    /// <param name="y">Y-values of the region.</param>
    /// <param name="centerEstimate">Initial guess for the peak centre.</param>
    /// <param name="widthEstimate">Initial guess for σ.</param>
    /// <param name="peakIndex">Index of the peak in the original full series (for bookkeeping). Default −1.</param>
    public static PeakResult FitGaussian(
        VectorN x, VectorN y,
        double centerEstimate, double widthEstimate,
        int peakIndex = -1)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("x and y must have the same length.");

        double ampEst = 0;
        int n = y.Length;
        for (int i = 0; i < n; i++)
            if (y[i] > ampEst) ampEst = y[i];

        double baseEst = Math.Min(y[0], y[n - 1]);
        ampEst -= baseEst;

        // model: p[0]=A, p[1]=μ, p[2]=σ, p[3]=baseline
        Func<double, VectorN, double> model = (xi, p) =>
            p[0] * Math.Exp(-0.5 * ((xi - p[1]) / p[2]) * ((xi - p[1]) / p[2])) + p[3];

        var init = new VectorN(new[] { ampEst, centerEstimate, widthEstimate, baseEst });
        FittingResult fit = NonlinearLeastSquaresFitter.Fit(model, x, y, init);

        return new PeakResult(
            peakIndex,
            fit.Coefficients[1],
            fit.Coefficients[0],
            Math.Abs(fit.Coefficients[2]),
            fit.Coefficients[3],
            PeakShape.Gaussian,
            fit.RSquared);
    }

    /// <summary>
    /// Fit a Lorentzian (Cauchy) profile: y = A · γ² / ((x−x₀)² + γ²) + C.
    /// </summary>
    public static PeakResult FitLorentzian(
        VectorN x, VectorN y,
        double centerEstimate, double widthEstimate,
        int peakIndex = -1)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("x and y must have the same length.");

        double ampEst = 0;
        int n = y.Length;
        for (int i = 0; i < n; i++)
            if (y[i] > ampEst) ampEst = y[i];

        double baseEst = Math.Min(y[0], y[n - 1]);
        ampEst -= baseEst;

        // model: p[0]=A, p[1]=x₀, p[2]=γ, p[3]=baseline
        Func<double, VectorN, double> model = (xi, p) =>
            p[0] * (p[2] * p[2]) / ((xi - p[1]) * (xi - p[1]) + p[2] * p[2]) + p[3];

        var init = new VectorN(new[] { ampEst, centerEstimate, widthEstimate, baseEst });
        FittingResult fit = NonlinearLeastSquaresFitter.Fit(model, x, y, init);

        return new PeakResult(
            peakIndex,
            fit.Coefficients[1],
            fit.Coefficients[0],
            Math.Abs(fit.Coefficients[2]),
            fit.Coefficients[3],
            PeakShape.Lorentzian,
            fit.RSquared);
    }

    /// <summary>
    /// Automatically find peaks and fit each one with the specified shape.
    /// A region of ±<paramref name="fitRadius"/> samples around each peak is used for fitting.
    /// </summary>
    /// <param name="x">Full x-axis values.</param>
    /// <param name="y">Full y-axis values.</param>
    /// <param name="prominence">Minimum peak prominence.</param>
    /// <param name="shape">Peak shape to fit.</param>
    /// <param name="fitRadius">Number of samples on each side of the peak to include in the fit. Default 10.</param>
    /// <param name="minDistance">Minimum distance between detected peaks.</param>
    public static List<PeakResult> FindAndFitPeaks(
        VectorN x, VectorN y,
        double prominence,
        PeakShape shape = PeakShape.Gaussian,
        int fitRadius = 10,
        int minDistance = 1)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("x and y must have the same length.");

        List<int> peaks = FindPeaks(y, prominence, minDistance);
        var results = new List<PeakResult>();
        int n = x.Length;
        double[] xArr = x.Values;
        double[] yArr = y.Values;

        foreach (int pk in peaks)
        {
            int lo = Math.Max(0, pk - fitRadius);
            int hi = Math.Min(n - 1, pk + fitRadius);
            int count = hi - lo + 1;
            if (count < 4) continue; // need enough points for 4 parameters

            double[] regionX = new double[count];
            double[] regionY = new double[count];
            for (int i = 0; i < count; i++)
            {
                regionX[i] = xArr[lo + i];
                regionY[i] = yArr[lo + i];
            }

            double centerEst = xArr[pk];
            // Width estimate: half-width at half-max approximation
            double halfMax = (yArr[pk] + Math.Min(regionY[0], regionY[count - 1])) / 2.0;
            double widthEst = 0;
            for (int i = pk - lo; i < count; i++)
            {
                if (regionY[i] < halfMax)
                {
                    widthEst = regionX[i] - centerEst;
                    break;
                }
            }
            if (widthEst <= 0) widthEst = (regionX[count - 1] - regionX[0]) / 4.0;

            PeakResult result = shape == PeakShape.Gaussian
                ? FitGaussian(new VectorN(regionX), new VectorN(regionY), centerEst, widthEst, pk)
                : FitLorentzian(new VectorN(regionX), new VectorN(regionY), centerEst, widthEst, pk);

            results.Add(result);
        }

        return results;
    }

    /// <summary>
    /// Compute prominence of a peak: height above the higher of the two 
    /// enclosing troughs found by walking left and right to either a higher peak or the boundary.
    /// </summary>
    private static double ComputeProminence(double[] y, int peakIdx, int n)
    {
        double peakVal = y[peakIdx];

        // Walk left to find the lowest point before reaching a higher peak or boundary
        double leftMin = peakVal;
        for (int i = peakIdx - 1; i >= 0; i--)
        {
            if (y[i] < leftMin) leftMin = y[i];
            if (y[i] > peakVal) break;
        }

        // Walk right
        double rightMin = peakVal;
        for (int i = peakIdx + 1; i < n; i++)
        {
            if (y[i] < rightMin) rightMin = y[i];
            if (y[i] > peakVal) break;
        }

        double referenceLine = Math.Max(leftMin, rightMin);
        return peakVal - referenceLine;
    }
}
