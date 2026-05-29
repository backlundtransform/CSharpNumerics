using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.SignalProcessing;

/// <summary>
/// Filter type for IIR filter design.
/// </summary>
public enum FilterType
{
    Lowpass,
    Highpass,
    Bandpass
}

/// <summary>
/// Static factory for designing digital IIR and FIR filters.
/// Uses the bilinear transform method for analog-to-digital conversion.
/// </summary>
public static class FilterDesign
{
    /// <summary>
    /// Designs a Butterworth lowpass filter.
    /// </summary>
    /// <param name="order">Filter order (number of poles). Must be ≥ 1.</param>
    /// <param name="cutoffFrequency">Cutoff frequency in Hz.</param>
    /// <param name="sampleRate">Sample rate in Hz.</param>
    /// <returns>A ButterworthFilter instance.</returns>
    public static ButterworthFilter DesignLowpass(int order, double cutoffFrequency, double sampleRate)
    {
        ValidateParameters(order, cutoffFrequency, sampleRate);
        double normalizedCutoff = cutoffFrequency / (sampleRate / 2.0);
        return DesignButterworth(order, normalizedCutoff, FilterType.Lowpass);
    }

    /// <summary>
    /// Designs a Butterworth highpass filter.
    /// </summary>
    /// <param name="order">Filter order (number of poles). Must be ≥ 1.</param>
    /// <param name="cutoffFrequency">Cutoff frequency in Hz.</param>
    /// <param name="sampleRate">Sample rate in Hz.</param>
    /// <returns>A ButterworthFilter instance.</returns>
    public static ButterworthFilter DesignHighpass(int order, double cutoffFrequency, double sampleRate)
    {
        ValidateParameters(order, cutoffFrequency, sampleRate);
        double normalizedCutoff = cutoffFrequency / (sampleRate / 2.0);
        return DesignButterworth(order, normalizedCutoff, FilterType.Highpass);
    }

    /// <summary>
    /// Designs a Butterworth bandpass filter.
    /// The resulting filter order is 2 × the specified order.
    /// </summary>
    /// <param name="order">Filter order per edge. Resulting filter is 2×order.</param>
    /// <param name="lowCutoff">Lower cutoff frequency in Hz.</param>
    /// <param name="highCutoff">Upper cutoff frequency in Hz.</param>
    /// <param name="sampleRate">Sample rate in Hz.</param>
    /// <returns>A ButterworthFilter instance.</returns>
    public static ButterworthFilter DesignBandpass(int order, double lowCutoff, double highCutoff, double sampleRate)
    {
        if (order < 1) throw new ArgumentException("Order must be at least 1.", nameof(order));
        if (lowCutoff <= 0 || highCutoff <= 0)
            throw new ArgumentException("Cutoff frequencies must be positive.");
        if (lowCutoff >= highCutoff)
            throw new ArgumentException("Low cutoff must be less than high cutoff.");
        double nyquist = sampleRate / 2.0;
        if (highCutoff >= nyquist)
            throw new ArgumentException("High cutoff must be less than Nyquist frequency.");

        double normalizedLow = lowCutoff / nyquist;
        double normalizedHigh = highCutoff / nyquist;
        return DesignButterworthBandpass(order, normalizedLow, normalizedHigh);
    }

    /// <summary>
    /// Designs a windowed-sinc FIR lowpass filter.
    /// </summary>
    /// <param name="numTaps">Number of filter taps (must be odd for Type I).</param>
    /// <param name="cutoffFrequency">Cutoff frequency in Hz.</param>
    /// <param name="sampleRate">Sample rate in Hz.</param>
    /// <returns>A FIRFilter instance.</returns>
    public static FIRFilter DesignFIRLowpass(int numTaps, double cutoffFrequency, double sampleRate)
    {
        if (numTaps < 3) throw new ArgumentException("Number of taps must be at least 3.", nameof(numTaps));
        if (numTaps % 2 == 0) numTaps++; // Ensure odd for symmetric filter

        double normalizedCutoff = cutoffFrequency / (sampleRate / 2.0);
        if (normalizedCutoff <= 0 || normalizedCutoff >= 1.0)
            throw new ArgumentException("Cutoff must be between 0 and Nyquist.");

        double fc = normalizedCutoff * 0.5; // Fraction of sample rate
        var coeffs = new double[numTaps];
        int halfTaps = numTaps / 2;

        for (int i = 0; i < numTaps; i++)
        {
            int n = i - halfTaps;
            if (n == 0)
            {
                coeffs[i] = 2.0 * fc;
            }
            else
            {
                coeffs[i] = Math.Sin(2.0 * Math.PI * fc * n) / (Math.PI * n);
            }

            // Apply Hamming window
            coeffs[i] *= 0.54 - 0.46 * Math.Cos(2.0 * Math.PI * i / (numTaps - 1));
        }

        // Normalize for unity gain at DC
        double sum = 0.0;
        for (int i = 0; i < numTaps; i++) sum += coeffs[i];
        for (int i = 0; i < numTaps; i++) coeffs[i] /= sum;

        return new FIRFilter(coeffs);
    }

    /// <summary>
    /// Designs a windowed-sinc FIR highpass filter.
    /// </summary>
    /// <param name="numTaps">Number of filter taps (must be odd).</param>
    /// <param name="cutoffFrequency">Cutoff frequency in Hz.</param>
    /// <param name="sampleRate">Sample rate in Hz.</param>
    /// <returns>A FIRFilter instance.</returns>
    public static FIRFilter DesignFIRHighpass(int numTaps, double cutoffFrequency, double sampleRate)
    {
        if (numTaps < 3) throw new ArgumentException("Number of taps must be at least 3.", nameof(numTaps));
        if (numTaps % 2 == 0) numTaps++;

        // Design lowpass then spectrally invert
        var lowpass = DesignFIRLowpass(numTaps, cutoffFrequency, sampleRate);
        var lpCoeffs = lowpass.Coefficients;
        var hpCoeffs = new double[lpCoeffs.Length];
        int halfTaps = lpCoeffs.Length / 2;

        for (int i = 0; i < lpCoeffs.Length; i++)
        {
            hpCoeffs[i] = -lpCoeffs[i];
        }
        hpCoeffs[halfTaps] += 1.0; // Spectral inversion

        return new FIRFilter(hpCoeffs);
    }

    private static ButterworthFilter DesignButterworth(int order, double normalizedCutoff, FilterType filterType)
    {
        // Pre-warp the cutoff frequency for bilinear transform
        double warpedCutoff = Math.Tan(Math.PI * normalizedCutoff * 0.5);

        var sections = new List<BiquadSection>();
        double gain = 1.0;

        int numSections = (order + 1) / 2;

        for (int k = 0; k < numSections; k++)
        {
            double theta;
            bool isSinglePole = (order % 2 == 1) && (k == numSections - 1);

            if (isSinglePole)
            {
                // First-order section (real pole)
                // Analog prototype: H(s) = 1 / (s + 1)
                // Bilinear transform with frequency warping
                double a = warpedCutoff;

                double b0, b1, b2, a1, a2;

                if (filterType == FilterType.Lowpass)
                {
                    // H(z) after bilinear transform of 1/(s/wc + 1)
                    double alpha = a / (1.0 + a);
                    b0 = alpha;
                    b1 = alpha;
                    b2 = 0.0;
                    a1 = -(1.0 - a) / (1.0 + a);
                    a2 = 0.0;
                }
                else // Highpass
                {
                    double alpha = 1.0 / (1.0 + a);
                    b0 = alpha;
                    b1 = -alpha;
                    b2 = 0.0;
                    a1 = -(1.0 - a) / (1.0 + a);
                    a2 = 0.0;
                }

                sections.Add(new BiquadSection(b0, b1, b2, a1, a2));
            }
            else
            {
                // Second-order section (complex conjugate pole pair)
                // Pole angle for Butterworth: θ_k = π(2k+1)/(2N)
                theta = Math.PI * (2.0 * k + 1.0) / (2.0 * order);
                double cosTheta = Math.Cos(theta);

                // Analog prototype pole: s = -sin(θ) ± j*cos(θ)
                // For normalized Butterworth: Q = 1/(2*sin(θ))
                double sinTheta = Math.Sin(theta);

                double wc = warpedCutoff;
                double wc2 = wc * wc;

                double b0, b1, b2, a1, a2;

                if (filterType == FilterType.Lowpass)
                {
                    // Bilinear transform of second-order analog lowpass section
                    double denom = 1.0 + 2.0 * sinTheta * wc + wc2;
                    b0 = wc2 / denom;
                    b1 = 2.0 * wc2 / denom;
                    b2 = wc2 / denom;
                    a1 = 2.0 * (wc2 - 1.0) / denom;
                    a2 = (1.0 - 2.0 * sinTheta * wc + wc2) / denom;
                }
                else // Highpass
                {
                    double denom = wc2 + 2.0 * sinTheta * wc + 1.0;
                    b0 = 1.0 / denom;
                    b1 = -2.0 / denom;
                    b2 = 1.0 / denom;
                    a1 = 2.0 * (1.0 - wc2) / denom;
                    a2 = (wc2 - 2.0 * sinTheta * wc + 1.0) / denom;
                }

                sections.Add(new BiquadSection(b0, b1, b2, a1, a2));
            }
        }

        return new ButterworthFilter(sections.ToArray(), gain);
    }

    private static ButterworthFilter DesignButterworthBandpass(int order, double normalizedLow, double normalizedHigh)
    {
        // Bandpass as cascade of highpass + lowpass
        // This is a simplified approach; for production, a proper bandpass transform would be used
        var lowpass = DesignButterworth(order, normalizedHigh, FilterType.Lowpass);
        var highpass = DesignButterworth(order, normalizedLow, FilterType.Highpass);

        // Combine sections from both filters
        var allSections = new List<BiquadSection>();
        allSections.AddRange(lowpass.Sections);
        allSections.AddRange(highpass.Sections);

        return new ButterworthFilter(allSections.ToArray(), 1.0);
    }

    private static void ValidateParameters(int order, double cutoffFrequency, double sampleRate)
    {
        if (order < 1) throw new ArgumentException("Order must be at least 1.", nameof(order));
        if (cutoffFrequency <= 0) throw new ArgumentException("Cutoff frequency must be positive.", nameof(cutoffFrequency));
        if (sampleRate <= 0) throw new ArgumentException("Sample rate must be positive.", nameof(sampleRate));
        if (cutoffFrequency >= sampleRate / 2.0)
            throw new ArgumentException("Cutoff frequency must be less than Nyquist frequency (sampleRate/2).");
    }
}
