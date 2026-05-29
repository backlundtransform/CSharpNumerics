using System;

namespace CSharpNumerics.Numerics.SignalProcessing;

/// <summary>
/// Zero-phase digital filtering via forward-backward (filtfilt) application.
/// Applies the filter forward and then backward, eliminating phase distortion.
/// The effective filter order is doubled. Suitable for offline analysis only.
/// </summary>
public static class ZeroPhaseFiltFilt
{
    /// <summary>
    /// Applies a Butterworth IIR filter with zero phase distortion.
    /// The signal is filtered forward and then backward through the biquad cascade.
    /// </summary>
    /// <param name="filter">The Butterworth filter to apply.</param>
    /// <param name="signal">Input signal.</param>
    /// <returns>Zero-phase filtered signal.</returns>
    public static double[] Apply(ButterworthFilter filter, double[] signal)
    {
        if (filter == null) throw new ArgumentNullException(nameof(filter));
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (signal.Length < 3) throw new ArgumentException("Signal must have at least 3 samples.", nameof(signal));

        int n = signal.Length;
        int padLen = Math.Min(3 * filter.Order, n - 1);

        // Extend signal with reflected boundaries to reduce transients
        var extended = new double[n + 2 * padLen];
        double first = signal[0];
        double last = signal[n - 1];

        // Left reflection: 2*signal[0] - signal[padLen], ..., 2*signal[0] - signal[1]
        for (int i = 0; i < padLen; i++)
        {
            extended[i] = 2.0 * first - signal[padLen - i];
        }

        // Original signal
        Array.Copy(signal, 0, extended, padLen, n);

        // Right reflection: 2*signal[n-1] - signal[n-2], ..., 2*signal[n-1] - signal[n-1-padLen]
        for (int i = 0; i < padLen; i++)
        {
            extended[n + padLen + i] = 2.0 * last - signal[n - 2 - i];
        }

        // Forward pass
        var forward = ApplyForward(filter, extended);

        // Reverse
        Array.Reverse(forward);

        // Backward pass
        var backward = ApplyForward(filter, forward);

        // Reverse back
        Array.Reverse(backward);

        // Extract original segment
        var result = new double[n];
        Array.Copy(backward, padLen, result, 0, n);

        return result;
    }

    /// <summary>
    /// Applies an FIR filter with zero phase distortion.
    /// </summary>
    /// <param name="filter">The FIR filter to apply.</param>
    /// <param name="signal">Input signal.</param>
    /// <returns>Zero-phase filtered signal.</returns>
    public static double[] Apply(FIRFilter filter, double[] signal)
    {
        if (filter == null) throw new ArgumentNullException(nameof(filter));
        if (signal == null) throw new ArgumentNullException(nameof(signal));
        if (signal.Length < 3) throw new ArgumentException("Signal must have at least 3 samples.", nameof(signal));

        int n = signal.Length;
        int padLen = Math.Min(3 * (filter.Order + 1), n - 1);

        // Extend signal with reflected boundaries
        var extended = new double[n + 2 * padLen];
        double first = signal[0];
        double last = signal[n - 1];

        for (int i = 0; i < padLen; i++)
        {
            extended[i] = 2.0 * first - signal[padLen - i];
        }

        Array.Copy(signal, 0, extended, padLen, n);

        for (int i = 0; i < padLen; i++)
        {
            extended[n + padLen + i] = 2.0 * last - signal[n - 2 - i];
        }

        // Forward pass
        var forward = filter.Apply(extended);

        // Reverse
        Array.Reverse(forward);

        // Backward pass
        var backward = filter.Apply(forward);

        // Reverse back
        Array.Reverse(backward);

        // Extract original segment
        var result = new double[n];
        Array.Copy(backward, padLen, result, 0, n);

        return result;
    }

    private static double[] ApplyForward(ButterworthFilter filter, double[] signal)
    {
        var output = new double[signal.Length];
        Array.Copy(signal, output, signal.Length);

        foreach (var section in filter.Sections)
        {
            double z1 = 0.0, z2 = 0.0;

            // Initialize state to steady-state for first sample to reduce transients
            double x0 = output[0];
            z1 = x0 * (section.B0 + section.B1 - section.A1 * section.B0);
            z2 = x0 * (section.B0 + section.B1 + section.B2 - (section.A1 + section.A2) * section.B0);

            // Simple approach: just start from zero state with extended signal
            z1 = 0.0;
            z2 = 0.0;

            for (int i = 0; i < output.Length; i++)
            {
                double input = output[i];
                double y = section.B0 * input + z1;
                z1 = section.B1 * input - section.A1 * y + z2;
                z2 = section.B2 * input - section.A2 * y;
                output[i] = y;
            }
        }

        return output;
    }
}
