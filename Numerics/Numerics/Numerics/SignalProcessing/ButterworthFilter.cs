using System;

namespace CSharpNumerics.Numerics.SignalProcessing;

/// <summary>
/// Represents a second-order section (biquad) for IIR filter implementation.
/// Transfer function: H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
/// </summary>
public class BiquadSection
{
    public double B0 { get; }
    public double B1 { get; }
    public double B2 { get; }
    public double A1 { get; }
    public double A2 { get; }

    public BiquadSection(double b0, double b1, double b2, double a1, double a2)
    {
        B0 = b0;
        B1 = b1;
        B2 = b2;
        A1 = a1;
        A2 = a2;
    }

    /// <summary>
    /// Process a single sample through this biquad section using Direct Form II Transposed.
    /// </summary>
    internal double Process(double input, ref double z1, ref double z2)
    {
        double output = B0 * input + z1;
        z1 = B1 * input - A1 * output + z2;
        z2 = B2 * input - A2 * output;
        return output;
    }
}

/// <summary>
/// Butterworth IIR filter implemented as a cascade of second-order sections (biquads).
/// Supports lowpass, highpass, and bandpass configurations.
/// Maximally flat magnitude response in the passband.
/// </summary>
public class ButterworthFilter
{
    private readonly BiquadSection[] _sections;
    private readonly double _gain;

    /// <summary>
    /// Creates a Butterworth filter from pre-computed biquad sections.
    /// Use <see cref="FilterDesign"/> to create instances.
    /// </summary>
    internal ButterworthFilter(BiquadSection[] sections, double gain)
    {
        _sections = sections;
        _gain = gain;
    }

    /// <summary>
    /// Gets the filter order.
    /// </summary>
    public int Order => _sections.Length * 2;

    /// <summary>
    /// Gets the second-order sections.
    /// </summary>
    public BiquadSection[] Sections => (BiquadSection[])_sections.Clone();

    /// <summary>
    /// Applies the filter to the input signal (causal, forward-only).
    /// Introduces phase delay.
    /// </summary>
    /// <param name="signal">Input signal.</param>
    /// <returns>Filtered signal of the same length.</returns>
    public double[] Apply(double[] signal)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));

        var output = new double[signal.Length];
        Array.Copy(signal, output, signal.Length);

        // Apply gain
        for (int i = 0; i < output.Length; i++)
            output[i] *= _gain;

        // Cascade through each biquad section
        foreach (var section in _sections)
        {
            double z1 = 0.0, z2 = 0.0;
            for (int i = 0; i < output.Length; i++)
            {
                output[i] = section.Process(output[i], ref z1, ref z2);
            }
        }

        return output;
    }

    /// <summary>
    /// Computes the magnitude response at the specified normalized frequencies.
    /// </summary>
    /// <param name="normalizedFrequencies">Frequencies in range [0, 0.5] (fraction of Nyquist).</param>
    /// <returns>Magnitude response in linear scale.</returns>
    public double[] FrequencyResponse(double[] normalizedFrequencies)
    {
        if (normalizedFrequencies == null) throw new ArgumentNullException(nameof(normalizedFrequencies));

        var response = new double[normalizedFrequencies.Length];

        for (int f = 0; f < normalizedFrequencies.Length; f++)
        {
            double omega = 2.0 * Math.PI * normalizedFrequencies[f];
            double cosW = Math.Cos(omega);
            double cos2W = Math.Cos(2.0 * omega);
            double sinW = Math.Sin(omega);
            double sin2W = Math.Sin(2.0 * omega);

            double magSquared = _gain * _gain;

            foreach (var s in _sections)
            {
                // Numerator: B0 + B1*e^{-jw} + B2*e^{-j2w}
                double numReal = s.B0 + s.B1 * cosW + s.B2 * cos2W;
                double numImag = -(s.B1 * sinW + s.B2 * sin2W);

                // Denominator: 1 + A1*e^{-jw} + A2*e^{-j2w}
                double denReal = 1.0 + s.A1 * cosW + s.A2 * cos2W;
                double denImag = -(s.A1 * sinW + s.A2 * sin2W);

                double numMagSq = numReal * numReal + numImag * numImag;
                double denMagSq = denReal * denReal + denImag * denImag;

                magSquared *= numMagSq / denMagSq;
            }

            response[f] = Math.Sqrt(magSquared);
        }

        return response;
    }
}
