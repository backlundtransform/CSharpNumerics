using System;

namespace CSharpNumerics.Numerics.SignalProcessing;

/// <summary>
/// Finite Impulse Response (FIR) filter.
/// Applies a set of filter coefficients via linear convolution.
/// FIR filters are inherently stable and can provide exact linear phase.
/// </summary>
public class FIRFilter
{
    private readonly double[] _coefficients;

    /// <summary>
    /// Creates an FIR filter with the specified coefficients.
    /// </summary>
    /// <param name="coefficients">Filter tap coefficients h[0], h[1], ..., h[N-1].</param>
    public FIRFilter(double[] coefficients)
    {
        if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));
        if (coefficients.Length == 0) throw new ArgumentException("Coefficients array must not be empty.", nameof(coefficients));

        _coefficients = (double[])coefficients.Clone();
    }

    /// <summary>
    /// Gets the filter order (number of taps - 1).
    /// </summary>
    public int Order => _coefficients.Length - 1;

    /// <summary>
    /// Gets a copy of the filter coefficients.
    /// </summary>
    public double[] Coefficients => (double[])_coefficients.Clone();

    /// <summary>
    /// Applies the FIR filter to the input signal.
    /// Output length equals input length (zero-padded at boundaries).
    /// </summary>
    /// <param name="signal">Input signal.</param>
    /// <returns>Filtered signal of the same length.</returns>
    public double[] Apply(double[] signal)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));

        int n = signal.Length;
        int m = _coefficients.Length;
        var output = new double[n];

        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < m; j++)
            {
                int idx = i - j;
                if (idx >= 0 && idx < n)
                {
                    sum += _coefficients[j] * signal[idx];
                }
            }
            output[i] = sum;
        }

        return output;
    }

    /// <summary>
    /// Applies the FIR filter with symmetric boundary extension (mirror).
    /// Reduces boundary artifacts.
    /// </summary>
    /// <param name="signal">Input signal.</param>
    /// <returns>Filtered signal of the same length.</returns>
    public double[] ApplySymmetric(double[] signal)
    {
        if (signal == null) throw new ArgumentNullException(nameof(signal));

        int n = signal.Length;
        int m = _coefficients.Length;
        var output = new double[n];

        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < m; j++)
            {
                int idx = i - j;
                // Mirror extension
                if (idx < 0) idx = -idx;
                if (idx >= n) idx = 2 * (n - 1) - idx;
                idx = Math.Max(0, Math.Min(n - 1, idx));
                sum += _coefficients[j] * signal[idx];
            }
            output[i] = sum;
        }

        return output;
    }

    /// <summary>
    /// Computes the magnitude response at the specified normalized frequencies.
    /// </summary>
    /// <param name="normalizedFrequencies">Frequencies in [0, 1] where 1 = Nyquist frequency.</param>
    /// <returns>Magnitude response in linear scale.</returns>
    public double[] FrequencyResponse(double[] normalizedFrequencies)
    {
        if (normalizedFrequencies == null) throw new ArgumentNullException(nameof(normalizedFrequencies));

        var response = new double[normalizedFrequencies.Length];
        int m = _coefficients.Length;

        for (int f = 0; f < normalizedFrequencies.Length; f++)
        {
            double omega = Math.PI * normalizedFrequencies[f];
            double real = 0.0, imag = 0.0;

            for (int k = 0; k < m; k++)
            {
                real += _coefficients[k] * Math.Cos(omega * k);
                imag -= _coefficients[k] * Math.Sin(omega * k);
            }

            response[f] = Math.Sqrt(real * real + imag * imag);
        }

        return response;
    }
}
