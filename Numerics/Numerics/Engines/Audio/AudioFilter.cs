using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Frequency-domain audio filters: LowPass, HighPass, BandPass.
    /// Uses FFT → mask → IFFT for precise cutoff control.
    /// </summary>
    public static class AudioFilter
    {
        /// <summary>Filter type.</summary>
        public enum FilterType
        {
            LowPass,
            HighPass,
            BandPass
        }

        /// <summary>
        /// Apply a frequency-domain filter to an AudioBuffer (mono).
        /// </summary>
        /// <param name="buffer">Input audio buffer (mono).</param>
        /// <param name="filterType">Type of filter.</param>
        /// <param name="cutoffLow">Low cutoff frequency in Hz (used by LowPass and BandPass).</param>
        /// <param name="cutoffHigh">High cutoff frequency in Hz (used by HighPass and BandPass).</param>
        /// <returns>New filtered AudioBuffer.</returns>
        public static AudioBuffer Apply(AudioBuffer buffer, FilterType filterType,
            double cutoffLow, double cutoffHigh = 0)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            int n = mono.FrameCount;

            // Pad to next power of 2
            int nFft = 1;
            while (nFft < n) nFft <<= 1;

            var samples = new List<double>(nFft);
            for (int i = 0; i < n; i++) samples.Add(mono.Samples[i]);
            while (samples.Count < nFft) samples.Add(0.0);

            // Forward FFT
            var spectrum = samples.FastFourierTransform();

            // Build mask
            double freqResolution = (double)buffer.SampleRate / nFft;
            for (int i = 0; i < spectrum.Count; i++)
            {
                double freq = i * freqResolution;
                if (i > nFft / 2) freq = (nFft - i) * freqResolution; // mirror

                bool pass;
                switch (filterType)
                {
                    case FilterType.LowPass:
                        pass = freq <= cutoffLow;
                        break;
                    case FilterType.HighPass:
                        pass = freq >= cutoffHigh;
                        break;
                    case FilterType.BandPass:
                        pass = freq >= cutoffLow && freq <= cutoffHigh;
                        break;
                    default:
                        pass = true;
                        break;
                }

                if (!pass)
                    spectrum[i] = new ComplexNumber(0, 0);
            }

            // Inverse FFT
            var filtered = spectrum.InverseFastFourierTransform();

            // Extract real part, trim to original length
            var result = new double[n];
            for (int i = 0; i < n; i++)
                result[i] = filtered[i].realPart;

            return new AudioBuffer(result, buffer.SampleRate, 1);
        }
    }
}
