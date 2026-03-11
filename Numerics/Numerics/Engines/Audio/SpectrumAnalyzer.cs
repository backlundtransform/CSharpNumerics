using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Real-time spectrum analyzer with windowing (Hann, Hamming, Blackman).
    /// Performs FFT on windowed audio frames and returns frequency/magnitude data.
    /// </summary>
    public class SpectrumAnalyzer
    {
        /// <summary>Available window functions.</summary>
        public enum WindowType
        {
            Rectangular,
            Hann,
            Hamming,
            Blackman
        }

        /// <summary>FFT size (must be power of 2).</summary>
        public int FftSize { get; }

        /// <summary>Window function to apply before FFT.</summary>
        public WindowType Window { get; set; }

        /// <summary>
        /// Creates a SpectrumAnalyzer with the given FFT size.
        /// </summary>
        /// <param name="fftSize">FFT size (power of 2). Default 1024.</param>
        /// <param name="window">Window type. Default Hann.</param>
        public SpectrumAnalyzer(int fftSize = 1024, WindowType window = WindowType.Hann)
        {
            // Ensure power of 2
            int p = 1;
            while (p < fftSize) p <<= 1;
            FftSize = p;
            Window = window;
        }

        /// <summary>
        /// Analyze a frame of samples and return magnitude spectrum as (frequency, magnitude) pairs.
        /// </summary>
        /// <param name="samples">Input samples (at least FftSize).</param>
        /// <param name="sampleRate">Sample rate in Hz.</param>
        /// <param name="offset">Starting offset into the samples array.</param>
        /// <returns>List of (frequency Hz, magnitude) up to Nyquist.</returns>
        public List<Serie> Analyze(double[] samples, int sampleRate, int offset = 0)
        {
            var windowed = new List<double>(FftSize);
            for (int i = 0; i < FftSize; i++)
            {
                int idx = offset + i;
                double s = idx < samples.Length ? samples[idx] : 0.0;
                windowed.Add(s * WindowFunction(i, FftSize));
            }

            var spectrum = windowed.FastFourierTransform();

            int halfN = FftSize / 2;
            double freqRes = (double)sampleRate / FftSize;
            var result = new List<Serie>(halfN);

            for (int i = 0; i <= halfN; i++)
            {
                double mag = spectrum[i].GetMagnitude() / FftSize;
                if (i > 0 && i < halfN) mag *= 2; // mirror compensation
                result.Add(new Serie { Index = i * freqRes, Value = mag });
            }

            return result;
        }

        /// <summary>
        /// Analyze an entire AudioBuffer and return the averaged magnitude spectrum.
        /// Overlapping frames (50% overlap) are averaged.
        /// </summary>
        public List<Serie> AnalyzeBuffer(AudioBuffer buffer)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            int hop = FftSize / 2;
            int frames = 0;
            double[] accumulator = null;

            for (int offset = 0; offset + FftSize <= mono.FrameCount; offset += hop)
            {
                var frame = Analyze(mono.Samples, mono.SampleRate, offset);
                if (accumulator == null)
                    accumulator = new double[frame.Count];
                for (int i = 0; i < frame.Count; i++)
                    accumulator[i] += frame[i].Value;
                frames++;
            }

            var result = new List<Serie>();
            if (accumulator != null && frames > 0)
            {
                double freqRes = (double)buffer.SampleRate / FftSize;
                for (int i = 0; i < accumulator.Length; i++)
                    result.Add(new Serie { Index = i * freqRes, Value = accumulator[i] / frames });
            }

            return result;
        }

        /// <summary>
        /// Generate a window function coefficient for sample index i in a window of size N.
        /// </summary>
        public static double ApplyWindow(WindowType type, int i, int n)
        {
            switch (type)
            {
                case WindowType.Hann:
                    return 0.5 * (1.0 - Math.Cos(2.0 * Math.PI * i / (n - 1)));
                case WindowType.Hamming:
                    return 0.54 - 0.46 * Math.Cos(2.0 * Math.PI * i / (n - 1));
                case WindowType.Blackman:
                    return 0.42 - 0.5 * Math.Cos(2.0 * Math.PI * i / (n - 1))
                           + 0.08 * Math.Cos(4.0 * Math.PI * i / (n - 1));
                default:
                    return 1.0;
            }
        }

        private double WindowFunction(int i, int n) => ApplyWindow(Window, i, n);
    }
}
