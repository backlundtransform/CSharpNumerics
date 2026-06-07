using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Detects the fundamental pitch (F0) of a signal using autocorrelation
    /// and Harmonic Product Spectrum (HPS).
    /// </summary>
    public class PitchDetector
    {
        /// <summary>Detection method.</summary>
        public enum Method
        {
            Autocorrelation,
            HarmonicProductSpectrum
        }

        /// <summary>FFT size for analysis.</summary>
        public int FftSize { get; }

        /// <summary>Minimum detectable frequency in Hz.</summary>
        public double MinFrequency { get; set; } = 50;

        /// <summary>Maximum detectable frequency in Hz.</summary>
        public double MaxFrequency { get; set; } = 2000;

        /// <summary>
        /// Creates a new PitchDetector.
        /// </summary>
        /// <param name="fftSize">FFT size (power of 2). Default 4096.</param>
        public PitchDetector(int fftSize = 4096)
        {
            int p = 1;
            while (p < fftSize) p <<= 1;
            FftSize = p;
        }

        /// <summary>
        /// Detect the fundamental frequency from an array of samples.
        /// </summary>
        /// <param name="samples">Audio samples.</param>
        /// <param name="sampleRate">Sample rate in Hz.</param>
        /// <param name="method">Detection method. Default: Autocorrelation.</param>
        /// <returns>Detected frequency in Hz, or 0 if no pitch found.</returns>
        public double Detect(double[] samples, int sampleRate, Method method = Method.Autocorrelation)
        {
            switch (method)
            {
                case Method.Autocorrelation:
                    return DetectAutocorrelation(samples, sampleRate);
                case Method.HarmonicProductSpectrum:
                    return DetectHPS(samples, sampleRate);
                default:
                    return 0;
            }
        }

        /// <summary>
        /// Detect pitch from an AudioBuffer.
        /// </summary>
        public double Detect(AudioBuffer buffer, Method method = Method.Autocorrelation)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            return Detect(mono.Samples, mono.SampleRate, method);
        }

        private double DetectAutocorrelation(double[] samples, int sampleRate)
        {
            int n = Math.Min(samples.Length, FftSize);
            int minLag = (int)(sampleRate / MaxFrequency);
            int maxLag = Math.Min(n / 2, (int)(sampleRate / MinFrequency));

            if (maxLag <= minLag) return 0;

            double bestCorr = double.NegativeInfinity;
            int bestLag = 0;

            for (int lag = minLag; lag <= maxLag; lag++)
            {
                double sum = 0;
                for (int i = 0; i < n - lag; i++)
                    sum += samples[i] * samples[i + lag];

                if (sum > bestCorr)
                {
                    bestCorr = sum;
                    bestLag = lag;
                }
            }

            if (bestLag == 0) return 0;
            return (double)sampleRate / bestLag;
        }

        private double DetectHPS(double[] samples, int sampleRate)
        {
            // Windowed FFT
            var windowed = new List<double>(FftSize);
            for (int i = 0; i < FftSize; i++)
            {
                double s = i < samples.Length ? samples[i] : 0.0;
                windowed.Add(s * SpectrumAnalyzer.ApplyWindow(SpectrumAnalyzer.WindowType.Hann, i, FftSize));
            }

            var spectrum = windowed.FastFourierTransform();
            int halfN = FftSize / 2;

            // Magnitude spectrum
            var mag = new double[halfN];
            for (int i = 0; i < halfN; i++)
                mag[i] = spectrum[i].GetMagnitude();

            // Harmonic product: multiply downsampled spectra (3 harmonics)
            int harmonics = 3;
            var hps = new double[halfN / harmonics];
            for (int i = 0; i < hps.Length; i++)
            {
                hps[i] = mag[i];
                for (int h = 2; h <= harmonics; h++)
                    hps[i] *= mag[i * h];
            }

            // Find peak in valid frequency range
            double freqRes = (double)sampleRate / FftSize;
            int minBin = Math.Max(1, (int)(MinFrequency / freqRes));
            int maxBin = Math.Min(hps.Length - 1, (int)(MaxFrequency / freqRes));

            double bestVal = 0;
            int bestBin = 0;
            for (int i = minBin; i <= maxBin; i++)
            {
                if (hps[i] > bestVal)
                {
                    bestVal = hps[i];
                    bestBin = i;
                }
            }

            return bestBin * freqRes;
        }
    }
}
