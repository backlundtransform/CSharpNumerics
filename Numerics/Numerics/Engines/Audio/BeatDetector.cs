using CSharpNumerics.Numerics;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Detects beats / onsets in an audio signal using energy-based onset detection.
    /// </summary>
    public class BeatDetector
    {
        /// <summary>FFT frame size for spectral analysis.</summary>
        public int FrameSize { get; }

        /// <summary>Hop size between frames (default: FrameSize / 2).</summary>
        public int HopSize { get; }

        /// <summary>
        /// Sensitivity threshold multiplier relative to mean spectral flux.
        /// Higher values → fewer detections. Default 1.5.
        /// </summary>
        public double Threshold { get; set; } = 1.5;

        /// <summary>
        /// Creates a new BeatDetector.
        /// </summary>
        /// <param name="frameSize">FFT frame size (power of 2). Default 1024.</param>
        /// <param name="hopSize">Hop size. Default: frameSize / 2.</param>
        public BeatDetector(int frameSize = 1024, int hopSize = 0)
        {
            int p = 1;
            while (p < frameSize) p <<= 1;
            FrameSize = p;
            HopSize = hopSize > 0 ? hopSize : FrameSize / 2;
        }

        /// <summary>
        /// Detect beat onset times in an AudioBuffer.
        /// </summary>
        /// <returns>List of onset times in seconds.</returns>
        public List<double> Detect(AudioBuffer buffer)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            return Detect(mono.Samples, mono.SampleRate);
        }

        /// <summary>
        /// Detect beat onset times from raw samples.
        /// </summary>
        public List<double> Detect(double[] samples, int sampleRate)
        {
            // Compute spectral flux for each frame
            var fluxValues = new List<double>();
            var fluxTimes = new List<double>();
            double[] prevMag = null;

            for (int offset = 0; offset + FrameSize <= samples.Length; offset += HopSize)
            {
                // Windowed FFT
                var windowed = new List<double>(FrameSize);
                for (int i = 0; i < FrameSize; i++)
                    windowed.Add(samples[offset + i] *
                        SpectrumAnalyzer.ApplyWindow(SpectrumAnalyzer.WindowType.Hann, i, FrameSize));

                var spectrum = windowed.FastFourierTransform();
                int halfN = FrameSize / 2;

                var mag = new double[halfN];
                for (int i = 0; i < halfN; i++)
                    mag[i] = spectrum[i].GetMagnitude();

                // Spectral flux: sum of positive differences
                double flux = 0;
                if (prevMag != null)
                {
                    for (int i = 0; i < halfN; i++)
                    {
                        double diff = mag[i] - prevMag[i];
                        if (diff > 0) flux += diff;
                    }
                }

                fluxValues.Add(flux);
                fluxTimes.Add((double)offset / sampleRate);
                prevMag = mag;
            }

            // Pick peaks above threshold * mean
            double mean = 0;
            foreach (var f in fluxValues) mean += f;
            if (fluxValues.Count > 0) mean /= fluxValues.Count;

            double thresh = mean * Threshold;
            var onsets = new List<double>();

            for (int i = 1; i < fluxValues.Count - 1; i++)
            {
                if (fluxValues[i] > thresh &&
                    fluxValues[i] > fluxValues[i - 1] &&
                    fluxValues[i] > fluxValues[i + 1])
                {
                    onsets.Add(fluxTimes[i]);
                }
            }

            return onsets;
        }

        /// <summary>
        /// Estimate tempo (BPM) from detected onsets.
        /// </summary>
        public double EstimateTempo(AudioBuffer buffer)
        {
            var onsets = Detect(buffer);
            return EstimateTempoFromOnsets(onsets);
        }

        /// <summary>
        /// Estimate tempo (BPM) from a list of onset times.
        /// </summary>
        public static double EstimateTempoFromOnsets(List<double> onsets)
        {
            if (onsets.Count < 2) return 0;

            // Average inter-onset interval
            double totalInterval = 0;
            for (int i = 1; i < onsets.Count; i++)
                totalInterval += onsets[i] - onsets[i - 1];

            double avgInterval = totalInterval / (onsets.Count - 1);
            if (avgInterval <= 0) return 0;

            return 60.0 / avgInterval;
        }
    }
}
