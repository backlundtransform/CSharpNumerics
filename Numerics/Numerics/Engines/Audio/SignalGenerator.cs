using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Generates standard audio waveforms into an <see cref="AudioBuffer"/>.
    /// Supported waveforms: Sine, Square, Sawtooth, Triangle, WhiteNoise.
    /// </summary>
    public static class SignalGenerator
    {
        /// <summary>Available waveform shapes.</summary>
        public enum Waveform
        {
            Sine,
            Square,
            Sawtooth,
            Triangle,
            WhiteNoise
        }

        /// <summary>
        /// Generate a mono AudioBuffer containing the specified waveform.
        /// </summary>
        /// <param name="waveform">Waveform type.</param>
        /// <param name="frequency">Frequency in Hz.</param>
        /// <param name="amplitude">Peak amplitude (0-1).</param>
        /// <param name="duration">Duration in seconds.</param>
        /// <param name="sampleRate">Sample rate in Hz.</param>
        /// <param name="phase">Initial phase in radians.</param>
        public static AudioBuffer Generate(Waveform waveform, double frequency, double amplitude,
            double duration, int sampleRate = 44100, double phase = 0.0)
        {
            var buffer = new AudioBuffer(sampleRate, 1, duration);
            FillBuffer(buffer.Samples, waveform, frequency, amplitude, sampleRate, phase);
            return buffer;
        }

        /// <summary>
        /// Fill an existing sample array with the specified waveform.
        /// </summary>
        public static void FillBuffer(double[] samples, Waveform waveform, double frequency,
            double amplitude, int sampleRate, double phase = 0.0)
        {
            var rng = waveform == Waveform.WhiteNoise ? new Random(42) : null;

            for (int i = 0; i < samples.Length; i++)
            {
                double t = (double)i / sampleRate;
                samples[i] = amplitude * Sample(waveform, frequency, t, phase, rng);
            }
        }

        private static double Sample(Waveform waveform, double freq, double t, double phase, Random rng)
        {
            double p = freq * t + phase / (2.0 * Math.PI);
            double frac = p - Math.Floor(p); // normalized phase [0,1)

            switch (waveform)
            {
                case Waveform.Sine:
                    return Math.Sin(2.0 * Math.PI * p);
                case Waveform.Square:
                    return frac < 0.5 ? 1.0 : -1.0;
                case Waveform.Sawtooth:
                    return 2.0 * frac - 1.0;
                case Waveform.Triangle:
                    return frac < 0.5
                        ? 4.0 * frac - 1.0
                        : 3.0 - 4.0 * frac;
                case Waveform.WhiteNoise:
                    return rng.NextDouble() * 2.0 - 1.0;
                default:
                    return 0.0;
            }
        }
    }
}
