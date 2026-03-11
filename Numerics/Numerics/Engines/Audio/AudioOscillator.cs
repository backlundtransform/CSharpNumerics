using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Real-time audio oscillator with frequency and amplitude modulation support.
    /// Maintains internal phase to ensure continuous, click-free output.
    /// </summary>
    public class AudioOscillator
    {
        private double _phase;

        /// <summary>Waveform shape.</summary>
        public SignalGenerator.Waveform Waveform { get; set; }

        /// <summary>Frequency in Hz.</summary>
        public double Frequency { get; set; }

        /// <summary>Peak amplitude (0-1).</summary>
        public double Amplitude { get; set; }

        /// <summary>Current oscillator phase in radians.</summary>
        public double Phase => _phase;

        /// <summary>
        /// Creates a new audio oscillator.
        /// </summary>
        public AudioOscillator(SignalGenerator.Waveform waveform, double frequency, double amplitude = 1.0)
        {
            Waveform = waveform;
            Frequency = frequency;
            Amplitude = amplitude;
            _phase = 0;
        }

        /// <summary>
        /// Generate the next sample at the given sample rate, advancing internal phase.
        /// </summary>
        public double NextSample(int sampleRate)
        {
            double value = Amplitude * EvaluateWaveform(_phase);
            _phase += 2.0 * Math.PI * Frequency / sampleRate;
            if (_phase > 2.0 * Math.PI)
                _phase -= 2.0 * Math.PI;
            return value;
        }

        /// <summary>
        /// Fill a buffer with samples from this oscillator, advancing phase continuously.
        /// </summary>
        public AudioBuffer GenerateBuffer(double duration, int sampleRate = 44100)
        {
            var buffer = new AudioBuffer(sampleRate, 1, duration);
            for (int i = 0; i < buffer.FrameCount; i++)
                buffer.Samples[i] = NextSample(sampleRate);
            return buffer;
        }

        /// <summary>Reset phase to zero.</summary>
        public void Reset() => _phase = 0;

        private double EvaluateWaveform(double phase)
        {
            double frac = phase / (2.0 * Math.PI);
            frac -= Math.Floor(frac);

            switch (Waveform)
            {
                case SignalGenerator.Waveform.Sine:
                    return Math.Sin(phase);
                case SignalGenerator.Waveform.Square:
                    return frac < 0.5 ? 1.0 : -1.0;
                case SignalGenerator.Waveform.Sawtooth:
                    return 2.0 * frac - 1.0;
                case SignalGenerator.Waveform.Triangle:
                    return frac < 0.5
                        ? 4.0 * frac - 1.0
                        : 3.0 - 4.0 * frac;
                default:
                    return 0.0;
            }
        }
    }
}
