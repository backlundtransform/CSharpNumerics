using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Represents a buffer of audio samples with metadata (sample rate, channels, duration).
    /// Supports mono and stereo. Samples are stored as doubles in [-1, 1].
    /// </summary>
    public class AudioBuffer
    {
        /// <summary>Raw sample data. For stereo, samples are interleaved: [L0, R0, L1, R1, …].</summary>
        public double[] Samples { get; }

        /// <summary>Sample rate in Hz (e.g. 44100).</summary>
        public int SampleRate { get; }

        /// <summary>Number of channels (1 = mono, 2 = stereo).</summary>
        public int Channels { get; }

        /// <summary>Duration in seconds.</summary>
        public double Duration => (double)FrameCount / SampleRate;

        /// <summary>Number of sample frames (= Samples.Length / Channels).</summary>
        public int FrameCount => Samples.Length / Channels;

        /// <summary>
        /// Creates a new AudioBuffer with the specified parameters.
        /// </summary>
        /// <param name="sampleRate">Sample rate in Hz.</param>
        /// <param name="channels">Number of channels (1 or 2).</param>
        /// <param name="durationSeconds">Duration in seconds.</param>
        public AudioBuffer(int sampleRate, int channels, double durationSeconds)
        {
            if (sampleRate <= 0) throw new ArgumentOutOfRangeException(nameof(sampleRate));
            if (channels < 1 || channels > 2) throw new ArgumentOutOfRangeException(nameof(channels));
            if (durationSeconds <= 0) throw new ArgumentOutOfRangeException(nameof(durationSeconds));

            SampleRate = sampleRate;
            Channels = channels;
            int frameCount = (int)(sampleRate * durationSeconds);
            Samples = new double[frameCount * channels];
        }

        /// <summary>
        /// Wraps existing sample data.
        /// </summary>
        public AudioBuffer(double[] samples, int sampleRate, int channels)
        {
            if (samples == null) throw new ArgumentNullException(nameof(samples));
            if (sampleRate <= 0) throw new ArgumentOutOfRangeException(nameof(sampleRate));
            if (channels < 1 || channels > 2) throw new ArgumentOutOfRangeException(nameof(channels));

            Samples = samples;
            SampleRate = sampleRate;
            Channels = channels;
        }

        /// <summary>
        /// Gets or sets a sample at the given frame and channel index.
        /// </summary>
        public double this[int frame, int channel]
        {
            get => Samples[frame * Channels + channel];
            set => Samples[frame * Channels + channel] = value;
        }

        /// <summary>
        /// Returns a mono mixdown. If already mono, returns this.
        /// </summary>
        public AudioBuffer ToMono()
        {
            if (Channels == 1) return this;
            var mono = new double[FrameCount];
            for (int i = 0; i < FrameCount; i++)
                mono[i] = (Samples[i * 2] + Samples[i * 2 + 1]) * 0.5;
            return new AudioBuffer(mono, SampleRate, 1);
        }

        /// <summary>
        /// Mix another buffer into this one (additive). Buffers must share sample rate and channels.
        /// </summary>
        public void MixIn(AudioBuffer other, double gain = 1.0)
        {
            if (other.SampleRate != SampleRate || other.Channels != Channels)
                throw new ArgumentException("Sample rate and channel count must match.");
            int len = Math.Min(Samples.Length, other.Samples.Length);
            for (int i = 0; i < len; i++)
                Samples[i] += other.Samples[i] * gain;
        }

        /// <summary>
        /// Clamp all samples to [-1, 1].
        /// </summary>
        public void Normalize()
        {
            double peak = 0;
            for (int i = 0; i < Samples.Length; i++)
            {
                double abs = Math.Abs(Samples[i]);
                if (abs > peak) peak = abs;
            }
            if (peak > 0 && peak != 1.0)
            {
                double scale = 1.0 / peak;
                for (int i = 0; i < Samples.Length; i++)
                    Samples[i] *= scale;
            }
        }
    }
}
