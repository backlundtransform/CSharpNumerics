using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Dynamic range compressor.
    /// Reduces gain when signal exceeds a threshold.
    /// </summary>
    public class Compressor
    {
        /// <summary>Threshold in linear amplitude (0-1). Signals above this are compressed.</summary>
        public double Threshold { get; set; }

        /// <summary>Compression ratio (e.g. 4 means 4:1 compression above threshold).</summary>
        public double Ratio { get; set; }

        /// <summary>Attack time in seconds (how fast gain reduction starts).</summary>
        public double Attack { get; set; }

        /// <summary>Release time in seconds (how fast gain returns to normal).</summary>
        public double Release { get; set; }

        /// <summary>Make-up gain applied after compression (linear).</summary>
        public double MakeupGain { get; set; }

        /// <summary>
        /// Creates a new compressor.
        /// </summary>
        public Compressor(double threshold = 0.5, double ratio = 4.0,
            double attack = 0.01, double release = 0.1, double makeupGain = 1.0)
        {
            Threshold = threshold;
            Ratio = ratio;
            Attack = attack;
            Release = release;
            MakeupGain = makeupGain;
        }

        /// <summary>
        /// Apply compression to an AudioBuffer (mono). Returns a new buffer.
        /// </summary>
        public AudioBuffer Process(AudioBuffer buffer)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            var input = mono.Samples;
            int n = input.Length;
            var output = new double[n];

            double envelope = 0;
            double attackCoeff = Attack > 0
                ? 1.0 - Math.Exp(-1.0 / (Attack * buffer.SampleRate))
                : 1.0;
            double releaseCoeff = Release > 0
                ? 1.0 - Math.Exp(-1.0 / (Release * buffer.SampleRate))
                : 1.0;

            for (int i = 0; i < n; i++)
            {
                double abs = Math.Abs(input[i]);

                // Smooth envelope follower
                if (abs > envelope)
                    envelope += attackCoeff * (abs - envelope);
                else
                    envelope += releaseCoeff * (abs - envelope);

                // Compute gain
                double gain;
                if (envelope <= Threshold)
                {
                    gain = 1.0;
                }
                else
                {
                    // Above threshold: compress
                    double excess = envelope - Threshold;
                    double compressed = Threshold + excess / Ratio;
                    gain = compressed / envelope;
                }

                output[i] = input[i] * gain * MakeupGain;
            }

            return new AudioBuffer(output, buffer.SampleRate, 1);
        }
    }
}
