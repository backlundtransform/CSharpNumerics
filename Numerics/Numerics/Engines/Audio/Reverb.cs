using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Schroeder reverb model using parallel comb filters and series all-pass filters.
    /// </summary>
    public class Reverb
    {
        /// <summary>Room size parameter (controls comb filter delays). Range 0-1.</summary>
        public double RoomSize { get; set; }

        /// <summary>Damping parameter (controls high-frequency absorption). Range 0-1.</summary>
        public double Damping { get; set; }

        /// <summary>Wet/dry mix (0 = dry, 1 = fully wet).</summary>
        public double WetMix { get; set; }

        // Comb filter delay lengths (in samples at 44100 Hz), Schroeder standard
        private static readonly int[] CombDelays = { 1557, 1617, 1491, 1422, 1277, 1356, 1188, 1116 };
        private static readonly int[] AllPassDelays = { 225, 556, 441, 341 };

        /// <summary>
        /// Creates a new reverb effect.
        /// </summary>
        /// <param name="roomSize">Room size (0-1). Default 0.5.</param>
        /// <param name="damping">Damping (0-1). Default 0.5.</param>
        /// <param name="wetMix">Wet/dry mix (0-1). Default 0.3.</param>
        public Reverb(double roomSize = 0.5, double damping = 0.5, double wetMix = 0.3)
        {
            RoomSize = roomSize;
            Damping = damping;
            WetMix = wetMix;
        }

        /// <summary>
        /// Apply reverb to an AudioBuffer (mono). Returns a new buffer.
        /// </summary>
        public AudioBuffer Process(AudioBuffer buffer)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            var input = mono.Samples;
            int n = input.Length;

            // Scale delays by sample rate ratio
            double srRatio = (double)buffer.SampleRate / 44100.0;

            // Parallel comb filters
            double[] combSum = new double[n];
            double feedback = 0.7 + 0.28 * RoomSize; // [0.7, 0.98]

            for (int c = 0; c < CombDelays.Length; c++)
            {
                int delay = Math.Max(1, (int)(CombDelays[c] * srRatio));
                var combOut = CombFilter(input, delay, feedback, Damping);
                for (int i = 0; i < n; i++)
                    combSum[i] += combOut[i];
            }

            // Normalize comb sum
            double combScale = 1.0 / CombDelays.Length;
            for (int i = 0; i < n; i++)
                combSum[i] *= combScale;

            // Series all-pass filters
            double[] current = combSum;
            for (int a = 0; a < AllPassDelays.Length; a++)
            {
                int delay = Math.Max(1, (int)(AllPassDelays[a] * srRatio));
                current = AllPassFilter(current, delay, 0.5);
            }

            // Mix wet/dry
            var output = new double[n];
            for (int i = 0; i < n; i++)
                output[i] = input[i] * (1.0 - WetMix) + current[i] * WetMix;

            return new AudioBuffer(output, buffer.SampleRate, 1);
        }

        private static double[] CombFilter(double[] input, int delay, double feedback, double damping)
        {
            int n = input.Length;
            var output = new double[n];
            var buf = new double[delay];
            double filterStore = 0;
            int bufIdx = 0;

            for (int i = 0; i < n; i++)
            {
                double bufOut = buf[bufIdx];
                filterStore = bufOut * (1.0 - damping) + filterStore * damping;
                buf[bufIdx] = input[i] + filterStore * feedback;
                output[i] = bufOut;
                bufIdx = (bufIdx + 1) % delay;
            }

            return output;
        }

        private static double[] AllPassFilter(double[] input, int delay, double gain)
        {
            int n = input.Length;
            var output = new double[n];
            var buf = new double[delay];
            int bufIdx = 0;

            for (int i = 0; i < n; i++)
            {
                double bufOut = buf[bufIdx];
                buf[bufIdx] = input[i] + bufOut * gain;
                output[i] = bufOut - input[i] * gain;
                bufIdx = (bufIdx + 1) % delay;
            }

            return output;
        }
    }
}
