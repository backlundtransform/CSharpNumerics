using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Delay effect with feedback using a circular buffer.
    /// </summary>
    public class Delay
    {
        /// <summary>Delay time in seconds.</summary>
        public double DelayTime { get; }

        /// <summary>Feedback gain (0 = single echo, 0.9 = long tail). Range [0, 1).</summary>
        public double Feedback { get; set; }

        /// <summary>Wet/dry mix (0 = dry, 1 = fully wet).</summary>
        public double WetMix { get; set; }

        /// <summary>
        /// Creates a new delay effect.
        /// </summary>
        /// <param name="delayTime">Delay time in seconds.</param>
        /// <param name="feedback">Feedback amount [0, 1).</param>
        /// <param name="wetMix">Wet/dry mix [0, 1].</param>
        public Delay(double delayTime, double feedback = 0.5, double wetMix = 0.5)
        {
            if (delayTime <= 0) throw new ArgumentOutOfRangeException(nameof(delayTime));
            if (feedback < 0 || feedback >= 1) throw new ArgumentOutOfRangeException(nameof(feedback));

            DelayTime = delayTime;
            Feedback = feedback;
            WetMix = wetMix;
        }

        /// <summary>
        /// Apply delay to an AudioBuffer (mono). Returns a new buffer.
        /// </summary>
        public AudioBuffer Process(AudioBuffer buffer)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            var input = mono.Samples;
            int n = input.Length;
            int delaySamples = (int)(DelayTime * buffer.SampleRate);
            if (delaySamples < 1) delaySamples = 1;

            var delayBuf = new double[delaySamples];
            int writeIdx = 0;
            var output = new double[n];

            for (int i = 0; i < n; i++)
            {
                double delayed = delayBuf[writeIdx];
                double inputWithFeedback = input[i] + delayed * Feedback;
                delayBuf[writeIdx] = inputWithFeedback;
                output[i] = input[i] * (1.0 - WetMix) + delayed * WetMix;
                writeIdx = (writeIdx + 1) % delaySamples;
            }

            return new AudioBuffer(output, buffer.SampleRate, 1);
        }
    }
}
