using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Additive synthesizer: combines multiple <see cref="AudioOscillator"/>s with
    /// individual gains, applies an optional <see cref="Envelope"/>, and renders to an <see cref="AudioBuffer"/>.
    /// </summary>
    public class Synthesizer
    {
        private readonly List<(AudioOscillator Oscillator, double Gain)> _voices = new();

        /// <summary>Optional ADSR envelope applied to the mixed output.</summary>
        public Envelope Envelope { get; set; }

        /// <summary>Master gain (default 1.0).</summary>
        public double MasterGain { get; set; } = 1.0;

        /// <summary>Sample rate in Hz.</summary>
        public int SampleRate { get; set; } = 44100;

        /// <summary>
        /// Add an oscillator voice with a relative gain.
        /// </summary>
        public void AddVoice(AudioOscillator oscillator, double gain = 1.0)
        {
            _voices.Add((oscillator, gain));
        }

        /// <summary>Number of active voices.</summary>
        public int VoiceCount => _voices.Count;

        /// <summary>
        /// Render the synthesizer output to an AudioBuffer.
        /// </summary>
        /// <param name="duration">Duration in seconds.</param>
        /// <param name="noteOffTime">
        /// Time in seconds when the note is released (for envelope).
        /// Default: no release (sustain forever).
        /// </param>
        public AudioBuffer Render(double duration, double noteOffTime = double.MaxValue)
        {
            var buffer = new AudioBuffer(SampleRate, 1, duration);

            for (int i = 0; i < buffer.FrameCount; i++)
            {
                double sample = 0;
                foreach (var (osc, gain) in _voices)
                    sample += osc.NextSample(SampleRate) * gain;

                buffer.Samples[i] = sample * MasterGain;
            }

            if (Envelope != null)
                Envelope.Apply(buffer, noteOffTime);

            return buffer;
        }

        /// <summary>Reset all oscillator phases.</summary>
        public void Reset()
        {
            foreach (var (osc, _) in _voices)
                osc.Reset();
        }
    }
}
