using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// ADSR (Attack-Decay-Sustain-Release) envelope generator.
    /// Returns an amplitude multiplier (0-1) based on elapsed time and note-off time.
    /// </summary>
    public class Envelope
    {
        /// <summary>Attack time in seconds.</summary>
        public double Attack { get; }

        /// <summary>Decay time in seconds.</summary>
        public double Decay { get; }

        /// <summary>Sustain level (0-1).</summary>
        public double Sustain { get; }

        /// <summary>Release time in seconds.</summary>
        public double Release { get; }

        /// <summary>
        /// Creates an ADSR envelope.
        /// </summary>
        public Envelope(double attack, double decay, double sustain, double release)
        {
            if (attack < 0) throw new ArgumentOutOfRangeException(nameof(attack));
            if (decay < 0) throw new ArgumentOutOfRangeException(nameof(decay));
            if (sustain < 0 || sustain > 1) throw new ArgumentOutOfRangeException(nameof(sustain));
            if (release < 0) throw new ArgumentOutOfRangeException(nameof(release));

            Attack = attack;
            Decay = decay;
            Sustain = sustain;
            Release = release;
        }

        /// <summary>Total duration of the ADS phases before release (Attack + Decay = minimum note-on time for full cycle).</summary>
        public double NoteOnDuration => Attack + Decay;

        /// <summary>
        /// Evaluate the envelope at time <paramref name="t"/> seconds after note-on.
        /// <paramref name="noteOffTime"/> is the time (seconds after note-on) when the key was released.
        /// Pass <see cref="double.MaxValue"/> if the note is still held.
        /// </summary>
        public double Evaluate(double t, double noteOffTime = double.MaxValue)
        {
            if (t < 0) return 0;

            double level;

            // ADS phase
            if (t < Attack)
            {
                // Attack: ramp from 0 to 1
                level = t / Attack;
            }
            else if (t < Attack + Decay)
            {
                // Decay: ramp from 1 to Sustain
                double decayProgress = (t - Attack) / Decay;
                level = 1.0 - (1.0 - Sustain) * decayProgress;
            }
            else
            {
                // Sustain
                level = Sustain;
            }

            // Release phase
            if (t >= noteOffTime)
            {
                double releaseElapsed = t - noteOffTime;
                if (Release <= 0 || releaseElapsed >= Release)
                    return 0;

                // Level at note-off
                double levelAtOff = EvaluateAtNoteOff(noteOffTime);
                level = levelAtOff * (1.0 - releaseElapsed / Release);
            }

            return Math.Max(0, level);
        }

        /// <summary>
        /// Applies the envelope to an AudioBuffer. The entire buffer represents the note duration;
        /// release begins at <paramref name="noteOffTime"/>.
        /// </summary>
        public void Apply(AudioBuffer buffer, double noteOffTime = double.MaxValue)
        {
            for (int i = 0; i < buffer.FrameCount; i++)
            {
                double t = (double)i / buffer.SampleRate;
                double env = Evaluate(t, noteOffTime);
                for (int ch = 0; ch < buffer.Channels; ch++)
                    buffer[i, ch] *= env;
            }
        }

        private double EvaluateAtNoteOff(double noteOffTime)
        {
            if (noteOffTime < Attack)
                return noteOffTime / Attack;
            if (noteOffTime < Attack + Decay)
            {
                double dp = (noteOffTime - Attack) / Decay;
                return 1.0 - (1.0 - Sustain) * dp;
            }
            return Sustain;
        }
    }
}
