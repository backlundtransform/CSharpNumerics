using System;

namespace CSharpNumerics.Engines.Audio
{
    /// <summary>
    /// Spatial audio processing: panning and distance attenuation.
    /// </summary>
    public static class SpatialAudio
    {
        /// <summary>
        /// Apply stereo panning to a mono buffer.
        /// </summary>
        /// <param name="buffer">Mono input buffer.</param>
        /// <param name="pan">Pan position: -1 = full left, 0 = center, +1 = full right.</param>
        /// <returns>Stereo AudioBuffer.</returns>
        public static AudioBuffer Pan(AudioBuffer buffer, double pan)
        {
            var mono = buffer.Channels == 1 ? buffer : buffer.ToMono();
            pan = Math.Max(-1, Math.Min(1, pan));

            // Constant-power panning
            double angle = (pan + 1.0) * Math.PI / 4.0; // [0, π/2]
            double gainL = Math.Cos(angle);
            double gainR = Math.Sin(angle);

            var stereo = new AudioBuffer(mono.SampleRate, 2, mono.Duration);
            for (int i = 0; i < mono.FrameCount; i++)
            {
                stereo[i, 0] = mono.Samples[i] * gainL;
                stereo[i, 1] = mono.Samples[i] * gainR;
            }

            return stereo;
        }

        /// <summary>
        /// Apply distance attenuation using inverse-distance law.
        /// </summary>
        /// <param name="buffer">Input buffer.</param>
        /// <param name="distance">Distance from source in arbitrary units.</param>
        /// <param name="referenceDistance">Reference distance where gain = 1.</param>
        /// <param name="maxDistance">Maximum distance beyond which attenuation is clamped.</param>
        /// <returns>New attenuated AudioBuffer.</returns>
        public static AudioBuffer AttenuateByDistance(AudioBuffer buffer, double distance,
            double referenceDistance = 1.0, double maxDistance = 100.0)
        {
            if (distance < 0) throw new ArgumentOutOfRangeException(nameof(distance));

            distance = Math.Max(distance, referenceDistance);
            distance = Math.Min(distance, maxDistance);

            // Inverse-distance attenuation
            double gain = referenceDistance / distance;

            var output = new double[buffer.Samples.Length];
            for (int i = 0; i < output.Length; i++)
                output[i] = buffer.Samples[i] * gain;

            return new AudioBuffer(output, buffer.SampleRate, buffer.Channels);
        }

        /// <summary>
        /// Combine panning and distance attenuation for a source at a given position.
        /// </summary>
        /// <param name="buffer">Mono input buffer.</param>
        /// <param name="sourceX">Source X position.</param>
        /// <param name="sourceY">Source Y position (depth).</param>
        /// <param name="listenerX">Listener X position.</param>
        /// <param name="listenerY">Listener Y position.</param>
        /// <param name="referenceDistance">Reference distance.</param>
        /// <returns>Stereo buffer with spatial positioning.</returns>
        public static AudioBuffer Spatialize(AudioBuffer buffer,
            double sourceX, double sourceY, double listenerX, double listenerY,
            double referenceDistance = 1.0)
        {
            double dx = sourceX - listenerX;
            double dy = sourceY - listenerY;
            double distance = Math.Sqrt(dx * dx + dy * dy);

            // Pan from relative X: normalize to [-1, 1] range
            double pan = distance > 0 ? dx / distance : 0;

            var panned = Pan(buffer, pan);
            return AttenuateByDistance(panned, distance, referenceDistance);
        }
    }
}
