using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Waves
{
    /// <summary>
    /// Analytically evaluates the superposition of travelling sinusoidal waves:
    /// u(x,t) = Σ Aᵢ sin(ωᵢt − kᵢx + φᵢ).
    /// <para>
    /// No PDE solver is required — pure closed-form evaluation.
    /// Useful for interference patterns, beat frequencies, and constructing initial conditions.
    /// </para>
    /// </summary>
    public class WaveSuperposition
    {
        #region Types

        private readonly struct Harmonic
        {
            public readonly double Amplitude;
            public readonly double AngularFrequency;
            public readonly double WaveNumber;
            public readonly double Phase;

            public Harmonic(double a, double w, double k, double phi)
            {
                Amplitude = a;
                AngularFrequency = w;
                WaveNumber = k;
                Phase = phi;
            }
        }

        #endregion

        #region Fields

        private readonly List<Harmonic> _harmonics = new List<Harmonic>();

        #endregion

        #region Adding Components

        /// <summary>
        /// Adds a harmonic component Asin(ωt − kx + φ).
        /// </summary>
        /// <param name="amplitude">Peak amplitude A.</param>
        /// <param name="angularFrequency">Angular frequency ω (rad/s).</param>
        /// <param name="waveNumber">Wave number k (rad/m).</param>
        /// <param name="phase">Phase offset φ (rad). Default 0.</param>
        public void AddHarmonic(double amplitude, double angularFrequency,
            double waveNumber, double phase = 0.0)
        {
            _harmonics.Add(new Harmonic(amplitude, angularFrequency, waveNumber, phase));
        }

        /// <summary>Number of harmonic components.</summary>
        public int Count => _harmonics.Count;

        #endregion

        #region Evaluation

        /// <summary>
        /// Evaluates the superposed field at a single point (x, t).
        /// </summary>
        public double Evaluate(double x, double t)
        {
            double sum = 0;
            for (int i = 0; i < _harmonics.Count; i++)
            {
                var h = _harmonics[i];
                sum += h.Amplitude * Math.Sin(h.AngularFrequency * t - h.WaveNumber * x + h.Phase);
            }
            return sum;
        }

        /// <summary>
        /// Evaluates the field at multiple x-positions at a fixed time.
        /// </summary>
        /// <param name="xValues">Spatial positions.</param>
        /// <param name="t">Time (s).</param>
        public List<Serie> EvaluateRange(double[] xValues, double t)
        {
            var result = new List<Serie>(xValues.Length);
            for (int i = 0; i < xValues.Length; i++)
                result.Add(new Serie { Index = xValues[i], Value = Evaluate(xValues[i], t) });
            return result;
        }

        #endregion

        #region Analysis

        /// <summary>
        /// Beat frequency |ω₁ − ω₂| / 2π for a two-component superposition.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown if fewer than 2 harmonics.</exception>
        public double BeatFrequency()
        {
            if (_harmonics.Count < 2)
                throw new InvalidOperationException("Beat frequency requires at least two harmonics.");

            return Math.Abs(_harmonics[0].AngularFrequency - _harmonics[1].AngularFrequency) / (2.0 * Math.PI);
        }

        /// <summary>
        /// Classifies each spatial sample as constructive (|u| > threshold × max amplitude sum)
        /// or destructive interference. Returns the field at each x with the classification in <see cref="Serie.Index"/> = x, <see cref="Serie.Value"/> = u(x,t).
        /// </summary>
        /// <param name="xValues">Spatial positions.</param>
        /// <param name="t">Time (s).</param>
        public List<Serie> InterferencePattern(double[] xValues, double t)
        {
            return EvaluateRange(xValues, t);
        }

        /// <summary>
        /// Computes the Fourier coefficients of the superposed field sampled at the given positions.
        /// </summary>
        /// <param name="xValues">Spatial sample positions (should be uniformly spaced).</param>
        /// <param name="t">Time (s).</param>
        public List<Serie> FourierCoefficients(double[] xValues, double t)
        {
            var samples = new List<double>(xValues.Length);
            for (int i = 0; i < xValues.Length; i++)
                samples.Add(Evaluate(xValues[i], t));

            // Pad to next power of two
            int nFft = 1;
            while (nFft < samples.Count) nFft <<= 1;
            while (samples.Count < nFft) samples.Add(0.0);

            double dx = xValues.Length > 1 ? xValues[1] - xValues[0] : 1.0;
            double samplingFrequency = 1.0 / dx;

            return samples.FastFourierTransform().ToFrequencyResolution(samplingFrequency);
        }

        #endregion
    }
}
