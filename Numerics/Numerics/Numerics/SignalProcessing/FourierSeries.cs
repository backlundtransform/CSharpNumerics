using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Numerics.SignalProcessing
{
    /// <summary>
    /// Fourier series analysis and synthesis for periodic functions.
    /// <para>
    /// Given a periodic function f(t) with period T, computes the real Fourier coefficients
    /// a₀, aₙ, bₙ and can reconstruct the function from a finite number of terms.
    /// Also provides power spectrum and Parseval energy.
    /// </para>
    /// </summary>
    public class FourierSeries
    {
        #region Fields

        private double _a0;
        private readonly List<double> _an = new List<double>();
        private readonly List<double> _bn = new List<double>();
        private double _period;
        private int _nTerms;

        #endregion

        #region Properties

        /// <summary>Period T of the analysed function.</summary>
        public double Period => _period;

        /// <summary>Number of harmonic terms computed.</summary>
        public int NTerms => _nTerms;

        /// <summary>DC component a₀ (average value).</summary>
        public double A0 => _a0;

        /// <summary>Cosine coefficients aₙ for n = 1 … NTerms.</summary>
        public IReadOnlyList<double> An => _an;

        /// <summary>Sine coefficients bₙ for n = 1 … NTerms.</summary>
        public IReadOnlyList<double> Bn => _bn;

        #endregion

        #region Analysis

        /// <summary>
        /// Analyses a periodic function by computing real Fourier coefficients
        /// via numerical integration (midpoint rule for speed).
        /// </summary>
        /// <param name="f">The periodic function to analyse.</param>
        /// <param name="period">Period T.</param>
        /// <param name="nTerms">Number of harmonic terms (n = 1 … nTerms).</param>
        /// <param name="nSamples">Number of quadrature samples. Default: 1024.</param>
        public void Analyze(Func<double, double> f, double period, int nTerms, int nSamples = 1024)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));
            if (period <= 0) throw new ArgumentException("Period must be positive.", nameof(period));
            if (nTerms < 1) throw new ArgumentException("Need at least 1 term.", nameof(nTerms));
            if (nSamples < 4) throw new ArgumentException("Need at least 4 samples.", nameof(nSamples));

            _period = period;
            _nTerms = nTerms;
            _an.Clear();
            _bn.Clear();

            double dt = period / nSamples;
            double omega = 2.0 * Math.PI / period;

            // a₀ = (2/T) ∫₀ᵀ f(t) dt — but note: some conventions use 1/T (we use 2/T for consistency)
            double sumA0 = 0;
            for (int i = 0; i < nSamples; i++)
            {
                double t = (i + 0.5) * dt;
                sumA0 += f(t);
            }
            _a0 = (2.0 / nSamples) * sumA0;

            for (int n = 1; n <= nTerms; n++)
            {
                double sumAn = 0, sumBn = 0;
                for (int i = 0; i < nSamples; i++)
                {
                    double t = (i + 0.5) * dt;
                    double val = f(t);
                    double angle = n * omega * t;
                    sumAn += val * Math.Cos(angle);
                    sumBn += val * Math.Sin(angle);
                }
                _an.Add((2.0 / nSamples) * sumAn);
                _bn.Add((2.0 / nSamples) * sumBn);
            }
        }

        #endregion

        #region Synthesis

        /// <summary>
        /// Reconstructs the function at time <paramref name="t"/> using the computed coefficients:
        /// f̂(t) = a₀/2 + Σ [aₙcos(nωt) + bₙsin(nωt)].
        /// </summary>
        /// <param name="t">Time (s).</param>
        /// <returns>Approximated function value.</returns>
        public double Synthesize(double t)
        {
            if (_nTerms == 0)
                throw new InvalidOperationException("Call Analyze before Synthesize.");

            double omega = 2.0 * Math.PI / _period;
            double sum = _a0 / 2.0;

            for (int n = 0; n < _nTerms; n++)
            {
                double angle = (n + 1) * omega * t;
                sum += _an[n] * Math.Cos(angle) + _bn[n] * Math.Sin(angle);
            }

            return sum;
        }

        /// <summary>
        /// Synthesizes the signal at multiple time points using at most <paramref name="maxTerms"/> terms.
        /// Useful for demonstrating Gibbs phenomenon by varying the term count.
        /// </summary>
        /// <param name="tValues">Time values to evaluate at.</param>
        /// <param name="maxTerms">Maximum number of terms to include (clamped to NTerms).</param>
        public List<Serie> SynthesizeRange(double[] tValues, int maxTerms = int.MaxValue)
        {
            if (_nTerms == 0)
                throw new InvalidOperationException("Call Analyze before Synthesize.");

            int terms = Math.Min(maxTerms, _nTerms);
            double omega = 2.0 * Math.PI / _period;

            var result = new List<Serie>(tValues.Length);
            for (int j = 0; j < tValues.Length; j++)
            {
                double t = tValues[j];
                double sum = _a0 / 2.0;
                for (int n = 0; n < terms; n++)
                {
                    double angle = (n + 1) * omega * t;
                    sum += _an[n] * Math.Cos(angle) + _bn[n] * Math.Sin(angle);
                }
                result.Add(new Serie { Index = t, Value = sum });
            }

            return result;
        }

        #endregion

        #region Spectrum & Energy

        /// <summary>
        /// Returns the power spectrum |cₙ|² = (aₙ² + bₙ²) / 4 for n = 1 … NTerms.
        /// Index is the harmonic number n.
        /// </summary>
        public List<Serie> PowerSpectrum()
        {
            if (_nTerms == 0)
                throw new InvalidOperationException("Call Analyze first.");

            var result = new List<Serie>(_nTerms);
            for (int n = 0; n < _nTerms; n++)
            {
                double power = (_an[n] * _an[n] + _bn[n] * _bn[n]) / 4.0;
                result.Add(new Serie { Index = n + 1, Value = power });
            }
            return result;
        }

        /// <summary>
        /// Parseval energy estimate in the frequency domain:
        /// E = (a₀²/4) + Σ (aₙ² + bₙ²) / 4.
        /// <para>
        /// Parseval's theorem states that this equals (1/T)∫₀ᵀ f(t)² dt,
        /// i.e. energy in time domain = energy in frequency domain.
        /// </para>
        /// </summary>
        public double ParsevalEnergy()
        {
            if (_nTerms == 0)
                throw new InvalidOperationException("Call Analyze first.");

            double energy = _a0 * _a0 / 4.0;
            for (int n = 0; n < _nTerms; n++)
                energy += (_an[n] * _an[n] + _bn[n] * _bn[n]) / 2.0;
            return energy;
        }

        /// <summary>
        /// Time-domain energy: (1/T) ∫₀ᵀ f(t)² dt, computed numerically.
        /// Should match <see cref="ParsevalEnergy"/> (Parseval's theorem).
        /// </summary>
        /// <param name="f">The original function.</param>
        /// <param name="nSamples">Number of quadrature samples.</param>
        public double TimeDomainEnergy(Func<double, double> f, int nSamples = 1024)
        {
            double dt = _period / nSamples;
            double sum = 0;
            for (int i = 0; i < nSamples; i++)
            {
                double t = (i + 0.5) * dt;
                double val = f(t);
                sum += val * val;
            }
            return sum / nSamples;
        }

        #endregion
    }
}
