using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Waves
{
    /// <summary>
    /// Propagates a Gaussian wave packet through a dispersive medium using FFT.
    /// <para>
    /// The initial wave packet in k-space is a Gaussian centred at <see cref="CenterK"/>
    /// with width σ. Propagation is performed spectrally:
    /// Û(k,t) = Û(k,0) · exp(−iω(k)t), then inverse FFT to recover u(x,t).
    /// </para>
    /// </summary>
    public class WavePacket
    {
        #region Fields

        private readonly double[] _x;
        private readonly double _centerK;
        private readonly double _sigma;
        private readonly Func<double, double> _dispersion;
        private readonly int _nFft;
        private readonly double _dx;
        private readonly List<ComplexNumber> _initialSpectrum;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a Gaussian wave packet.
        /// </summary>
        /// <param name="x">Spatial grid points (must be uniformly spaced).</param>
        /// <param name="centerK">Central wave number k₀ (rad/m).</param>
        /// <param name="sigma">Spatial width σ of the Gaussian envelope (m).</param>
        /// <param name="dispersion">Dispersion relation ω(k).</param>
        public WavePacket(double[] x, double centerK, double sigma,
            Func<double, double> dispersion)
        {
            if (x == null || x.Length < 4)
                throw new ArgumentException("Need at least 4 grid points.", nameof(x));
            if (sigma <= 0) throw new ArgumentException("Sigma must be positive.", nameof(sigma));
            if (dispersion == null) throw new ArgumentNullException(nameof(dispersion));

            _x = x;
            _centerK = centerK;
            _sigma = sigma;
            _dispersion = dispersion;
            _dx = x[1] - x[0];

            // Pad to power of 2
            _nFft = 1;
            while (_nFft < x.Length) _nFft <<= 1;

            // Build initial spatial signal: Gaussian envelope × carrier
            var signal = new List<ComplexNumber>(_nFft);
            for (int i = 0; i < x.Length; i++)
            {
                double xi = x[i];
                double xc = x[x.Length / 2]; // centre of spatial domain
                double env = Math.Exp(-0.5 * ((xi - xc) / sigma) * ((xi - xc) / sigma));
                double re = env * Math.Cos(centerK * xi);
                double im = env * Math.Sin(centerK * xi);
                signal.Add(new ComplexNumber(re, im));
            }

            // Zero-pad
            while (signal.Count < _nFft)
                signal.Add(new ComplexNumber(0, 0));

            // Forward FFT to get Û(k, 0)
            _initialSpectrum = signal.FastFourierTransform();
        }

        #endregion

        #region Properties

        /// <summary>Central wave number k₀ (rad/m).</summary>
        public double CenterK => _centerK;

        /// <summary>Spatial width of the Gaussian envelope (m).</summary>
        public double Sigma => _sigma;

        /// <summary>
        /// Phase velocity at the central wave number: v_p = ω(k₀) / k₀.
        /// </summary>
        public double PhaseVelocity =>
            _centerK == 0 ? 0 : _dispersion(_centerK) / _centerK;

        /// <summary>
        /// Group velocity at the central wave number: v_g = dω/dk evaluated at k₀.
        /// </summary>
        public double GroupVelocity
        {
            get
            {
                // Central difference: dω/dk ≈ (ω(k+h) - ω(k-h)) / 2h
                double h = 1e-5;
                return (_dispersion(_centerK + h) - _dispersion(_centerK - h)) / (2.0 * h);
            }
        }

        #endregion

        #region Propagation

        /// <summary>
        /// Propagates the wave packet to time <paramref name="t"/> and returns the displacement
        /// at each original grid point.
        /// </summary>
        /// <param name="t">Time (s).</param>
        public List<Serie> Propagate(double t)
        {
            double dk = 2.0 * Math.PI / (_nFft * _dx);

            // Apply phase rotation: Û(k,t) = Û(k,0) · exp(-iω(k)t)
            var spectrum = new List<ComplexNumber>(_nFft);
            for (int j = 0; j < _nFft; j++)
            {
                // Map index to wave number (accounting for FFT ordering)
                double k = j <= _nFft / 2 ? j * dk : (j - _nFft) * dk;
                double omega = _dispersion(k);
                double phase = -omega * t;

                var rotation = new ComplexNumber(Math.Cos(phase), Math.Sin(phase));
                spectrum.Add(_initialSpectrum[j] * rotation);
            }

            // Inverse FFT
            var spatial = spectrum.InverseFastFourierTransform();

            // Extract real part at original grid points
            var result = new List<Serie>(_x.Length);
            for (int i = 0; i < _x.Length; i++)
                result.Add(new Serie { Index = _x[i], Value = spatial[i].realPart });

            return result;
        }

        /// <summary>
        /// Measures the spatial width (standard deviation) of |u(x,t)|² at time <paramref name="t"/>.
        /// </summary>
        /// <param name="t">Time (s).</param>
        public double Width(double t)
        {
            var field = Propagate(t);

            // Compute intensity-weighted mean and variance
            double totalWeight = 0;
            double meanX = 0;

            for (int i = 0; i < field.Count; i++)
            {
                double intensity = field[i].Value * field[i].Value;
                totalWeight += intensity;
                meanX += field[i].Index * intensity;
            }

            if (totalWeight < 1e-30) return 0;
            meanX /= totalWeight;

            double variance = 0;
            for (int i = 0; i < field.Count; i++)
            {
                double intensity = field[i].Value * field[i].Value;
                double dx = field[i].Index - meanX;
                variance += dx * dx * intensity;
            }

            return Math.Sqrt(variance / totalWeight);
        }

        /// <summary>
        /// Estimates the rate at which the packet spreads: (σ(t₂) − σ(t₁)) / (t₂ − t₁).
        /// </summary>
        /// <param name="t1">First time sample (s).</param>
        /// <param name="t2">Second time sample (s).</param>
        public double SpreadRate(double t1, double t2)
        {
            if (Math.Abs(t2 - t1) < 1e-15)
                throw new ArgumentException("Times must differ.");
            return (Width(t2) - Width(t1)) / (t2 - t1);
        }

        #endregion
    }
}
