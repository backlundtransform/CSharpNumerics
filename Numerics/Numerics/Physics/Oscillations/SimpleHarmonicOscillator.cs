using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Oscillations
{
    /// <summary>
    /// Models a simple (undamped) harmonic oscillator: ẍ + ω₀²x = 0.
    /// <para>
    /// The system is integrated with Velocity Verlet (symplectic, O(dt²)),
    /// which conserves energy over long simulations — ideal for undamped systems.
    /// </para>
    /// <para>
    /// All analytic helpers are provided for verification:
    /// exact solution, energy, frequency spectrum (via FFT), and phase portrait.
    /// </para>
    /// </summary>
    public class SimpleHarmonicOscillator : IOscillator
    {
        #region Fields

        private readonly double _mass;
        private readonly double _stiffness;
        private readonly double _x0;
        private readonly double _v0;

        private double _x;
        private double _v;
        private double _t;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a simple harmonic oscillator with the given mass and spring stiffness.
        /// </summary>
        /// <param name="mass">Mass in kg (must be &gt; 0).</param>
        /// <param name="stiffness">Spring constant k in N/m (must be &gt; 0).</param>
        /// <param name="initialPosition">Initial displacement from equilibrium in metres.</param>
        /// <param name="initialVelocity">Initial velocity in m/s.</param>
        public SimpleHarmonicOscillator(double mass, double stiffness,
            double initialPosition = 1.0, double initialVelocity = 0.0)
        {
            if (mass <= 0) throw new ArgumentException("Mass must be greater than zero.", nameof(mass));
            if (stiffness <= 0) throw new ArgumentException("Stiffness must be greater than zero.", nameof(stiffness));

            _mass = mass;
            _stiffness = stiffness;
            _x0 = initialPosition;
            _v0 = initialVelocity;

            _x = _x0;
            _v = _v0;
            _t = 0.0;
        }

        #endregion

        #region IOscillator — State

        /// <inheritdoc />
        public double Position => _x;

        /// <inheritdoc />
        public double Velocity => _v;

        /// <inheritdoc />
        public double Time => _t;

        /// <inheritdoc />
        public double TotalEnergy => KineticEnergy + PotentialEnergy;

        /// <inheritdoc />
        public double KineticEnergy => 0.5 * _mass * _v * _v;

        /// <inheritdoc />
        public double PotentialEnergy => 0.5 * _stiffness * _x * _x;

        #endregion

        #region Physical Properties

        /// <summary>Mass in kg.</summary>
        public double Mass => _mass;

        /// <summary>Spring constant k in N/m.</summary>
        public double Stiffness => _stiffness;

        /// <summary>Natural angular frequency ω₀ = √(k/m) in rad/s.</summary>
        public double AngularFrequency => Math.Sqrt(_stiffness / _mass);

        /// <summary>Natural frequency f₀ = ω₀ / 2π in Hz.</summary>
        public double Frequency => AngularFrequency / (2.0 * Math.PI);

        /// <summary>Period T = 1 / f₀ in seconds.</summary>
        public double Period => 1.0 / Frequency;

        /// <summary>
        /// Amplitude A = √(x₀² + (v₀/ω₀)²), derived from the initial conditions.
        /// </summary>
        public double Amplitude
        {
            get
            {
                var w = AngularFrequency;
                return Math.Sqrt(_x0 * _x0 + (_v0 / w) * (_v0 / w));
            }
        }

        /// <summary>
        /// Phase offset φ such that x(t) = A cos(ω₀t + φ).
        /// Computed from initial conditions: φ = atan2(-v₀/ω₀, x₀).
        /// </summary>
        public double PhaseOffset => Math.Atan2(-_v0 / AngularFrequency, _x0);

        #endregion

        #region Analytic Solution

        /// <summary>
        /// Returns the exact displacement at time <paramref name="t"/>:
        /// x(t) = A cos(ω₀t + φ).
        /// </summary>
        /// <param name="t">Time in seconds.</param>
        public double AnalyticPosition(double t)
        {
            return Amplitude * Math.Cos(AngularFrequency * t + PhaseOffset);
        }

        /// <summary>
        /// Returns the exact velocity at time <paramref name="t"/>:
        /// v(t) = -Aω₀ sin(ω₀t + φ).
        /// </summary>
        /// <param name="t">Time in seconds.</param>
        public double AnalyticVelocity(double t)
        {
            return -Amplitude * AngularFrequency * Math.Sin(AngularFrequency * t + PhaseOffset);
        }

        #endregion

        #region Simulation

        /// <inheritdoc />
        public void Step(double dt)
        {
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            // Velocity Verlet: symplectic, O(dt²), excellent energy conservation.
            // a(t) = -(k/m) * x(t)
            double a = -(_stiffness / _mass) * _x;

            // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
            _x += _v * dt + 0.5 * a * dt * dt;

            // a(t+dt) = -(k/m) * x(t+dt)
            double aNew = -(_stiffness / _mass) * _x;

            // v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
            _v += 0.5 * (a + aNew) * dt;

            _t += dt;
        }

        /// <inheritdoc />
        public void Reset()
        {
            _x = _x0;
            _v = _v0;
            _t = 0.0;
        }

        /// <inheritdoc />
        public List<Serie> Trajectory(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            var result = new List<Serie> { new Serie { Index = _t, Value = _x } };

            while (_t < tEnd)
            {
                Step(dt);
                result.Add(new Serie { Index = _t, Value = _x });
            }

            Reset();
            return result;
        }

        /// <inheritdoc />
        public List<Serie> PhasePortrait(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            var result = new List<Serie> { new Serie { Index = _x, Value = _v } };

            while (_t < tEnd)
            {
                Step(dt);
                result.Add(new Serie { Index = _x, Value = _v });
            }

            Reset();
            return result;
        }

        #endregion

        #region Energy

        /// <summary>
        /// Returns the total energy at each time step from t = 0 to <paramref name="tEnd"/>.
        /// For an undamped SHM this should stay constant.
        /// </summary>
        /// <param name="tEnd">End time in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>Total energy vs time.</returns>
        public List<Serie> EnergyOverTime(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            var result = new List<Serie> { new Serie { Index = _t, Value = TotalEnergy } };

            while (_t < tEnd)
            {
                Step(dt);
                result.Add(new Serie { Index = _t, Value = TotalEnergy });
            }

            Reset();
            return result;
        }

        #endregion

        #region Frequency Spectrum

        /// <summary>
        /// Computes the amplitude spectrum of the oscillator displacement via FFT.
        /// Returns frequency (Hz) vs normalised magnitude.
        /// </summary>
        /// <param name="tEnd">Simulation duration in seconds (determines frequency resolution).</param>
        /// <param name="dt">Time step in seconds (determines Nyquist frequency).</param>
        /// <returns>Frequency bins with magnitudes.</returns>
        public List<Serie> FrequencySpectrum(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            var trajectory = Trajectory(tEnd, dt);
            var samples = trajectory.Select(s => s.Value).ToList();

            // Pad to next power of two for FFT
            int n = 1;
            while (n < samples.Count) n <<= 1;
            while (samples.Count < n) samples.Add(0.0);

            double samplingFrequency = 1.0 / dt;
            var fftResult = samples.FastFourierTransform();

            return fftResult.ToFrequencyResolution(samplingFrequency);
        }

        #endregion

        /// <summary>
        /// Returns a human-readable summary of the oscillator parameters.
        /// </summary>
        public override string ToString()
        {
            return $"SHM: m={_mass:G4} kg, k={_stiffness:G4} N/m, " +
                   $"ω₀={AngularFrequency:G4} rad/s, f={Frequency:G4} Hz, " +
                   $"T={Period:G4} s, A={Amplitude:G4} m";
        }
    }
}
