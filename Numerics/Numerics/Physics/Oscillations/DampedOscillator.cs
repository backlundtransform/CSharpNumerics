using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Oscillations
{
    /// <summary>
    /// Models a damped harmonic oscillator: ẍ + 2γẋ + ω₀²x = 0.
    /// <para>
    /// The damping coefficient γ determines the regime:
    /// <list type="bullet">
    ///   <item><b>Underdamped</b> (γ &lt; ω₀): oscillates with exponentially decaying amplitude.</item>
    ///   <item><b>Critically damped</b> (γ = ω₀): returns to equilibrium as fast as possible without oscillating.</item>
    ///   <item><b>Overdamped</b> (γ &gt; ω₀): decays exponentially without oscillation, slower than critical.</item>
    /// </list>
    /// </para>
    /// <para>
    /// Integration uses RK4 (4th-order Runge-Kutta), which is appropriate for dissipative
    /// systems where symplectic (energy-conserving) integrators are not required.
    /// </para>
    /// </summary>
    public class DampedOscillator : IOscillator
    {
        #region Fields

        private readonly double _mass;
        private readonly double _stiffness;
        private readonly double _damping;
        private readonly double _x0;
        private readonly double _v0;

        private double _x;
        private double _v;
        private double _t;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a damped harmonic oscillator.
        /// </summary>
        /// <param name="mass">Mass in kg (must be &gt; 0).</param>
        /// <param name="stiffness">Spring constant k in N/m (must be &gt; 0).</param>
        /// <param name="damping">Damping coefficient c in N·s/m (must be ≥ 0). γ = c / (2m).</param>
        /// <param name="initialPosition">Initial displacement from equilibrium in metres.</param>
        /// <param name="initialVelocity">Initial velocity in m/s.</param>
        public DampedOscillator(double mass, double stiffness, double damping,
            double initialPosition = 1.0, double initialVelocity = 0.0)
        {
            if (mass <= 0) throw new ArgumentException("Mass must be greater than zero.", nameof(mass));
            if (stiffness <= 0) throw new ArgumentException("Stiffness must be greater than zero.", nameof(stiffness));
            if (damping < 0) throw new ArgumentException("Damping coefficient must be non-negative.", nameof(damping));

            _mass = mass;
            _stiffness = stiffness;
            _damping = damping;
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

        /// <summary>Damping coefficient c in N·s/m.</summary>
        public double Damping => _damping;

        /// <summary>
        /// Damping ratio γ = c / (2m) in rad/s.
        /// </summary>
        public double Gamma => _damping / (2.0 * _mass);

        /// <summary>Natural angular frequency ω₀ = √(k/m) in rad/s (undamped).</summary>
        public double NaturalFrequency => Math.Sqrt(_stiffness / _mass);

        /// <summary>
        /// Damped angular frequency ω_d = √(ω₀² − γ²) in rad/s.
        /// Only meaningful for underdamped systems; returns 0 for critically/overdamped.
        /// </summary>
        public double DampedFrequency
        {
            get
            {
                double w0 = NaturalFrequency;
                double g = Gamma;
                double disc = w0 * w0 - g * g;
                return disc > 0 ? Math.Sqrt(disc) : 0.0;
            }
        }

        /// <summary>
        /// The current damping regime based on the relationship between γ and ω₀.
        /// Uses a relative tolerance of 1e-10 for the critical-damping boundary.
        /// </summary>
        public DampingRegime Regime
        {
            get
            {
                double g = Gamma;
                double w0 = NaturalFrequency;
                double ratio = g / w0;

                if (Math.Abs(ratio - 1.0) < 1e-10)
                    return DampingRegime.CriticallyDamped;
                return g < w0 ? DampingRegime.Underdamped : DampingRegime.Overdamped;
            }
        }

        /// <summary>
        /// Quality factor Q = ω₀ / (2γ). Higher Q means less energy loss per cycle.
        /// Returns <see cref="double.PositiveInfinity"/> if γ = 0 (undamped).
        /// </summary>
        public double QualityFactor
        {
            get
            {
                double g = Gamma;
                return g == 0 ? double.PositiveInfinity : NaturalFrequency / (2.0 * g);
            }
        }

        /// <summary>
        /// Logarithmic decrement δ = γ · T_d, where T_d = 2π/ω_d.
        /// Measures the rate of amplitude decay per cycle.
        /// Only meaningful for underdamped systems; returns <see cref="double.NaN"/> otherwise.
        /// </summary>
        public double LogarithmicDecrement
        {
            get
            {
                if (Regime != DampingRegime.Underdamped) return double.NaN;
                double wd = DampedFrequency;
                double Td = 2.0 * Math.PI / wd;
                return Gamma * Td;
            }
        }

        /// <summary>
        /// Damped period T_d = 2π/ω_d in seconds.
        /// Only meaningful for underdamped systems; returns <see cref="double.NaN"/> otherwise.
        /// </summary>
        public double DampedPeriod
        {
            get
            {
                if (Regime != DampingRegime.Underdamped) return double.NaN;
                return 2.0 * Math.PI / DampedFrequency;
            }
        }

        #endregion

        #region Analytic Solution

        /// <summary>
        /// Returns the exact displacement at time <paramref name="t"/> for all three damping regimes.
        /// <list type="bullet">
        ///   <item>Underdamped: x(t) = A·e^(−γt)·cos(ω_d·t + φ)</item>
        ///   <item>Critically damped: x(t) = (C₁ + C₂·t)·e^(−γt)</item>
        ///   <item>Overdamped: x(t) = C₁·e^(r₁t) + C₂·e^(r₂t)</item>
        /// </list>
        /// </summary>
        /// <param name="t">Time in seconds.</param>
        public double AnalyticPosition(double t)
        {
            double g = Gamma;
            double w0 = NaturalFrequency;

            switch (Regime)
            {
                case DampingRegime.Underdamped:
                {
                    double wd = DampedFrequency;
                    // A and φ from initial conditions: x₀ = A cos(φ), v₀ = -γA cos(φ) + ωdA(-sin(φ))
                    double A = Math.Sqrt(_x0 * _x0 + ((_v0 + g * _x0) / wd) * ((_v0 + g * _x0) / wd));
                    double phi = Math.Atan2(-(_v0 + g * _x0) / wd, _x0);
                    return A * Math.Exp(-g * t) * Math.Cos(wd * t + phi);
                }

                case DampingRegime.CriticallyDamped:
                {
                    // x(t) = (C₁ + C₂t) e^(-γt), C₁ = x₀, C₂ = v₀ + γx₀
                    double C1 = _x0;
                    double C2 = _v0 + g * _x0;
                    return (C1 + C2 * t) * Math.Exp(-g * t);
                }

                case DampingRegime.Overdamped:
                {
                    // Roots: r₁,₂ = -γ ± √(γ² − ω₀²)
                    double s = Math.Sqrt(g * g - w0 * w0);
                    double r1 = -g + s;
                    double r2 = -g - s;
                    // C₁ and C₂ from: x₀ = C₁ + C₂, v₀ = C₁r₁ + C₂r₂
                    double C1 = (_v0 - r2 * _x0) / (r1 - r2);
                    double C2 = (r1 * _x0 - _v0) / (r1 - r2);
                    return C1 * Math.Exp(r1 * t) + C2 * Math.Exp(r2 * t);
                }

                default:
                    throw new InvalidOperationException("Unknown damping regime.");
            }
        }

        /// <summary>
        /// Returns the exact velocity at time <paramref name="t"/> (time derivative of the analytic position).
        /// </summary>
        /// <param name="t">Time in seconds.</param>
        public double AnalyticVelocity(double t)
        {
            double g = Gamma;
            double w0 = NaturalFrequency;

            switch (Regime)
            {
                case DampingRegime.Underdamped:
                {
                    double wd = DampedFrequency;
                    double A = Math.Sqrt(_x0 * _x0 + ((_v0 + g * _x0) / wd) * ((_v0 + g * _x0) / wd));
                    double phi = Math.Atan2(-(_v0 + g * _x0) / wd, _x0);
                    // v(t) = d/dt [A e^(-γt) cos(ωd t + φ)]
                    //       = A e^(-γt) [-γ cos(ωd t + φ) - ωd sin(ωd t + φ)]
                    return A * Math.Exp(-g * t) * (-g * Math.Cos(wd * t + phi) - wd * Math.Sin(wd * t + phi));
                }

                case DampingRegime.CriticallyDamped:
                {
                    double C1 = _x0;
                    double C2 = _v0 + g * _x0;
                    // v(t) = d/dt [(C₁ + C₂t) e^(-γt)] = (C₂ - γ(C₁ + C₂t)) e^(-γt)
                    return (C2 - g * (C1 + C2 * t)) * Math.Exp(-g * t);
                }

                case DampingRegime.Overdamped:
                {
                    double s = Math.Sqrt(g * g - w0 * w0);
                    double r1 = -g + s;
                    double r2 = -g - s;
                    double C1 = (_v0 - r2 * _x0) / (r1 - r2);
                    double C2 = (r1 * _x0 - _v0) / (r1 - r2);
                    return C1 * r1 * Math.Exp(r1 * t) + C2 * r2 * Math.Exp(r2 * t);
                }

                default:
                    throw new InvalidOperationException("Unknown damping regime.");
            }
        }

        #endregion

        #region Envelope

        /// <summary>
        /// Returns the exponential envelope A·e^(−γt) at time <paramref name="t"/>.
        /// For underdamped systems, the displacement stays within ±envelope.
        /// </summary>
        /// <param name="t">Time in seconds.</param>
        public double Envelope(double t)
        {
            double g = Gamma;
            double wd = DampedFrequency;
            if (wd == 0)
            {
                // For critically/overdamped, envelope concept is less meaningful;
                // return |analytic position| as a fallback.
                return Math.Abs(AnalyticPosition(t));
            }

            double A = Math.Sqrt(_x0 * _x0 + ((_v0 + g * _x0) / wd) * ((_v0 + g * _x0) / wd));
            return A * Math.Exp(-g * t);
        }

        /// <summary>
        /// Returns the envelope curve from t = 0 to <paramref name="tEnd"/>.
        /// </summary>
        /// <param name="tEnd">End time in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>Envelope amplitude vs time.</returns>
        public List<Serie> EnvelopeCurve(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            var result = new List<Serie>();
            for (double t = 0; t <= tEnd; t += dt)
            {
                result.Add(new Serie { Index = t, Value = Envelope(t) });
            }
            return result;
        }

        #endregion

        #region Simulation

        /// <inheritdoc />
        public void Step(double dt)
        {
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            // RK4 applied to the system:
            //   dx/dt = v
            //   dv/dt = -(k/m)x - (c/m)v = -ω₀²x - 2γv
            double w0sq = _stiffness / _mass;
            double twoGamma = _damping / _mass; // 2γ = c/m

            double x = _x;
            double v = _v;

            // Stage 1
            double k1x = v;
            double k1v = -w0sq * x - twoGamma * v;

            // Stage 2
            double x2 = x + 0.5 * dt * k1x;
            double v2 = v + 0.5 * dt * k1v;
            double k2x = v2;
            double k2v = -w0sq * x2 - twoGamma * v2;

            // Stage 3
            double x3 = x + 0.5 * dt * k2x;
            double v3 = v + 0.5 * dt * k2v;
            double k3x = v3;
            double k3v = -w0sq * x3 - twoGamma * v3;

            // Stage 4
            double x4 = x + dt * k3x;
            double v4 = v + dt * k3v;
            double k4x = v4;
            double k4v = -w0sq * x4 - twoGamma * v4;

            _x += dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x);
            _v += dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v);
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
        /// For a damped oscillator the total mechanical energy decreases monotonically.
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

        /// <summary>
        /// Returns the energy dissipated by damping from t = 0 to <paramref name="tEnd"/>.
        /// Dissipated energy = E(0) − E(t) at each step.
        /// </summary>
        /// <param name="tEnd">End time in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>Dissipated energy vs time.</returns>
        public List<Serie> EnergyDissipation(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            double E0 = TotalEnergy;
            var result = new List<Serie> { new Serie { Index = _t, Value = 0.0 } };

            while (_t < tEnd)
            {
                Step(dt);
                result.Add(new Serie { Index = _t, Value = E0 - TotalEnergy });
            }

            Reset();
            return result;
        }

        #endregion

        #region Frequency Spectrum

        /// <summary>
        /// Computes the amplitude spectrum of the displacement via FFT.
        /// For underdamped oscillators, the peak is at ω_d with a spectral width proportional to γ.
        /// </summary>
        /// <param name="tEnd">Simulation duration in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>Frequency (Hz) vs normalised magnitude.</returns>
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
            return $"DampedOscillator: m={_mass:G4} kg, k={_stiffness:G4} N/m, " +
                   $"c={_damping:G4} N·s/m, γ={Gamma:G4} rad/s, " +
                   $"ω₀={NaturalFrequency:G4} rad/s, regime={Regime}, Q={QualityFactor:G4}";
        }
    }
}
