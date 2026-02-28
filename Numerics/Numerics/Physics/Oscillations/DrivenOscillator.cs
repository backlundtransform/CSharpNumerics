using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Oscillations
{
    /// <summary>
    /// Models a driven (forced) damped harmonic oscillator:
    /// ẍ + 2γẋ + ω₀²x = (F₀/m)cos(ω_d t).
    /// <para>
    /// Extends <see cref="DampedOscillator"/> behaviour with a harmonic driving force.
    /// The system exhibits resonance when the drive frequency approaches the natural
    /// frequency of the system.
    /// </para>
    /// <para>
    /// Integration uses RK4 (4th-order Runge-Kutta) with the time-dependent driving term.
    /// </para>
    /// </summary>
    public class DrivenOscillator : IOscillator
    {
        #region Fields

        private readonly double _mass;
        private readonly double _stiffness;
        private readonly double _damping;
        private readonly double _driveAmplitude; // F₀ (force amplitude in N)
        private readonly double _driveFrequency; // ω_d (angular frequency in rad/s)
        private readonly double _x0;
        private readonly double _v0;

        private double _x;
        private double _v;
        private double _t;

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a driven damped harmonic oscillator.
        /// </summary>
        /// <param name="mass">Mass in kg (must be &gt; 0).</param>
        /// <param name="stiffness">Spring constant k in N/m (must be &gt; 0).</param>
        /// <param name="damping">Damping coefficient c in N·s/m (must be ≥ 0). γ = c / (2m).</param>
        /// <param name="driveAmplitude">Driving force amplitude F₀ in N (must be ≥ 0).</param>
        /// <param name="driveFrequency">Driving angular frequency ω_d in rad/s (must be ≥ 0).</param>
        /// <param name="initialPosition">Initial displacement from equilibrium in metres.</param>
        /// <param name="initialVelocity">Initial velocity in m/s.</param>
        public DrivenOscillator(double mass, double stiffness, double damping,
            double driveAmplitude, double driveFrequency,
            double initialPosition = 0.0, double initialVelocity = 0.0)
        {
            if (mass <= 0) throw new ArgumentException("Mass must be greater than zero.", nameof(mass));
            if (stiffness <= 0) throw new ArgumentException("Stiffness must be greater than zero.", nameof(stiffness));
            if (damping < 0) throw new ArgumentException("Damping coefficient must be non-negative.", nameof(damping));
            if (driveAmplitude < 0) throw new ArgumentException("Drive amplitude must be non-negative.", nameof(driveAmplitude));
            if (driveFrequency < 0) throw new ArgumentException("Drive frequency must be non-negative.", nameof(driveFrequency));

            _mass = mass;
            _stiffness = stiffness;
            _damping = damping;
            _driveAmplitude = driveAmplitude;
            _driveFrequency = driveFrequency;
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

        /// <summary>Driving force amplitude F₀ in N.</summary>
        public double DriveAmplitude => _driveAmplitude;

        /// <summary>Driving angular frequency ω_d in rad/s.</summary>
        public double DriveFrequency => _driveFrequency;

        /// <summary>Damping ratio γ = c / (2m) in rad/s.</summary>
        public double Gamma => _damping / (2.0 * _mass);

        /// <summary>Natural angular frequency ω₀ = √(k/m) in rad/s (undamped, undriven).</summary>
        public double NaturalFrequency => Math.Sqrt(_stiffness / _mass);

        /// <summary>
        /// Damped angular frequency ω_d_free = √(ω₀² − γ²) in rad/s.
        /// The free oscillation frequency (without driving). Returns 0 if overdamped/critically damped.
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
        /// The damping regime of the underlying (undriven) system.
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
        /// Quality factor Q = ω₀ / (2γ).
        /// Returns <see cref="double.PositiveInfinity"/> if γ = 0.
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
        /// Resonance angular frequency ω_r = √(ω₀² − 2γ²).
        /// The drive frequency at which the steady-state amplitude is maximised.
        /// Returns 0 if ω₀² ≤ 2γ² (heavily damped — no resonance peak above DC).
        /// </summary>
        public double ResonanceFrequency
        {
            get
            {
                double w0 = NaturalFrequency;
                double g = Gamma;
                double disc = w0 * w0 - 2.0 * g * g;
                return disc > 0 ? Math.Sqrt(disc) : 0.0;
            }
        }

        /// <summary>
        /// Bandwidth Δω — the full width at half-maximum-power (FWHM) of the resonance curve.
        /// For a lightly damped system: Δω ≈ 2γ.
        /// </summary>
        public double Bandwidth => 2.0 * Gamma;

        #endregion

        #region Steady-State Amplitude & Phase

        /// <summary>
        /// Returns the steady-state displacement amplitude at a given drive frequency ω_d:
        /// A(ω_d) = (F₀/m) / √((ω₀² − ω_d²)² + (2γω_d)²).
        /// </summary>
        /// <param name="wd">Drive angular frequency in rad/s.</param>
        public double SteadyStateAmplitude(double wd)
        {
            double w0sq = _stiffness / _mass;
            double g = Gamma;
            double f0m = _driveAmplitude / _mass;

            double denom2 = (w0sq - wd * wd) * (w0sq - wd * wd) + (2.0 * g * wd) * (2.0 * g * wd);
            return f0m / Math.Sqrt(denom2);
        }

        /// <summary>
        /// Returns the steady-state displacement amplitude at the configured drive frequency.
        /// </summary>
        public double SteadyStateAmplitudeAtDrive => SteadyStateAmplitude(_driveFrequency);

        /// <summary>
        /// Returns the steady-state phase lag φ at a given drive frequency ω_d:
        /// φ(ω_d) = −atan2(2γω_d, ω₀² − ω_d²).
        /// The phase is in the range [−π, 0]: the response always lags the driving force.
        /// </summary>
        /// <param name="wd">Drive angular frequency in rad/s.</param>
        public double SteadyStatePhase(double wd)
        {
            double w0sq = _stiffness / _mass;
            double g = Gamma;
            return -Math.Atan2(2.0 * g * wd, w0sq - wd * wd);
        }

        /// <summary>
        /// Returns the steady-state phase lag at the configured drive frequency.
        /// </summary>
        public double SteadyStatePhaseAtDrive => SteadyStatePhase(_driveFrequency);

        #endregion

        #region Resonance Curve & Phase Response

        /// <summary>
        /// Computes the resonance (amplitude response) curve: steady-state amplitude vs drive frequency.
        /// </summary>
        /// <param name="wMin">Minimum angular frequency in rad/s.</param>
        /// <param name="wMax">Maximum angular frequency in rad/s.</param>
        /// <param name="steps">Number of frequency steps.</param>
        /// <returns>List of (ω_d, A(ω_d)) data points.</returns>
        public List<Serie> ResonanceCurve(double wMin, double wMax, int steps)
        {
            if (steps < 2) throw new ArgumentException("Must have at least 2 steps.", nameof(steps));
            if (wMax <= wMin) throw new ArgumentException("wMax must be greater than wMin.");

            var result = new List<Serie>(steps);
            double dw = (wMax - wMin) / (steps - 1);

            for (int i = 0; i < steps; i++)
            {
                double wd = wMin + i * dw;
                result.Add(new Serie { Index = wd, Value = SteadyStateAmplitude(wd) });
            }

            return result;
        }

        /// <summary>
        /// Computes the phase response curve: steady-state phase lag vs drive frequency.
        /// </summary>
        /// <param name="wMin">Minimum angular frequency in rad/s.</param>
        /// <param name="wMax">Maximum angular frequency in rad/s.</param>
        /// <param name="steps">Number of frequency steps.</param>
        /// <returns>List of (ω_d, φ(ω_d)) data points in radians.</returns>
        public List<Serie> PhaseResponse(double wMin, double wMax, int steps)
        {
            if (steps < 2) throw new ArgumentException("Must have at least 2 steps.", nameof(steps));
            if (wMax <= wMin) throw new ArgumentException("wMax must be greater than wMin.");

            var result = new List<Serie>(steps);
            double dw = (wMax - wMin) / (steps - 1);

            for (int i = 0; i < steps; i++)
            {
                double wd = wMin + i * dw;
                result.Add(new Serie { Index = wd, Value = SteadyStatePhase(wd) });
            }

            return result;
        }

        #endregion

        #region Transfer Function

        /// <summary>
        /// Evaluates the transfer function H(s) = 1 / (s² + 2γs + ω₀²) at a complex frequency s.
        /// This is the Laplace-domain representation of the oscillator's response.
        /// </summary>
        /// <param name="s">Complex Laplace variable s = σ + iω.</param>
        /// <returns>Complex transfer function value H(s).</returns>
        public ComplexNumber TransferFunction(ComplexNumber s)
        {
            double w0sq = _stiffness / _mass;
            double twoGamma = _damping / _mass;

            // s² + 2γs + ω₀²
            ComplexNumber s2 = s * s;
            ComplexNumber denom = s2 + new ComplexNumber(twoGamma, 0) * s + new ComplexNumber(w0sq, 0);

            return new ComplexNumber(1, 0) / denom;
        }

        /// <summary>
        /// Evaluates the frequency response H(iω) at a real frequency ω.
        /// Equivalent to <c>TransferFunction(new ComplexNumber(0, ω))</c>.
        /// </summary>
        /// <param name="w">Angular frequency in rad/s.</param>
        /// <returns>Complex frequency response.</returns>
        public ComplexNumber FrequencyResponse(double w)
        {
            return TransferFunction(new ComplexNumber(0, w));
        }

        #endregion

        #region Simulation

        /// <inheritdoc />
        public void Step(double dt)
        {
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            // RK4 applied to the system:
            //   dx/dt = v
            //   dv/dt = -(k/m)x - (c/m)v + (F₀/m)cos(ω_d t)
            double w0sq = _stiffness / _mass;
            double twoGamma = _damping / _mass;
            double f0m = _driveAmplitude / _mass;
            double wd = _driveFrequency;

            double x = _x;
            double v = _v;
            double t = _t;

            // Stage 1
            double k1x = v;
            double k1v = -w0sq * x - twoGamma * v + f0m * Math.Cos(wd * t);

            // Stage 2
            double t2 = t + 0.5 * dt;
            double x2 = x + 0.5 * dt * k1x;
            double v2 = v + 0.5 * dt * k1v;
            double k2x = v2;
            double k2v = -w0sq * x2 - twoGamma * v2 + f0m * Math.Cos(wd * t2);

            // Stage 3
            double x3 = x + 0.5 * dt * k2x;
            double v3 = v + 0.5 * dt * k2v;
            double k3x = v3;
            double k3v = -w0sq * x3 - twoGamma * v3 + f0m * Math.Cos(wd * t2);

            // Stage 4
            double t4 = t + dt;
            double x4 = x + dt * k3x;
            double v4 = v + dt * k3v;
            double k4x = v4;
            double k4v = -w0sq * x4 - twoGamma * v4 + f0m * Math.Cos(wd * t4);

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
        /// Returns the total mechanical energy (KE + PE) at each time step.
        /// For a driven oscillator, energy can increase as the driving force does work on the system.
        /// In steady-state, the time-averaged energy input equals the time-averaged dissipation.
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
        /// Returns the instantaneous power input from the driving force: P(t) = F(t)·v(t).
        /// The time-averaged power in steady-state equals the energy dissipated by damping.
        /// </summary>
        /// <param name="tEnd">End time in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <returns>Instantaneous power vs time.</returns>
        public List<Serie> PowerInput(double tEnd, double dt)
        {
            if (tEnd <= 0) throw new ArgumentException("End time must be positive.", nameof(tEnd));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            Reset();
            var result = new List<Serie>();
            result.Add(new Serie
            {
                Index = _t,
                Value = _driveAmplitude * Math.Cos(_driveFrequency * _t) * _v
            });

            while (_t < tEnd)
            {
                Step(dt);
                double force = _driveAmplitude * Math.Cos(_driveFrequency * _t);
                result.Add(new Serie { Index = _t, Value = force * _v });
            }

            Reset();
            return result;
        }

        #endregion

        #region Frequency Spectrum

        /// <summary>
        /// Computes the amplitude spectrum of the displacement via FFT.
        /// For a driven oscillator, the dominant peak appears at the driving frequency.
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

        #region Steady-State Detection

        /// <summary>
        /// Simulates the oscillator and returns the trajectory only after transients have decayed,
        /// keeping only the last <paramref name="steadyDuration"/> seconds of a longer simulation.
        /// A good rule of thumb is to discard the first 5/γ seconds (five time constants).
        /// </summary>
        /// <param name="steadyDuration">Duration of steady-state data to collect in seconds.</param>
        /// <param name="dt">Time step in seconds.</param>
        /// <param name="transientDuration">
        /// Duration to simulate before collecting data (seconds). 
        /// If null, defaults to 5/γ (or 50·T₀ if undamped).
        /// </param>
        /// <returns>Displacement vs time (starting from t = transientDuration).</returns>
        public List<Serie> SteadyStateTrajectory(double steadyDuration, double dt,
            double? transientDuration = null)
        {
            if (steadyDuration <= 0) throw new ArgumentException("Duration must be positive.", nameof(steadyDuration));
            if (dt <= 0) throw new ArgumentException("Time step must be positive.", nameof(dt));

            double g = Gamma;
            double transient = transientDuration ??
                (g > 0 ? 5.0 / g : 50.0 * 2.0 * Math.PI / NaturalFrequency);

            Reset();

            // Skip transient
            while (_t < transient)
                Step(dt);

            // Collect steady-state
            double tStart = _t;
            var result = new List<Serie> { new Serie { Index = _t, Value = _x } };

            while (_t - tStart < steadyDuration)
            {
                Step(dt);
                result.Add(new Serie { Index = _t, Value = _x });
            }

            Reset();
            return result;
        }

        /// <summary>
        /// Measures the numerical steady-state amplitude by finding the maximum displacement
        /// after transients have decayed.
        /// </summary>
        /// <param name="dt">Time step in seconds (default: T₀/200 for good resolution).</param>
        /// <returns>Measured steady-state amplitude.</returns>
        public double MeasuredSteadyStateAmplitude(double dt = 0)
        {
            double T0 = 2.0 * Math.PI / NaturalFrequency;
            if (dt <= 0) dt = T0 / 200.0;

            double drivePeriod = _driveFrequency > 0
                ? 2.0 * Math.PI / _driveFrequency
                : T0;

            // Collect several drive periods in steady state
            var steadyData = SteadyStateTrajectory(drivePeriod * 10, dt);
            return steadyData.Max(s => Math.Abs(s.Value));
        }

        #endregion

        /// <summary>
        /// Returns a human-readable summary of the oscillator parameters.
        /// </summary>
        public override string ToString()
        {
            return $"DrivenOscillator: m={_mass:G4} kg, k={_stiffness:G4} N/m, " +
                   $"c={_damping:G4} N·s/m, F₀={_driveAmplitude:G4} N, " +
                   $"ω_d={_driveFrequency:G4} rad/s, ω₀={NaturalFrequency:G4} rad/s, " +
                   $"regime={Regime}, Q={QualityFactor:G4}";
        }
    }
}
