using CSharpNumerics.Numerics;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Statistics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CSharpNumerics.Physics.Oscillations
{
    /// <summary>
    /// N-mass spring chain with fixed walls on both ends:
    /// <code>wall —k₀— m₁ —k₁— m₂ —k₂— … —mₙ —kₙ— wall</code>
    /// <para>
    /// Equations of motion: M·ẍ + C·ẋ + K·x = 0 where M is a diagonal
    /// mass matrix, K is the tridiagonal stiffness matrix, and C is a diagonal
    /// damping matrix.
    /// </para>
    /// <para>
    /// Normal modes are computed via the Jacobi eigenvalue algorithm on the
    /// symmetrised dynamical matrix D = L⁻¹KL⁻¹ where L = diag(√m_i),
    /// giving eigenvalues ω² and orthonormal mode shapes.
    /// </para>
    /// <para>
    /// Integration uses RK4 (4th-order Runge-Kutta).
    /// </para>
    /// </summary>
    public class CoupledOscillators
    {
        private readonly double[] _masses;
        private readonly double[] _stiffnesses;
        private readonly double[] _dampings;
        private readonly double[] _x0;
        private readonly double[] _v0;

        private double[] _x;
        private double[] _v;
        private double _t;

        private readonly int _n;

        // Cached eigendecomposition (lazy)
        private double[] _cachedEigenvalues;
        private double[,] _cachedEigenvectors;

        #region Constructors

        /// <summary>
        /// General constructor: each mass and spring can have different values.
        /// </summary>
        /// <param name="masses">Array of N masses (all must be &gt; 0).</param>
        /// <param name="stiffnesses">Array of N+1 spring constants (all must be ≥ 0).</param>
        /// <param name="dampings">Optional array of N damping coefficients (must be ≥ 0).</param>
        /// <param name="initialPositions">Optional initial displacements (length N).</param>
        /// <param name="initialVelocities">Optional initial velocities (length N).</param>
        public CoupledOscillators(
            double[] masses,
            double[] stiffnesses,
            double[] dampings = null,
            double[] initialPositions = null,
            double[] initialVelocities = null)
        {
            if (masses == null || masses.Length == 0)
                throw new ArgumentException("At least one mass is required.", nameof(masses));
            if (stiffnesses == null || stiffnesses.Length != masses.Length + 1)
                throw new ArgumentException(
                    $"Expected {masses.Length + 1} stiffness values for {masses.Length} masses.",
                    nameof(stiffnesses));
            if (masses.Any(m => m <= 0))
                throw new ArgumentException("All masses must be positive.", nameof(masses));
            if (stiffnesses.Any(k => k < 0))
                throw new ArgumentException("All stiffnesses must be non-negative.", nameof(stiffnesses));

            _n = masses.Length;
            _masses = (double[])masses.Clone();
            _stiffnesses = (double[])stiffnesses.Clone();

            if (dampings != null)
            {
                if (dampings.Length != _n)
                    throw new ArgumentException(
                        $"Damping array must have {_n} elements.", nameof(dampings));
                if (dampings.Any(c => c < 0))
                    throw new ArgumentException(
                        "All damping coefficients must be non-negative.", nameof(dampings));
                _dampings = (double[])dampings.Clone();
            }
            else
            {
                _dampings = new double[_n];
            }

            _x0 = initialPositions != null
                ? (double[])initialPositions.Clone()
                : new double[_n];
            _v0 = initialVelocities != null
                ? (double[])initialVelocities.Clone()
                : new double[_n];

            if (_x0.Length != _n)
                throw new ArgumentException(
                    $"Initial positions must have {_n} elements.", nameof(initialPositions));
            if (_v0.Length != _n)
                throw new ArgumentException(
                    $"Initial velocities must have {_n} elements.", nameof(initialVelocities));

            _x = (double[])_x0.Clone();
            _v = (double[])_v0.Clone();
            _t = 0;
        }

        /// <summary>
        /// Uniform chain: all masses equal, all springs equal.
        /// </summary>
        public CoupledOscillators(
            int count,
            double mass,
            double stiffness,
            double damping = 0,
            double[] initialPositions = null,
            double[] initialVelocities = null)
            : this(
                Enumerable.Repeat(mass, count).ToArray(),
                Enumerable.Repeat(stiffness, count + 1).ToArray(),
                Enumerable.Repeat(damping, count).ToArray(),
                initialPositions,
                initialVelocities)
        { }

        #endregion

        #region State Properties

        /// <summary>Number of masses in the chain.</summary>
        public int Count => _n;

        /// <summary>Current simulation time.</summary>
        public double Time => _t;

        /// <summary>Copy of all current positions.</summary>
        public double[] Positions => (double[])_x.Clone();

        /// <summary>Copy of all current velocities.</summary>
        public double[] Velocities => (double[])_v.Clone();

        /// <summary>Position of a specific mass.</summary>
        public double Position(int index)
        {
            if (index < 0 || index >= _n)
                throw new ArgumentOutOfRangeException(nameof(index));
            return _x[index];
        }

        /// <summary>Velocity of a specific mass.</summary>
        public double Velocity(int index)
        {
            if (index < 0 || index >= _n)
                throw new ArgumentOutOfRangeException(nameof(index));
            return _v[index];
        }

        #endregion

        #region Parameter Properties

        /// <summary>Copy of mass values.</summary>
        public double[] Masses => (double[])_masses.Clone();

        /// <summary>Copy of spring stiffness values.</summary>
        public double[] Stiffnesses => (double[])_stiffnesses.Clone();

        /// <summary>Copy of damping coefficients.</summary>
        public double[] Dampings => (double[])_dampings.Clone();

        #endregion

        #region Energy

        /// <summary>Total kinetic energy of all masses.</summary>
        public double KineticEnergy
        {
            get
            {
                double ke = 0;
                for (int i = 0; i < _n; i++)
                    ke += 0.5 * _masses[i] * _v[i] * _v[i];
                return ke;
            }
        }

        /// <summary>
        /// Total potential energy stored in all springs.
        /// Includes wall-to-mass springs and inter-mass springs.
        /// </summary>
        public double PotentialEnergy
        {
            get
            {
                double pe = 0;
                // Left wall to first mass
                pe += 0.5 * _stiffnesses[0] * _x[0] * _x[0];
                // Between adjacent masses
                for (int i = 0; i < _n - 1; i++)
                {
                    double dx = _x[i + 1] - _x[i];
                    pe += 0.5 * _stiffnesses[i + 1] * dx * dx;
                }
                // Last mass to right wall
                pe += 0.5 * _stiffnesses[_n] * _x[_n - 1] * _x[_n - 1];
                return pe;
            }
        }

        /// <summary>Total mechanical energy (kinetic + potential).</summary>
        public double TotalEnergy => KineticEnergy + PotentialEnergy;

        #endregion

        #region Matrices

        /// <summary>
        /// Returns the NxN stiffness matrix K (tridiagonal, symmetric).
        /// K[i,i] = k_i + k_{i+1}; K[i,i±1] = −k.
        /// </summary>
        public Matrix StiffnessMatrix()
        {
            var K = new Matrix(_n, _n);
            for (int i = 0; i < _n; i++)
            {
                K.values[i, i] = _stiffnesses[i] + _stiffnesses[i + 1];
                if (i > 0)
                    K.values[i, i - 1] = -_stiffnesses[i];
                if (i < _n - 1)
                    K.values[i, i + 1] = -_stiffnesses[i + 1];
            }
            return K;
        }

        /// <summary>
        /// Returns the NxN diagonal mass matrix M.
        /// </summary>
        public Matrix MassMatrix()
        {
            var M = new Matrix(_n, _n);
            for (int i = 0; i < _n; i++)
                M.values[i, i] = _masses[i];
            return M;
        }

        #endregion

        #region Eigen-decomposition (Jacobi)

        /// <summary>
        /// Builds the symmetrised dynamical matrix D = L⁻¹ K L⁻¹
        /// where L = diag(√m_i). D is symmetric and its eigenvalues
        /// are ω². The physical mode shapes are φ_i = v_i / √m_i.
        /// </summary>
        private double[,] SymmetricDynamicalMatrix()
        {
            var K = StiffnessMatrix();
            var D = new double[_n, _n];
            for (int i = 0; i < _n; i++)
                for (int j = 0; j < _n; j++)
                    D[i, j] = K.values[i, j] / Math.Sqrt(_masses[i] * _masses[j]);
            return D;
        }

        /// <summary>
        /// Jacobi eigenvalue algorithm for real symmetric matrices.
        /// Returns eigenvalues in ascending order and orthonormal eigenvectors
        /// as columns of the eigenvectors matrix.
        /// </summary>
        private void JacobiEigen(double[,] A, int n, out double[] eigenvalues, out double[,] eigenvectors)
        {
            // Work on a copy
            var S = new double[n, n];
            Array.Copy(A, S, A.Length);

            // Eigenvector accumulator (starts as identity)
            var V = new double[n, n];
            for (int i = 0; i < n; i++)
                V[i, i] = 1.0;

            int maxIterations = 100 * n * n;
            double tol = 1e-12;

            for (int iter = 0; iter < maxIterations; iter++)
            {
                // Find largest off-diagonal element
                int p = 0, q = 1;
                double maxVal = 0;
                for (int i = 0; i < n; i++)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        double absVal = Math.Abs(S[i, j]);
                        if (absVal > maxVal)
                        {
                            maxVal = absVal;
                            p = i;
                            q = j;
                        }
                    }
                }

                if (maxVal < tol) break;

                // Compute rotation angle
                double theta;
                if (Math.Abs(S[p, p] - S[q, q]) < 1e-15)
                {
                    theta = Math.PI / 4.0;
                }
                else
                {
                    theta = 0.5 * Math.Atan2(2.0 * S[p, q], S[p, p] - S[q, q]);
                }

                double c = Math.Cos(theta);
                double s = Math.Sin(theta);

                // Apply Jacobi rotation: S' = G^T S G
                // Update rows/columns p and q
                var Sp = new double[n];
                var Sq = new double[n];
                for (int i = 0; i < n; i++)
                {
                    Sp[i] = c * S[p, i] + s * S[q, i];
                    Sq[i] = -s * S[p, i] + c * S[q, i];
                }
                for (int i = 0; i < n; i++)
                {
                    S[p, i] = Sp[i];
                    S[q, i] = Sq[i];
                    S[i, p] = Sp[i];
                    S[i, q] = Sq[i];
                }
                // Fix the 2x2 block
                double Spp = c * Sp[p] + s * Sp[q];
                double Sqq = -s * Sq[p] + c * Sq[q];
                double Spq = -s * Sp[p] + c * Sp[q];
                S[p, p] = Spp;
                S[q, q] = Sqq;
                S[p, q] = 0;
                S[q, p] = 0;

                // Accumulate eigenvectors: V = V * G
                for (int i = 0; i < n; i++)
                {
                    double vip = V[i, p];
                    double viq = V[i, q];
                    V[i, p] = c * vip + s * viq;
                    V[i, q] = -s * vip + c * viq;
                }
            }

            // Extract eigenvalues and sort by ascending order
            var indices = Enumerable.Range(0, n)
                .OrderBy(i => S[i, i])
                .ToArray();

            eigenvalues = new double[n];
            eigenvectors = new double[n, n];
            for (int k = 0; k < n; k++)
            {
                eigenvalues[k] = S[indices[k], indices[k]];
                for (int i = 0; i < n; i++)
                    eigenvectors[i, k] = V[i, indices[k]];
            }
        }

        /// <summary>
        /// Ensures the eigendecomposition is computed and cached.
        /// </summary>
        private void EnsureEigendecomposition()
        {
            if (_cachedEigenvalues != null) return;

            var D = SymmetricDynamicalMatrix();
            JacobiEigen(D, _n, out var rawEigenvalues, out var rawEigenvectors);

            // rawEigenvectors are for the symmetrised problem D.
            // Physical mode shapes: φ_i = v_i / √m_i
            _cachedEigenvalues = rawEigenvalues;
            _cachedEigenvectors = new double[_n, _n];
            for (int mode = 0; mode < _n; mode++)
            {
                double norm = 0;
                for (int i = 0; i < _n; i++)
                {
                    double val = rawEigenvectors[i, mode] / Math.Sqrt(_masses[i]);
                    _cachedEigenvectors[i, mode] = val;
                    norm += val * val;
                }
                // Normalise (Euclidean)
                norm = Math.Sqrt(norm);
                if (norm > 1e-15)
                {
                    for (int i = 0; i < _n; i++)
                        _cachedEigenvectors[i, mode] /= norm;
                }
            }
        }

        #endregion

        #region Normal Modes

        /// <summary>
        /// Returns the normal mode angular frequencies ω (sorted ascending).
        /// Eigenvalues of M⁻¹K give ω².
        /// </summary>
        public List<double> NormalModes()
        {
            EnsureEigendecomposition();
            return _cachedEigenvalues!
                .Select(ev => Math.Sqrt(Math.Max(ev, 0)))
                .ToList();
        }

        /// <summary>
        /// Returns the mode shapes (eigenvectors), one per mode, sorted by
        /// increasing frequency. Each VectorN has length N.
        /// </summary>
        public List<VectorN> ModeShapes()
        {
            EnsureEigendecomposition();
            var shapes = new List<VectorN>(_n);
            for (int mode = 0; mode < _n; mode++)
            {
                var v = new double[_n];
                for (int i = 0; i < _n; i++)
                    v[i] = _cachedEigenvectors![i, mode];
                shapes.Add(new VectorN(v));
            }
            return shapes;
        }

        #endregion

        #region Dispersion Relation

        /// <summary>
        /// Theoretical dispersion relation for a uniform periodic chain:
        /// ω(k) = 2√(K/m)|sin(ka/2)| where a is the lattice spacing.
        /// Uses the first mass and spring constant as representative values.
        /// </summary>
        public List<Serie> DispersionRelation(double[] kValues, double latticeSpacing = 1.0)
        {
            double kSpring = _stiffnesses[0];
            double m = _masses[0];
            double prefactor = 2.0 * Math.Sqrt(kSpring / m);

            return kValues.Select(kv => new Serie
            {
                Index = kv,
                Value = prefactor * Math.Abs(Math.Sin(kv * latticeSpacing / 2.0))
            }).ToList();
        }

        /// <summary>
        /// Phase velocity v_p = ω/k for the uniform periodic chain.
        /// At k = 0 returns the long-wavelength limit a√(K/m).
        /// </summary>
        public double PhaseVelocity(double k, double latticeSpacing = 1.0)
        {
            double kSpring = _stiffnesses[0];
            double m = _masses[0];

            if (Math.Abs(k) < 1e-15)
                return latticeSpacing * Math.Sqrt(kSpring / m);

            double omega = 2.0 * Math.Sqrt(kSpring / m)
                           * Math.Abs(Math.Sin(k * latticeSpacing / 2.0));
            return omega / Math.Abs(k);
        }

        /// <summary>
        /// Group velocity v_g = dω/dk for the uniform periodic chain,
        /// computed via numerical differentiation.
        /// </summary>
        public double GroupVelocity(double k, double latticeSpacing = 1.0)
        {
            double kSpring = _stiffnesses[0];
            double m = _masses[0];
            double a = latticeSpacing;
            Func<double, double> omega = kv =>
                2.0 * Math.Sqrt(kSpring / m) * Math.Abs(Math.Sin(kv * a / 2.0));
            return omega.Derivate(k);
        }

        #endregion

        #region Integration (RK4)

        /// <summary>
        /// Compute accelerations for a given state.
        /// a_i = (1/m_i)[k_{i-1}(x_{i-1} − x_i) + k_i(x_{i+1} − x_i) − c_i·v_i]
        /// where x_0 = x_{N+1} = 0 (fixed walls).
        /// </summary>
        private double[] Accelerations(double[] x, double[] v)
        {
            var a = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                double xPrev = i > 0 ? x[i - 1] : 0.0;
                double xNext = i < _n - 1 ? x[i + 1] : 0.0;

                double springForce = _stiffnesses[i] * (xPrev - x[i])
                                   + _stiffnesses[i + 1] * (xNext - x[i]);
                double dampingForce = -_dampings[i] * v[i];

                a[i] = (springForce + dampingForce) / _masses[i];
            }
            return a;
        }

        /// <summary>
        /// Advance the system by dt using 4th-order Runge-Kutta.
        /// </summary>
        public void Step(double dt)
        {
            // Stage 1
            var a1 = Accelerations(_x, _v);

            // Stage 2
            var x2 = new double[_n];
            var v2 = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                x2[i] = _x[i] + 0.5 * dt * _v[i];
                v2[i] = _v[i] + 0.5 * dt * a1[i];
            }
            var a2 = Accelerations(x2, v2);

            // Stage 3
            var x3 = new double[_n];
            var v3 = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                x3[i] = _x[i] + 0.5 * dt * v2[i];
                v3[i] = _v[i] + 0.5 * dt * a2[i];
            }
            var a3 = Accelerations(x3, v3);

            // Stage 4
            var x4 = new double[_n];
            var v4 = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                x4[i] = _x[i] + dt * v3[i];
                v4[i] = _v[i] + dt * a3[i];
            }
            var a4 = Accelerations(x4, v4);

            // Combine
            for (int i = 0; i < _n; i++)
            {
                _x[i] += dt / 6.0 * (_v[i] + 2 * v2[i] + 2 * v3[i] + v4[i]);
                _v[i] += dt / 6.0 * (a1[i] + 2 * a2[i] + 2 * a3[i] + a4[i]);
            }
            _t += dt;
        }

        /// <summary>
        /// Reset to initial conditions.
        /// </summary>
        public void Reset()
        {
            _x = (double[])_x0.Clone();
            _v = (double[])_v0.Clone();
            _t = 0;
        }

        #endregion

        #region Trajectory & Phase Portrait

        /// <summary>
        /// Simulates the system and returns the trajectory of every mass.
        /// Outer list: one entry per mass. Inner list: Serie(time, position).
        /// State is reset before and after simulation.
        /// </summary>
        public List<List<Serie>> Trajectory(double tEnd, double dt)
        {
            Reset();
            int steps = (int)(tEnd / dt);

            var result = new List<List<Serie>>(_n);
            for (int i = 0; i < _n; i++)
                result.Add(new List<Serie>(steps + 1));

            // Record initial state
            for (int i = 0; i < _n; i++)
                result[i].Add(new Serie { Index = _t, Value = _x[i] });

            for (int s = 0; s < steps; s++)
            {
                Step(dt);
                for (int i = 0; i < _n; i++)
                    result[i].Add(new Serie { Index = _t, Value = _x[i] });
            }

            Reset();
            return result;
        }

        /// <summary>
        /// Phase portrait for a specific mass: Serie(position, velocity).
        /// State is reset before and after simulation.
        /// </summary>
        public List<Serie> PhasePortrait(int massIndex, double tEnd, double dt)
        {
            if (massIndex < 0 || massIndex >= _n)
                throw new ArgumentOutOfRangeException(nameof(massIndex));

            Reset();
            int steps = (int)(tEnd / dt);
            var result = new List<Serie>(steps + 1);

            result.Add(new Serie { Index = _x[massIndex], Value = _v[massIndex] });
            for (int s = 0; s < steps; s++)
            {
                Step(dt);
                result.Add(new Serie { Index = _x[massIndex], Value = _v[massIndex] });
            }

            Reset();
            return result;
        }

        #endregion

        #region Energy Analysis

        /// <summary>
        /// Total energy over time: Serie(time, totalEnergy).
        /// State is reset before and after simulation.
        /// </summary>
        public List<Serie> EnergyOverTime(double tEnd, double dt)
        {
            Reset();
            int steps = (int)(tEnd / dt);
            var result = new List<Serie>(steps + 1);

            result.Add(new Serie { Index = _t, Value = TotalEnergy });
            for (int s = 0; s < steps; s++)
            {
                Step(dt);
                result.Add(new Serie { Index = _t, Value = TotalEnergy });
            }

            Reset();
            return result;
        }

        /// <summary>
        /// Energy in a specific normal mode over time.
        /// Uses mass-weighted projection onto the mode shape.
        /// E_r = ½ m_r (q̇_r² + ω_r² q_r²) where m_r = φ_r^T M φ_r.
        /// State is reset before and after simulation.
        /// </summary>
        public List<Serie> ModalEnergy(int modeIndex, double tEnd, double dt)
        {
            if (modeIndex < 0 || modeIndex >= _n)
                throw new ArgumentOutOfRangeException(nameof(modeIndex));

            EnsureEigendecomposition();
            var omega = Math.Sqrt(Math.Max(_cachedEigenvalues![modeIndex], 0));

            // Get mode shape (already Euclidean-normalised)
            var phi = new double[_n];
            for (int i = 0; i < _n; i++)
                phi[i] = _cachedEigenvectors![i, modeIndex];

            // Generalised mass: m_r = Σ m_i φ_ri²
            double mNorm = 0;
            for (int i = 0; i < _n; i++)
                mNorm += _masses[i] * phi[i] * phi[i];

            Reset();
            int steps = (int)(tEnd / dt);
            var result = new List<Serie>(steps + 1);

            for (int s = 0; s <= steps; s++)
            {
                double qr = 0, qdotr = 0;
                for (int i = 0; i < _n; i++)
                {
                    qr += _masses[i] * phi[i] * _x[i];
                    qdotr += _masses[i] * phi[i] * _v[i];
                }
                qr /= mNorm;
                qdotr /= mNorm;

                double energy = 0.5 * mNorm * (qdotr * qdotr + omega * omega * qr * qr);
                result.Add(new Serie { Index = _t, Value = energy });

                if (s < steps) Step(dt);
            }

            Reset();
            return result;
        }

        #endregion

        #region Frequency Spectrum

        /// <summary>
        /// FFT amplitude spectrum of a specific mass's trajectory.
        /// State is reset before and after simulation.
        /// </summary>
        public List<Serie> FrequencySpectrum(int massIndex, double tEnd, double dt)
        {
            if (massIndex < 0 || massIndex >= _n)
                throw new ArgumentOutOfRangeException(nameof(massIndex));

            Reset();
            int steps = (int)(tEnd / dt);
            var samples = new List<double>(steps + 1);

            samples.Add(_x[massIndex]);
            for (int s = 0; s < steps; s++)
            {
                Step(dt);
                samples.Add(_x[massIndex]);
            }

            // Pad to next power of 2
            int nPad = 1;
            while (nPad < samples.Count) nPad <<= 1;
            while (samples.Count < nPad) samples.Add(0.0);

            double samplingFrequency = 1.0 / dt;
            var fftResult = samples.FastFourierTransform();

            Reset();
            return fftResult.ToFrequencyResolution(samplingFrequency);
        }

        #endregion

        #region ToString

        /// <inheritdoc/>
        public override string ToString()
        {
            string masses = string.Join(", ", _masses.Select(m => $"{m:G4}"));
            string springs = string.Join(", ", _stiffnesses.Select(k => $"{k:G4}"));
            bool hasDamping = _dampings.Any(c => c > 0);

            return $"CoupledOscillators(N={_n}, " +
                   $"m=[{masses}] kg, " +
                   $"k=[{springs}] N/m" +
                   (hasDamping
                       ? $", c=[{string.Join(", ", _dampings.Select(c => $"{c:G4}"))}] Ns/m"
                       : "") +
                   ")";
        }

        #endregion
    }
}
