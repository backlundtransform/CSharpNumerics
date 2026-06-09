using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// Extension methods for the one-dimensional Schrödinger equation.
///
/// Time-independent (TISE):   −(ħ²/2m) d²ψ/dx² + V(x)ψ = Eψ
/// Time-dependent  (TDSE):    iħ ∂Ψ/∂t = −(ħ²/2m) ∂²Ψ/∂x² + V(x)Ψ
///
/// The stationary states are found by discretising the Hamiltonian on a uniform grid with a
/// three-point Laplacian (Dirichlet walls) and diagonalising the resulting symmetric matrix.
/// Defaults use natural units (ħ = 1, m = 1); pass explicit constants for SI calculations.
/// </summary>
public static class SchrodingerExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Time-independent Schrödinger equation (eigenvalue problem)
    // ═══════════════════════════════════════════════════════════════

    #region Stationary states

    /// <summary>
    /// Solves the time-independent Schrödinger equation for a potential <c>V(x)</c> on
    /// <c>[xMin, xMax]</c> using a finite-difference Hamiltonian with Dirichlet boundaries
    /// (ψ = 0 at the walls). Returns the lowest <paramref name="states"/> energy eigenvalues and
    /// their normalised wavefunctions.
    /// </summary>
    /// <param name="potential">Potential energy function V(x).</param>
    /// <param name="xMin">Left domain boundary.</param>
    /// <param name="xMax">Right domain boundary.</param>
    /// <param name="points">Number of interior grid points (matrix dimension).</param>
    /// <param name="mass">Particle mass (default 1).</param>
    /// <param name="hbar">Reduced Planck constant (default 1, natural units).</param>
    /// <param name="states">Number of low-lying states to return (default: all).</param>
    public static StationaryStates SolveStationaryStates(
        this Func<double, double> potential,
        double xMin, double xMax, int points,
        double mass = 1.0, double hbar = 1.0, int states = -1)
    {
        if (potential == null) throw new ArgumentNullException(nameof(potential));
        if (points < 3) throw new ArgumentOutOfRangeException(nameof(points), "At least 3 grid points are required.");
        if (xMax <= xMin) throw new ArgumentException("xMax must be greater than xMin.");
        if (mass <= 0) throw new ArgumentOutOfRangeException(nameof(mass));
        if (hbar <= 0) throw new ArgumentOutOfRangeException(nameof(hbar));

        double dx = (xMax - xMin) / (points + 1);
        var grid = new double[points];
        for (int i = 0; i < points; i++)
            grid[i] = xMin + ((i + 1) * dx);

        // Kinetic prefactor: −ħ²/2m · d²/dx² with the 3-point stencil [−1, 2, −1]/dx².
        double kinetic = (hbar * hbar) / (2.0 * mass * dx * dx);

        var h = new double[points, points];
        for (int i = 0; i < points; i++)
        {
            h[i, i] = (2.0 * kinetic) + potential(grid[i]);
            if (i > 0) h[i, i - 1] = -kinetic;
            if (i < points - 1) h[i, i + 1] = -kinetic;
        }

        var (energies, vectors) = SymmetricEigenSolver.Solve(h);

        int count = states < 0 ? points : Math.Min(states, points);
        var selectedEnergies = new double[count];
        var waveFunctions = new double[count][];

        for (int k = 0; k < count; k++)
        {
            selectedEnergies[k] = energies[k];
            waveFunctions[k] = NormaliseReal(vectors[k], dx);
        }

        return new StationaryStates(selectedEnergies, waveFunctions, grid, dx);
    }

    private static double[] NormaliseReal(double[] vector, double dx)
    {
        double norm = 0.0;
        foreach (double v in vector) norm += v * v;
        norm = Math.Sqrt(norm * dx);

        // Fix the global sign so the largest-magnitude lobe is positive (deterministic output).
        int peak = 0;
        for (int i = 1; i < vector.Length; i++)
            if (Math.Abs(vector[i]) > Math.Abs(vector[peak])) peak = i;
        double sign = vector[peak] < 0 ? -1.0 : 1.0;

        var result = new double[vector.Length];
        for (int i = 0; i < vector.Length; i++)
            result[i] = sign * vector[i] / norm;
        return result;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Analytic energy levels
    // ═══════════════════════════════════════════════════════════════

    #region Analytic energies

    /// <summary>
    /// Energy of the n-th level of an infinite square well of width L:
    /// <c>Eₙ = n²π²ħ² / (2mL²)</c>, n = 1, 2, 3, …
    /// </summary>
    public static double InfiniteSquareWellEnergy(this int n, double width, double mass = 1.0, double hbar = 1.0)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), "Quantum number must be ≥ 1.");
        if (width <= 0) throw new ArgumentOutOfRangeException(nameof(width));
        return (n * n * Math.PI * Math.PI * hbar * hbar) / (2.0 * mass * width * width);
    }

    /// <summary>
    /// Energy of the n-th level of a quantum harmonic oscillator:
    /// <c>Eₙ = ħω(n + ½)</c>, n = 0, 1, 2, …
    /// </summary>
    public static double HarmonicOscillatorEnergy(this int n, double angularFrequency, double hbar = 1.0)
    {
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n), "Quantum number must be ≥ 0.");
        return hbar * angularFrequency * (n + 0.5);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Wave–particle relations
    // ═══════════════════════════════════════════════════════════════

    #region de Broglie & dispersion

    /// <summary>
    /// de Broglie wavelength of a particle with momentum p: <c>λ = h / p</c> (SI units, J·s / (kg·m/s)).
    /// </summary>
    public static double DeBroglieWavelength(this double momentum)
    {
        if (momentum <= 0) throw new ArgumentOutOfRangeException(nameof(momentum), "Momentum must be positive.");
        return PhysicsConstants.PlancksConstant / momentum;
    }

    /// <summary>
    /// de Broglie wavelength from non-relativistic kinetic energy E and mass m:
    /// <c>λ = h / √(2mE)</c>.
    /// </summary>
    public static double DeBroglieWavelengthFromEnergy(this double kineticEnergy, double mass)
    {
        if (kineticEnergy <= 0) throw new ArgumentOutOfRangeException(nameof(kineticEnergy));
        if (mass <= 0) throw new ArgumentOutOfRangeException(nameof(mass));
        return PhysicsConstants.PlancksConstant / Math.Sqrt(2.0 * mass * kineticEnergy);
    }

    /// <summary>
    /// Energy of a photon (or quantum) of angular frequency ω: <c>E = ħω</c>.
    /// </summary>
    public static double QuantumEnergy(this double angularFrequency, double hbar = 1.0)
        => hbar * angularFrequency;

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Time evolution (spectral method)
    // ═══════════════════════════════════════════════════════════════

    #region Time evolution

    /// <summary>
    /// Evolves a superposition of stationary states to time t using the exact spectral solution of
    /// the TDSE: <c>Ψ(x, t) = Σₙ cₙ φₙ(x) e^(−iEₙt/ħ)</c>. The coefficients weight each eigenstate;
    /// if they satisfy <c>Σ|cₙ|² = 1</c> the result is normalised and stays so for all t.
    /// </summary>
    /// <param name="states">The stationary states (basis).</param>
    /// <param name="coefficients">Expansion coefficients cₙ (one per state, real-valued).</param>
    /// <param name="time">Evolution time t.</param>
    /// <param name="hbar">Reduced Planck constant (default 1).</param>
    /// <returns>The complex wavefunction Ψ(x, t) sampled on the grid.</returns>
    public static ComplexNumber[] Evolve(
        this StationaryStates states, double[] coefficients, double time, double hbar = 1.0)
    {
        if (states == null) throw new ArgumentNullException(nameof(states));
        if (coefficients == null) throw new ArgumentNullException(nameof(coefficients));
        if (coefficients.Length != states.Count)
            throw new ArgumentException("One coefficient is required per stationary state.", nameof(coefficients));

        int gridSize = states.Grid.Length;
        var psi = new ComplexNumber[gridSize];
        for (int i = 0; i < gridSize; i++)
            psi[i] = new ComplexNumber(0.0, 0.0);

        for (int k = 0; k < states.Count; k++)
        {
            double phase = -states.Energies[k] * time / hbar;
            var rotation = new ComplexNumber(Math.Cos(phase), Math.Sin(phase));
            double[] phi = states.WaveFunctions[k];
            double c = coefficients[k];

            for (int i = 0; i < gridSize; i++)
                psi[i] += (c * phi[i]) * rotation;
        }

        return psi;
    }

    #endregion
}
