using System;

namespace CSharpNumerics.Physics.Quantum;

/// <summary>
/// The bound stationary states of a 1-D Schrödinger Hamiltonian: the energy eigenvalues
/// (ascending) and their normalised real wavefunctions sampled on the spatial grid. Produced by
/// <see cref="SchrodingerExtensions.SolveStationaryStates"/>.
/// </summary>
public class StationaryStates
{
    public StationaryStates(double[] energies, double[][] waveFunctions, double[] grid, double gridSpacing)
    {
        Energies = energies ?? throw new ArgumentNullException(nameof(energies));
        WaveFunctions = waveFunctions ?? throw new ArgumentNullException(nameof(waveFunctions));
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        GridSpacing = gridSpacing;
    }

    /// <summary>Energy eigenvalues E₀ ≤ E₁ ≤ … (ground state first).</summary>
    public double[] Energies { get; }

    /// <summary>Normalised wavefunctions; <c>WaveFunctions[k][i]</c> is φₖ at grid point i (∫|φₖ|²dx = 1).</summary>
    public double[][] WaveFunctions { get; }

    /// <summary>Spatial grid points (interior nodes; the wavefunction vanishes at the boundaries).</summary>
    public double[] Grid { get; }

    /// <summary>Uniform grid spacing Δx.</summary>
    public double GridSpacing { get; }

    /// <summary>Number of computed states.</summary>
    public int Count => Energies.Length;

    /// <summary>Ground-state energy E₀.</summary>
    public double GroundStateEnergy => Energies[0];
}
