using System;

namespace CSharpNumerics.Physics.Thermodynamics;

/// <summary>
/// One-dimensional transient heat conduction through the thickness of a solid
/// slab — a deck plate, a bulkhead, a wall — with convective exchange on both
/// faces:
/// <para>ρ·c_p·∂T/∂t = k·∂²T/∂x², with −k·∂T/∂x = h·(T_gas − T_surface) at each face.</para>
/// <para>
/// This is the standard fire-engineering treatment of structure: the plate is
/// thin compared to everything around it, so conduction along it is ignored and
/// conduction through it is resolved on its own fine grid. Explicit finite
/// differences with internal sub-stepping, so the caller may advance by any time
/// step and stability is handled here.
/// </para>
/// <para>
/// The slab holds its temperature profile between calls: create one per piece of
/// structure being tracked and call <see cref="Step"/> each simulation step.
/// Units are SI throughout — metres, watts, joules, kelvin.
/// </para>
/// <para>
/// Radiation is not included. Near flames radiation dominates the heating of
/// structure, so results here are a floor on how fast a plate heats, not a
/// ceiling.
/// </para>
/// </summary>
public class HeatSlab
{
    private readonly double[] _temps;
    private readonly double _dx;

    /// <summary>Slab thickness in metres.</summary>
    public double Thickness { get; }

    /// <summary>Thermal conductivity k in W/(m·K).</summary>
    public double Conductivity { get; }

    /// <summary>Density ρ in kg/m³.</summary>
    public double Density { get; }

    /// <summary>Specific heat c_p in J/(kg·K).</summary>
    public double SpecificHeat { get; }

    /// <summary>Number of nodes through the thickness.</summary>
    public int NodeCount { get; }

    /// <summary>Thermal diffusivity α = k / (ρ·c_p) in m²/s.</summary>
    public double ThermalDiffusivity => Conductivity / (Density * SpecificHeat);

    /// <summary>Temperature of the exposed (hot-side) face in kelvin.</summary>
    public double HotSideTemperature => _temps[0];

    /// <summary>Temperature of the unexposed (cold-side) face in kelvin.</summary>
    public double ColdSideTemperature => _temps[NodeCount - 1];

    /// <summary>Mean temperature through the thickness in kelvin.</summary>
    public double AverageTemperature
    {
        get
        {
            double sum = 0;
            for (int i = 0; i < NodeCount; i++) sum += _temps[i];
            return sum / NodeCount;
        }
    }

    /// <summary>
    /// Creates a slab at a uniform initial temperature.
    /// </summary>
    /// <param name="thickness">Thickness in metres. 0.010 is a typical deck plate.</param>
    /// <param name="conductivity">Thermal conductivity in W/(m·K). Steel is about 45.</param>
    /// <param name="density">Density in kg/m³. Steel is about 7850.</param>
    /// <param name="specificHeat">Specific heat in J/(kg·K). Steel is about 490.</param>
    /// <param name="initialTemperature">Uniform starting temperature in kelvin.</param>
    /// <param name="nodeCount">Nodes through the thickness. Five resolves a thin metal plate well.</param>
    public HeatSlab(
        double thickness,
        double conductivity,
        double density,
        double specificHeat,
        double initialTemperature = 293.15,
        int nodeCount = 5)
    {
        if (thickness <= 0) throw new ArgumentException("Thickness must be positive.", nameof(thickness));
        if (conductivity <= 0) throw new ArgumentException("Conductivity must be positive.", nameof(conductivity));
        if (density <= 0) throw new ArgumentException("Density must be positive.", nameof(density));
        if (specificHeat <= 0) throw new ArgumentException("Specific heat must be positive.", nameof(specificHeat));
        if (nodeCount < 2) throw new ArgumentException("At least two nodes are required.", nameof(nodeCount));

        Thickness = thickness;
        Conductivity = conductivity;
        Density = density;
        SpecificHeat = specificHeat;
        NodeCount = nodeCount;
        _dx = thickness / (nodeCount - 1);

        _temps = new double[nodeCount];
        for (int i = 0; i < nodeCount; i++)
            _temps[i] = initialTemperature;
    }

    /// <summary>
    /// Advances the slab by <paramref name="dtSeconds"/>, exchanging heat with
    /// the gas on each side. Internally sub-steps as needed for stability, so
    /// any positive time step is acceptable.
    /// <para>
    /// Returns the heat absorbed from the hot-side gas during the step, in
    /// joules per square metre of plate — positive when the slab is taking heat
    /// out of the gas. Callers coupling a gas model can remove that energy from
    /// the adjacent gas cell to keep the exchange two-way.
    /// </para>
    /// </summary>
    /// <param name="dtSeconds">Time to advance, in seconds.</param>
    /// <param name="hotGasTemperature">Gas temperature at the exposed face, in kelvin.</param>
    /// <param name="hotHeatTransferCoefficient">Convective coefficient at the exposed face in W/(m²·K). Around 25–50 in a fire compartment.</param>
    /// <param name="coldGasTemperature">Gas temperature at the unexposed face, in kelvin.</param>
    /// <param name="coldHeatTransferCoefficient">Convective coefficient at the unexposed face in W/(m²·K). Around 4–10 for still air.</param>
    public double Step(
        double dtSeconds,
        double hotGasTemperature,
        double hotHeatTransferCoefficient,
        double coldGasTemperature,
        double coldHeatTransferCoefficient)
    {
        if (dtSeconds <= 0) throw new ArgumentException("Time step must be positive.", nameof(dtSeconds));
        if (hotHeatTransferCoefficient < 0) throw new ArgumentException("Heat transfer coefficient must be non-negative.", nameof(hotHeatTransferCoefficient));
        if (coldHeatTransferCoefficient < 0) throw new ArgumentException("Heat transfer coefficient must be non-negative.", nameof(coldHeatTransferCoefficient));

        double alpha = ThermalDiffusivity;
        double rhoC = Density * SpecificHeat;

        // Stability. Interior nodes need dt ≤ dx²/(2α); the half-thickness
        // surface nodes additionally feel the convective conductance, which
        // shortens their time constant. 0.4 leaves margin on both.
        double interiorLimit = 0.5 * _dx * _dx / alpha;
        double hMax = Math.Max(hotHeatTransferCoefficient, coldHeatTransferCoefficient);
        double surfaceLimit = rhoC * (_dx / 2) / (hMax + Conductivity / _dx);
        double subLimit = 0.4 * Math.Min(interiorLimit, surfaceLimit);

        int subSteps = Math.Max(1, (int)Math.Ceiling(dtSeconds / subLimit));
        double dt = dtSeconds / subSteps;

        int last = NodeCount - 1;
        var next = new double[NodeCount];
        double absorbed = 0;

        for (int s = 0; s < subSteps; s++)
        {
            // Surface nodes carry half a cell of mass; their energy balance is
            // convection on one side and conduction to the neighbour on the other.
            double hotFlux = hotHeatTransferCoefficient * (hotGasTemperature - _temps[0]);
            next[0] = _temps[0] + dt * (hotFlux + Conductivity * (_temps[1] - _temps[0]) / _dx)
                                     / (rhoC * _dx / 2);

            for (int i = 1; i < last; i++)
                next[i] = _temps[i] + dt * alpha * (_temps[i + 1] - 2 * _temps[i] + _temps[i - 1]) / (_dx * _dx);

            double coldFlux = coldHeatTransferCoefficient * (coldGasTemperature - _temps[last]);
            next[last] = _temps[last] + dt * (coldFlux + Conductivity * (_temps[last - 1] - _temps[last]) / _dx)
                                           / (rhoC * _dx / 2);

            absorbed += hotFlux * dt;

            Array.Copy(next, _temps, NodeCount);
        }

        return absorbed;
    }
}
