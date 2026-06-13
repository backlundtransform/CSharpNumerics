using CSharpNumerics.Physics.Materials.Water;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination2D;

/// <summary>
/// Runtime configuration for a 2D water contamination spread simulation.
/// Unlike the 1D version, this operates on the full 2D grid with anisotropic
/// diffusion and a prescribed 2D velocity field — not restricted to a river network.
/// </summary>
public class WaterContamination2DParameters
{
    /// <summary>
    /// Point-source injection sites. Each tuple specifies: cell (ix, iy),
    /// injection concentration (mg/L), and duration of injection (seconds).
    /// After the duration expires the source stops injecting.
    /// </summary>
    public IReadOnlyList<(int ix, int iy, double concentrationMgL, double durationSeconds)> Sources { get; }

    /// <summary>
    /// Contaminant descriptor (decay half-life, partition coefficient, toxicity thresholds).
    /// </summary>
    public AquaticContaminant Contaminant { get; }

    /// <summary>
    /// Horizontal diffusion coefficient Dx (m²/s) in x-direction.
    /// Default 1.0 — typical for lakes.
    /// </summary>
    public double DiffusionX { get; }

    /// <summary>
    /// Horizontal diffusion coefficient Dy (m²/s) in y-direction.
    /// Default equals <see cref="DiffusionX"/> (isotropic).
    /// </summary>
    public double DiffusionY { get; }

    /// <summary>
    /// Prescribed 2D velocity field. Returns (vx, vy) in m/s at grid cell (ix, iy).
    /// If null, velocity is zero everywhere (pure diffusion + decay).
    /// </summary>
    public Func<int, int, (double vx, double vy)> VelocityField { get; }

    /// <summary>
    /// Land mask. If provided, cells where LandMask[iy * Nx + ix] = true are treated
    /// as impermeable barriers (zero concentration, zero flux).
    /// </summary>
    public bool[] LandMask { get; }

    /// <summary>
    /// Creates 2D water contamination simulation parameters.
    /// </summary>
    public WaterContamination2DParameters(
        IReadOnlyList<(int ix, int iy, double concentrationMgL, double durationSeconds)> sources,
        AquaticContaminant contaminant,
        double diffusionX = 1.0,
        double diffusionY = -1,
        Func<int, int, (double vx, double vy)> velocityField = null,
        bool[] landMask = null)
    {
        if (sources == null || sources.Count == 0)
            throw new ArgumentException("At least one contamination source is required.", nameof(sources));

        Sources = sources;
        Contaminant = contaminant;
        DiffusionX = diffusionX;
        DiffusionY = diffusionY < 0 ? diffusionX : diffusionY;
        VelocityField = velocityField;
        LandMask = landMask;
    }
}
