using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;

/// <summary>
/// Runtime configuration for a 3D volumetric water contamination simulation.
/// Models contaminant transport in open water bodies with depth (lakes, oceans, estuaries).
/// <para>
/// Supports anisotropic diffusion (horizontal ≠ vertical), prescribed 3D current
/// velocity field, first-order decay, and land/seabed masking.
/// </para>
/// </summary>
public class VolumetricContaminationParameters
{
    /// <summary>
    /// Point-source injection sites. Each tuple specifies: cell (ix, iy, iz),
    /// injection concentration (mg/L), and duration of injection (seconds).
    /// </summary>
    public IReadOnlyList<(int ix, int iy, int iz, double concentrationMgL, double durationSeconds)> Sources { get; }

    /// <summary>
    /// Horizontal diffusivity Dh (m²/s). Applied to x and y directions.
    /// Default ~1.0 — typical for lakes.
    /// </summary>
    public double HorizontalDiffusivity { get; }

    /// <summary>
    /// Vertical diffusivity Dv (m²/s). Typically Dv ≪ Dh in stratified water.
    /// Default ~0.01.
    /// </summary>
    public double VerticalDiffusivity { get; }

    /// <summary>
    /// Prescribed 3D velocity field. Returns (vx, vy, vz) in m/s at grid cell (ix, iy, iz).
    /// If null, velocity is zero everywhere (pure diffusion + decay).
    /// </summary>
    public Func<int, int, int, (double vx, double vy, double vz)> VelocityField { get; }

    /// <summary>
    /// First-order decay rate λ (s⁻¹). Zero means conservative tracer.
    /// </summary>
    public double DecayRate { get; }

    /// <summary>
    /// Land/seabed mask. If provided, cells where LandMask[flatIndex] = true
    /// are treated as impermeable barriers (zero concentration, zero flux).
    /// Index order: iz * (Nx*Ny) + iy * Nx + ix.
    /// </summary>
    public bool[] LandMask { get; }

    /// <summary>
    /// Toxicity threshold (mg/L) used to classify cells as contaminated.
    /// </summary>
    public double ToxicityThresholdMgL { get; }

    /// <summary>Creates 3D volumetric contamination simulation parameters.</summary>
    public VolumetricContaminationParameters(
        IReadOnlyList<(int ix, int iy, int iz, double concentrationMgL, double durationSeconds)> sources,
        double horizontalDiffusivity = 1.0,
        double verticalDiffusivity = 0.01,
        Func<int, int, int, (double vx, double vy, double vz)> velocityField = null,
        double decayRate = 0,
        bool[] landMask = null,
        double toxicityThresholdMgL = 0.01)
    {
        if (sources == null || sources.Count == 0)
            throw new ArgumentException("At least one contamination source is required.", nameof(sources));
        if (horizontalDiffusivity < 0)
            throw new ArgumentException("Diffusivity must be non-negative.", nameof(horizontalDiffusivity));
        if (verticalDiffusivity < 0)
            throw new ArgumentException("Diffusivity must be non-negative.", nameof(verticalDiffusivity));

        Sources = sources;
        HorizontalDiffusivity = horizontalDiffusivity;
        VerticalDiffusivity = verticalDiffusivity;
        VelocityField = velocityField;
        DecayRate = decayRate;
        LandMask = landMask;
        ToxicityThresholdMgL = toxicityThresholdMgL;
    }
}
