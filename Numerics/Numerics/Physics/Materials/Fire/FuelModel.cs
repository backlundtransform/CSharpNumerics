using System;
using System.Globalization;
using CSharpNumerics.Physics.Materials.Fire.Enums;

namespace CSharpNumerics.Physics.Materials.Fire;

/// <summary>
/// Immutable descriptor of a Rothermel surface fire fuel model carrying
/// the physical parameters needed for fire spread calculations.
/// <para>
/// Standard Anderson 13 models are available as static fields:
/// <c>FuelModel.ShortGrass</c>, <c>FuelModel.Chaparral</c>, etc.
/// </para>
/// </summary>
public readonly struct FuelModel : IEquatable<FuelModel>
{
    /// <summary>Anderson 13 fuel model type identifier.</summary>
    public FuelModelType Type { get; }

    /// <summary>Descriptive name (e.g. "Short Grass").</summary>
    public string Name { get; }

    /// <summary>Characteristic surface-area-to-volume ratio σ (1/m).</summary>
    public double SurfaceAreaToVolumeRatio { get; }

    /// <summary>Fuel bed depth δ (m).</summary>
    public double FuelBedDepth { get; }

    /// <summary>Ovendry fuel load w₀ (kg/m²).</summary>
    public double OvendryFuelLoad { get; }

    /// <summary>Dead fuel moisture of extinction Mₓ (fraction, 0–1).</summary>
    public double MoistureOfExtinction { get; }

    /// <summary>Low heat content h (kJ/kg).</summary>
    public double LowHeatContent { get; }

    /// <summary>Ovendry particle density ρp (kg/m³). Default 513 for wood.</summary>
    public double ParticleDensity { get; }

    /// <summary>
    /// Creates a fuel model descriptor.
    /// </summary>
    public FuelModel(
        FuelModelType type,
        string name,
        double surfaceAreaToVolumeRatio,
        double fuelBedDepth,
        double ovendryFuelLoad,
        double moistureOfExtinction,
        double lowHeatContent,
        double particleDensity = 513.0)
    {
        Type = type;
        Name = name ?? throw new ArgumentNullException(nameof(name));
        SurfaceAreaToVolumeRatio = surfaceAreaToVolumeRatio;
        FuelBedDepth = fuelBedDepth;
        OvendryFuelLoad = ovendryFuelLoad;
        MoistureOfExtinction = moistureOfExtinction;
        LowHeatContent = lowHeatContent;
        ParticleDensity = particleDensity;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Anderson 13 built-in fuel models
    //  Values from Anderson (1982) "Aids to Determining Fuel Models
    //  For Estimating Fire Behavior", converted to SI where needed.
    //  σ: 1/ft → 1/m  (× 3.28084)
    //  δ: ft → m       (× 0.3048)
    //  w₀: ton/acre → kg/m²  (× 0.224170)
    //  h:  8000 Btu/lb = 18 608 kJ/kg (standard for all models)
    //  Mₓ: fraction
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Model 1 — Short grass (1 ft). Fast-spreading grassland fires.</summary>
    public static readonly FuelModel ShortGrass = new FuelModel(
        FuelModelType.ShortGrass, "Short Grass",
        surfaceAreaToVolumeRatio: 11483,   // 3500 1/ft
        fuelBedDepth: 0.305,               // 1.0 ft
        ovendryFuelLoad: 0.166,            // 0.74 ton/acre
        moistureOfExtinction: 0.12,
        lowHeatContent: 18608);

    /// <summary>Model 2 — Timber (grass and understory).</summary>
    public static readonly FuelModel TimberGrassUnderstory = new FuelModel(
        FuelModelType.TimberGrassUnderstory, "Timber Grass Understory",
        surfaceAreaToVolumeRatio: 9842,    // 3000 1/ft
        fuelBedDepth: 0.305,               // 1.0 ft
        ovendryFuelLoad: 0.896,            // 4.0 ton/acre
        moistureOfExtinction: 0.15,
        lowHeatContent: 18608);

    /// <summary>Model 3 — Tall grass (2.5 ft). Highest spread rate in grass.</summary>
    public static readonly FuelModel TallGrass = new FuelModel(
        FuelModelType.TallGrass, "Tall Grass",
        surfaceAreaToVolumeRatio: 4921,    // 1500 1/ft
        fuelBedDepth: 0.762,               // 2.5 ft
        ovendryFuelLoad: 0.672,            // 3.0 ton/acre
        moistureOfExtinction: 0.25,
        lowHeatContent: 18608);

    /// <summary>Model 4 — Chaparral (6 ft). High intensity brush fires.</summary>
    public static readonly FuelModel Chaparral = new FuelModel(
        FuelModelType.Chaparral, "Chaparral",
        surfaceAreaToVolumeRatio: 6562,    // 2000 1/ft
        fuelBedDepth: 1.829,               // 6.0 ft
        ovendryFuelLoad: 11.230,           // 50.1 ton/acre  (total: 1hr+live)
        moistureOfExtinction: 0.20,
        lowHeatContent: 18608);

    /// <summary>Model 5 — Brush (2 ft). Low brush with some grass.</summary>
    public static readonly FuelModel Brush = new FuelModel(
        FuelModelType.Brush, "Brush",
        surfaceAreaToVolumeRatio: 6562,    // 2000 1/ft
        fuelBedDepth: 0.610,               // 2.0 ft
        ovendryFuelLoad: 0.784,            // 3.5 ton/acre
        moistureOfExtinction: 0.20,
        lowHeatContent: 18608);

    /// <summary>Model 6 — Dormant brush, hardwood slash.</summary>
    public static readonly FuelModel DormantBrush = new FuelModel(
        FuelModelType.DormantBrush, "Dormant Brush",
        surfaceAreaToVolumeRatio: 5741,    // 1750 1/ft
        fuelBedDepth: 0.762,               // 2.5 ft
        ovendryFuelLoad: 2.694,            // 12.0 ton/acre
        moistureOfExtinction: 0.25,
        lowHeatContent: 18608);

    /// <summary>Model 7 — Southern rough. Palmetto-gallberry understory.</summary>
    public static readonly FuelModel SouthernRough = new FuelModel(
        FuelModelType.SouthernRough, "Southern Rough",
        surfaceAreaToVolumeRatio: 5741,    // 1750 1/ft
        fuelBedDepth: 0.762,               // 2.5 ft
        ovendryFuelLoad: 2.471,            // 11.0 ton/acre
        moistureOfExtinction: 0.40,
        lowHeatContent: 18608);

    /// <summary>Model 8 — Closed timber litter. Short-needle compact litter.</summary>
    public static readonly FuelModel ClosedTimberLitter = new FuelModel(
        FuelModelType.ClosedTimberLitter, "Closed Timber Litter",
        surfaceAreaToVolumeRatio: 6562,    // 2000 1/ft
        fuelBedDepth: 0.061,               // 0.2 ft
        ovendryFuelLoad: 2.471,            // 11.0 ton/acre
        moistureOfExtinction: 0.30,
        lowHeatContent: 18608);

    /// <summary>Model 9 — Hardwood litter. Loosely compacted long-needle/hardwood.</summary>
    public static readonly FuelModel HardwoodLitter = new FuelModel(
        FuelModelType.HardwoodLitter, "Hardwood Litter",
        surfaceAreaToVolumeRatio: 8202,    // 2500 1/ft
        fuelBedDepth: 0.061,               // 0.2 ft
        ovendryFuelLoad: 0.651,            // 2.9 ton/acre
        moistureOfExtinction: 0.25,
        lowHeatContent: 18608);

    /// <summary>Model 10 — Timber (litter and understory). Heavy dead/down.</summary>
    public static readonly FuelModel TimberLitterUnderstory = new FuelModel(
        FuelModelType.TimberLitterUnderstory, "Timber Litter Understory",
        surfaceAreaToVolumeRatio: 8202,    // 2500 1/ft
        fuelBedDepth: 0.305,               // 1.0 ft
        ovendryFuelLoad: 2.694,            // 12.0 ton/acre
        moistureOfExtinction: 0.25,
        lowHeatContent: 18608);

    /// <summary>Model 11 — Light logging slash.</summary>
    public static readonly FuelModel LightLoggingSlash = new FuelModel(
        FuelModelType.LightLoggingSlash, "Light Logging Slash",
        surfaceAreaToVolumeRatio: 4921,    // 1500 1/ft
        fuelBedDepth: 0.305,               // 1.0 ft
        ovendryFuelLoad: 2.471,            // 11.0 ton/acre
        moistureOfExtinction: 0.15,
        lowHeatContent: 18608);

    /// <summary>Model 12 — Medium logging slash.</summary>
    public static readonly FuelModel MediumLoggingSlash = new FuelModel(
        FuelModelType.MediumLoggingSlash, "Medium Logging Slash",
        surfaceAreaToVolumeRatio: 4921,    // 1500 1/ft
        fuelBedDepth: 0.701,               // 2.3 ft
        ovendryFuelLoad: 8.085,            // 34.6 ton/acre
        moistureOfExtinction: 0.20,
        lowHeatContent: 18608);

    /// <summary>Model 13 — Heavy logging slash. Highest fuel load.</summary>
    public static readonly FuelModel HeavyLoggingSlash = new FuelModel(
        FuelModelType.HeavyLoggingSlash, "Heavy Logging Slash",
        surfaceAreaToVolumeRatio: 4921,    // 1500 1/ft
        fuelBedDepth: 0.914,               // 3.0 ft
        ovendryFuelLoad: 15.680,           // 69.9 ton/acre
        moistureOfExtinction: 0.25,
        lowHeatContent: 18608);

    // ═══════════════════════════════════════════════════════════════
    //  Equality & display
    // ═══════════════════════════════════════════════════════════════

    public bool Equals(FuelModel other) => Type == other.Type;

    public override bool Equals(object obj) => obj is FuelModel fm && Equals(fm);

    public override int GetHashCode() => Type.GetHashCode();

    public static bool operator ==(FuelModel a, FuelModel b) => a.Equals(b);
    public static bool operator !=(FuelModel a, FuelModel b) => !a.Equals(b);

    public override string ToString() =>
        string.Format(CultureInfo.InvariantCulture,
            "FuelModel {0} — {1} (σ={2:F0} 1/m, δ={3:F3} m, w₀={4:F3} kg/m²)",
            (int)Type, Name, SurfaceAreaToVolumeRatio, FuelBedDepth, OvendryFuelLoad);
}
