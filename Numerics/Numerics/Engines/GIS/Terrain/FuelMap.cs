using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Physics.Materials.Fire;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Terrain;

/// <summary>
/// Per-cell fuel model and moisture assignment on a <see cref="GeoGrid"/>.
/// Each (ix, iy) cell is assigned a <see cref="FuelModel"/> and a dead-fuel
/// moisture content used by the Rothermel spread model.
/// </summary>
public class FuelMap
{
    private readonly FuelModelType[] _fuelTypes;
    private readonly double[] _moisture;

    /// <summary>The underlying spatial grid.</summary>
    public GeoGrid Grid { get; }

    /// <summary>Number of cells along x.</summary>
    public int Nx => Grid.Nx;

    /// <summary>Number of cells along y.</summary>
    public int Ny => Grid.Ny;

    /// <summary>Default moisture content for new cells (fraction, 0–1).</summary>
    public const double DefaultMoisture = 0.08;

    /// <summary>
    /// Creates a fuel map for the given grid. All cells default to
    /// <see cref="FuelModelType.NoFuel"/> with <see cref="DefaultMoisture"/>.
    /// </summary>
    public FuelMap(GeoGrid grid)
    {
        Grid = grid ?? throw new ArgumentNullException(nameof(grid));
        int count = grid.Nx * grid.Ny;
        _fuelTypes = new FuelModelType[count];
        _moisture = new double[count];
        for (int i = 0; i < count; i++)
            _moisture[i] = DefaultMoisture;
    }

    // ═══════════════════════════════════════════════════════════════
    //  Fuel assignment
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Sets the fuel model type for cell (ix, iy).</summary>
    public void SetFuel(int ix, int iy, FuelModelType type) =>
        _fuelTypes[iy * Nx + ix] = type;

    /// <summary>Sets all cells to the same fuel model type.</summary>
    public void SetUniformFuel(FuelModelType type)
    {
        for (int i = 0; i < _fuelTypes.Length; i++)
            _fuelTypes[i] = type;
    }

    /// <summary>
    /// Assigns fuel models based on elevation bands. Each range maps
    /// an elevation interval [minElevation, maxElevation) to a fuel type.
    /// Requires a <see cref="TerrainGrid"/> for elevation lookup.
    /// </summary>
    /// <param name="terrain">Terrain grid providing elevation values.</param>
    /// <param name="ranges">
    /// List of (minElevation, maxElevation, fuelType) tuples.
    /// Ranges are evaluated in order; last match wins.
    /// </param>
    public void SetFuelByElevation(TerrainGrid terrain, IReadOnlyList<(double minElevation, double maxElevation, FuelModelType type)> ranges)
    {
        if (terrain == null) throw new ArgumentNullException(nameof(terrain));
        if (ranges == null) throw new ArgumentNullException(nameof(ranges));

        for (int iy = 0; iy < Ny; iy++)
        {
            for (int ix = 0; ix < Nx; ix++)
            {
                double elev = terrain[ix, iy];
                for (int r = ranges.Count - 1; r >= 0; r--)
                {
                    var (min, max, type) = ranges[r];
                    if (elev >= min && elev < max)
                    {
                        SetFuel(ix, iy, type);
                        break;
                    }
                }
            }
        }
    }

    /// <summary>
    /// Returns the <see cref="FuelModel"/> for cell (ix, iy).
    /// Looks up the fuel type in <see cref="FuelLibrary"/>.
    /// </summary>
    public FuelModel GetFuel(int ix, int iy)
    {
        var type = _fuelTypes[iy * Nx + ix];
        if (type == FuelModelType.NoFuel)
            return default;
        return FuelLibrary.Get(type);
    }

    /// <summary>Returns the <see cref="FuelModelType"/> for cell (ix, iy).</summary>
    public FuelModelType GetFuelType(int ix, int iy) =>
        _fuelTypes[iy * Nx + ix];

    // ═══════════════════════════════════════════════════════════════
    //  Moisture
    // ═══════════════════════════════════════════════════════════════

    /// <summary>Gets the dead fuel moisture content (fraction) for cell (ix, iy).</summary>
    public double GetMoisture(int ix, int iy) =>
        _moisture[iy * Nx + ix];

    /// <summary>Sets the dead fuel moisture content (fraction) for cell (ix, iy).</summary>
    public void SetMoisture(int ix, int iy, double moisture) =>
        _moisture[iy * Nx + ix] = moisture;

    /// <summary>Sets all cells to the same moisture content.</summary>
    public void SetUniformMoisture(double moisture)
    {
        for (int i = 0; i < _moisture.Length; i++)
            _moisture[i] = moisture;
    }
}
