using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Terrain;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread;

/// <summary>
/// Interface for any terrain-aware spread simulation (wildfire, flood, etc.).
/// </summary>
public interface ISpreadSimulator
{
    /// <summary>
    /// Runs the spread simulation and returns one <see cref="SpreadSnapshot"/>
    /// per time step.
    /// </summary>
    IReadOnlyList<SpreadSnapshot> Run(
        GeoGrid grid,
        TerrainGrid terrain,
        FuelMap fuelMap,
        TimeFrame timeFrame);
}
