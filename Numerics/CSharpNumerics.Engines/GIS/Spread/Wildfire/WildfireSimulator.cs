using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.Wildfire.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Environmental.Fire;
using CSharpNumerics.Physics.Materials.Fire;
using CSharpNumerics.Physics.Materials.Fire.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.Wildfire;

/// <summary>
/// Cellular-automaton wildfire spread simulator using the
/// Rothermel (1972) surface fire model.
/// <para>
/// 8-neighbour spread on a <see cref="GeoGrid"/>: for each Burning cell,
/// the rate of spread toward each Unburned neighbour is computed via
/// <see cref="RothermelModel.RateOfSpread"/> with directional slope and
/// wind components. A cell ignites when accumulated travel time exceeds
/// the distance / ROS threshold.
/// </para>
/// </summary>
public class WildfireSimulator : ISpreadSimulator
{
    // 8-neighbour offsets: N, NE, E, SE, S, SW, W, NW
    private static readonly (int dx, int dy)[] _neighbours = new[]
    {
        ( 0,  1), ( 1,  1), ( 1,  0), ( 1, -1),
        ( 0, -1), (-1, -1), (-1,  0), (-1,  1)
    };

    /// <summary>Simulation parameters (ignition, wind, burn duration).</summary>
    public WildfireParameters Parameters { get; }

    /// <summary>
    /// Creates a wildfire simulator with the given parameters.
    /// </summary>
    public WildfireSimulator(WildfireParameters parameters)
    {
        Parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
    }

    /// <inheritdoc/>
    public IReadOnlyList<SpreadSnapshot> Run(
        GeoGrid grid,
        TerrainGrid terrain,
        FuelMap fuelMap,
        TimeFrame timeFrame)
    {
        int nx = grid.Nx;
        int ny = grid.Ny;
        int cellCount = nx * ny;

        // State arrays
        var state = new CellBurnState[cellCount];
        var ignitionTime = new double[cellCount]; // time (s) when cell started burning
        var rosArr = new double[cellCount];       // last computed ROS per cell
        var flameLenArr = new double[cellCount];  // flame length per cell

        // Initialize all to max (unignited)
        for (int i = 0; i < cellCount; i++)
            ignitionTime[i] = double.MaxValue;

        // Apply firebreak cells for NoFuel
        for (int iy = 0; iy < ny; iy++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                if (fuelMap.GetFuelType(ix, iy) == FuelModelType.NoFuel)
                    state[iy * nx + ix] = CellBurnState.Firebreak;
            }
        }

        // Ignite initial cells at t = start
        double tStart = timeFrame.Start;
        foreach (var (ix, iy) in Parameters.IgnitionPoints)
        {
            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) continue;
            int idx = iy * nx + ix;
            if (state[idx] == CellBurnState.Firebreak) continue;
            state[idx] = CellBurnState.Burning;
            ignitionTime[idx] = tStart;

            // Compute initial ROS and flame length at ignition cell
            var fuel = fuelMap.GetFuel(ix, iy);
            if (fuel.OvendryFuelLoad > 0)
            {
                double moisture = fuelMap.GetMoisture(ix, iy);
                double ros = RothermelModel.RateOfSpread(fuel, moisture, Parameters.MidflameWindSpeed, 0);
                rosArr[idx] = ros;
                double ir = RothermelModel.ReactionIntensity(fuel, moisture);
                flameLenArr[idx] = RothermelModel.FlameLength(ir, ros);
            }
        }

        // Wind unit vector (horizontal)
        var windDir = Parameters.WindDirection;
        double windMag = Math.Sqrt(windDir.x * windDir.x + windDir.y * windDir.y);
        Vector windUnit = windMag > 1e-10
            ? new Vector(windDir.x / windMag, windDir.y / windMag, 0)
            : new Vector(0, 0, 0);

        double step = grid.Step;
        var snapshots = new List<SpreadSnapshot>();

        // Time-step loop
        for (int ti = 0; ti < timeFrame.Count; ti++)
        {
            double t = timeFrame.TimeAt(ti);
            double dt = timeFrame.StepSeconds;

            // --- Phase A: transition Burning → Burned ---
            for (int i = 0; i < cellCount; i++)
            {
                if (state[i] == CellBurnState.Burning &&
                    (t - ignitionTime[i]) >= Parameters.BurnDuration)
                {
                    state[i] = CellBurnState.Burned;
                    rosArr[i] = 0;
                    flameLenArr[i] = 0;
                }
            }

            // --- Phase B: spread from Burning cells to Unburned neighbours ---
            // Collect new ignitions to apply after scanning all cells
            var newIgnitions = new List<int>();

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = iy * nx + ix;
                    if (state[idx] != CellBurnState.Burning) continue;

                    double timeBurning = t - ignitionTime[idx];

                    foreach (var (dx, dy) in _neighbours)
                    {
                        int nix = ix + dx;
                        int niy = iy + dy;
                        if (nix < 0 || nix >= nx || niy < 0 || niy >= ny) continue;

                        int nIdx = niy * nx + nix;
                        if (state[nIdx] != CellBurnState.Unburned) continue;

                        // Neighbour fuel
                        var nFuel = fuelMap.GetFuel(nix, niy);
                        if (nFuel.OvendryFuelLoad <= 0) continue; // no fuel

                        double nMoisture = fuelMap.GetMoisture(nix, niy);

                        // Direction vector from burning cell to neighbour
                        var dir = new Vector(dx, dy, 0);
                        double dist = Math.Sqrt(dx * dx + dy * dy) * step;

                        // Slope in spread direction
                        double slopeRad = terrain.SlopeInDirection(ix, iy, dir);
                        // Only uphill slope accelerates fire; negative slope → use 0
                        double effectiveSlope = Math.Max(0, slopeRad);

                        // Wind component along spread direction
                        double dirLen = Math.Sqrt(dx * dx + dy * dy);
                        double dirUnitX = dx / dirLen;
                        double dirUnitY = dy / dirLen;
                        double windComponent = windUnit.x * dirUnitX + windUnit.y * dirUnitY;
                        // Effective wind speed along this direction (0 if headwind)
                        double effectiveWind = Math.Max(0, windComponent * Parameters.MidflameWindSpeed);

                        // Compute ROS toward this neighbour
                        double ros = RothermelModel.RateOfSpread(nFuel, nMoisture, effectiveWind, effectiveSlope);
                        if (ros <= 0) continue;

                        // Travel time = distance / ROS (ROS is m/min, dist is m → time in min)
                        double travelTimeMin = dist / ros;
                        double travelTimeSec = travelTimeMin * 60.0;

                        // Does the fire reach this neighbour?
                        if (timeBurning >= travelTimeSec)
                        {
                            newIgnitions.Add(nIdx);

                            // Store ROS and flame length for the newly ignited cell
                            rosArr[nIdx] = ros;
                            double ir = RothermelModel.ReactionIntensity(nFuel, nMoisture);
                            flameLenArr[nIdx] = RothermelModel.FlameLength(ir, ros);
                        }
                    }
                }
            }

            // Apply new ignitions
            foreach (int nIdx in newIgnitions)
            {
                if (state[nIdx] == CellBurnState.Unburned)
                {
                    state[nIdx] = CellBurnState.Burning;
                    ignitionTime[nIdx] = t;
                }
            }

            // --- Phase C: record snapshot ---
            var bsArr = new double[cellCount];
            var flArr = new double[cellCount];
            var rosOut = new double[cellCount];
            var btArr = new double[cellCount];

            for (int i = 0; i < cellCount; i++)
            {
                bsArr[i] = (double)state[i];
                flArr[i] = flameLenArr[i];
                rosOut[i] = rosArr[i];
                btArr[i] = state[i] == CellBurnState.Burning || state[i] == CellBurnState.Burned
                    ? t - ignitionTime[i]
                    : 0;
            }

            snapshots.Add(new SpreadSnapshot(grid, t, ti, bsArr, flArr, rosOut, btArr));
        }

        return snapshots.AsReadOnly();
    }
}
