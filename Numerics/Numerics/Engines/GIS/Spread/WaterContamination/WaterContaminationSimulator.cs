using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.WaterContamination.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Environmental.Water;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination;

/// <summary>
/// Eulerian finite-difference contaminant transport simulator for river networks.
/// Implements 1-D advection–diffusion with first-order decay, retardation,
/// and confluence mixing along a <see cref="RiverNetwork"/>.
/// <para>
/// Advection uses a first-order upwind scheme. Dispersion uses explicit
/// central differences. Processing order follows topological sort
/// (upstream → downstream) for correct mass balance at confluences.
/// </para>
/// </summary>
public class WaterContaminationSimulator : ISpreadSimulator
{
    /// <summary>Simulation parameters (sources, contaminant, discharge).</summary>
    public WaterContaminationParameters Parameters { get; }

    /// <summary>River network (flow graph).</summary>
    public RiverNetwork Network { get; }

    /// <summary>Channel hydraulic properties.</summary>
    public ChannelMap Channels { get; }

    /// <summary>
    /// If true, a CFL violation was detected during the last run.
    /// </summary>
    public bool CflViolationDetected { get; private set; }

    /// <summary>
    /// Creates a water contamination simulator.
    /// </summary>
    public WaterContaminationSimulator(
        WaterContaminationParameters parameters,
        RiverNetwork network,
        ChannelMap channels)
    {
        Parameters = parameters ?? throw new ArgumentNullException(nameof(parameters));
        Network = network ?? throw new ArgumentNullException(nameof(network));
        Channels = channels ?? throw new ArgumentNullException(nameof(channels));
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
        double dx = grid.Step;
        double dt = timeFrame.StepSeconds;

        var contaminant = Parameters.Contaminant;
        double lambda = contaminant.DecayConstant;
        double Rf = LongitudinalDispersion.RetardationFactor(
            Parameters.BedBulkDensity,
            contaminant.PartitionCoefficient,
            Parameters.BedPorosity);

        // State arrays
        var concentration = new double[cellCount];
        var exposureTime = new double[cellCount];
        var peakConcentration = new double[cellCount];

        // Pre-compute per-cell velocity and dispersion coefficient
        var velocity = new double[cellCount];
        var dispersionEL = new double[cellCount];
        var reachCells = Network.GetReachCells();

        foreach (var (ix, iy) in reachCells)
        {
            int idx = iy * nx + ix;
            double u = Channels.GetVelocity(ix, iy, Network, terrain);
            velocity[idx] = u;

            double w = Channels.GetWidth(ix, iy);
            double d = Channels.GetDepth(ix, iy);
            double slope = Channels.GetBedSlope(ix, iy, Network, terrain);
            double Rh = ManningEquation.RectangularHydraulicRadius(w, d);
            double uStar = LongitudinalDispersion.ShearVelocity(Rh, slope);
            dispersionEL[idx] = LongitudinalDispersion.FischerCoefficient(u, w, d, uStar);
        }

        // CFL stability check
        CflViolationDetected = false;
        foreach (var (ix, iy) in reachCells)
        {
            int idx = iy * nx + ix;
            double courant = velocity[idx] * dt / dx;
            if (courant > 1.0)
            {
                CflViolationDetected = true;
                break;
            }
        }

        // Build source lookup: flat index → (concentrationMgL, durationSeconds)
        var sourceMap = new Dictionary<int, (double concMgL, double durS)>();
        foreach (var (six, siy, conc, dur) in Parameters.Sources)
        {
            if (six >= 0 && six < nx && siy >= 0 && siy < ny)
                sourceMap[siy * nx + six] = (conc, dur);
        }

        // Build upstream index lookup for fast access
        var upstreamIdx = new Dictionary<int, List<int>>();
        foreach (var (ix, iy) in reachCells)
        {
            int idx = iy * nx + ix;
            var ups = Network.GetUpstream(ix, iy);
            if (ups.Count > 0)
            {
                var list = new List<int>(ups.Count);
                foreach (var u in ups)
                    list.Add(u.iy * nx + u.ix);
                upstreamIdx[idx] = list;
            }
        }

        var snapshots = new List<SpreadSnapshot>();

        // Time-step loop
        for (int ti = 0; ti < timeFrame.Count; ti++)
        {
            double t = timeFrame.TimeAt(ti);

            // Work on a new concentration array for this step
            var newConc = new double[cellCount];
            Array.Copy(concentration, newConc, cellCount);

            // Process cells in topological order (upstream → downstream)
            foreach (var (ix, iy) in reachCells)
            {
                int idx = iy * nx + ix;

                // ── Step 1: Source injection ──
                if (sourceMap.TryGetValue(idx, out var src) && t < src.durS)
                {
                    newConc[idx] = src.concMgL;
                }

                // ── Step 2: Advection (upwind scheme) ──
                // C_i^{n+1} = C_i^n - (u·dt)/(Rf·dx) · (C_i^n - C_upstream^n)
                double u = velocity[idx];
                if (u > 0 && Rf > 0)
                {
                    var ups = Network.GetUpstream(ix, iy);
                    double upstreamConc = 0;

                    if (ups.Count == 1)
                    {
                        upstreamConc = concentration[ups[0].iy * nx + ups[0].ix];
                    }
                    else if (ups.Count > 1)
                    {
                        // Confluence mixing: mass-weighted average from upstream branches
                        double totalQ = 0;
                        double totalMass = 0;
                        foreach (var up in ups)
                        {
                            int upIdx = up.iy * nx + up.ix;
                            double upU = velocity[upIdx];
                            double upW = Channels.GetWidth(up.ix, up.iy);
                            double upD = Channels.GetDepth(up.ix, up.iy);
                            double upQ = upU * upW * upD;
                            totalQ += upQ;
                            totalMass += concentration[upIdx] * upQ;
                        }
                        upstreamConc = totalQ > 0 ? totalMass / totalQ : 0;
                    }

                    // Only apply advection to non-active-source cells
                    if (!(sourceMap.TryGetValue(idx, out var srcCheck) && t < srcCheck.durS))
                    {
                        double courant = u * dt / (Rf * dx);
                        newConc[idx] = concentration[idx] - courant * (concentration[idx] - upstreamConc);
                    }
                }

                // ── Step 3: Dispersion (central difference) ──
                double EL = dispersionEL[idx];
                if (EL > 0 && Rf > 0)
                {
                    var ups = Network.GetUpstream(ix, iy);
                    var downs = Network.GetDownstream(ix, iy);

                    double cUp = 0, cDown = 0;
                    if (ups.Count > 0)
                        cUp = concentration[ups[0].iy * nx + ups[0].ix];
                    else
                        cUp = concentration[idx]; // boundary: zero gradient

                    if (downs.Count > 0)
                        cDown = concentration[downs[0].iy * nx + downs[0].ix];
                    else
                        cDown = concentration[idx]; // boundary: zero gradient

                    double dispTerm = EL * dt / (Rf * dx * dx) * (cDown - 2 * concentration[idx] + cUp);

                    // Only apply to non-active-source cells
                    if (!(sourceMap.TryGetValue(idx, out var srcCheck2) && t < srcCheck2.durS))
                    {
                        newConc[idx] += dispTerm;
                    }
                }

                // ── Step 4: Decay ──
                if (lambda > 0 && newConc[idx] > 0)
                {
                    newConc[idx] *= Math.Exp(-lambda * dt);
                }

                // Clamp to non-negative
                if (newConc[idx] < 0) newConc[idx] = 0;
            }

            // Update concentration state
            Array.Copy(newConc, concentration, cellCount);

            // Update exposure time and peak tracking
            double toxThreshold = contaminant.ToxicityThresholdMgL;
            for (int i = 0; i < cellCount; i++)
            {
                if (concentration[i] > toxThreshold)
                    exposureTime[i] += dt;
                if (concentration[i] > peakConcentration[i])
                    peakConcentration[i] = concentration[i];
            }

            // ── Step 6: Record snapshot ──
            var stateArr = new double[cellCount];
            var concArr = new double[cellCount];
            var velArr = new double[cellCount];
            var expArr = new double[cellCount];

            for (int i = 0; i < cellCount; i++)
            {
                concArr[i] = concentration[i];
                velArr[i] = velocity[i];
                expArr[i] = exposureTime[i];

                // Determine state
                if (sourceMap.ContainsKey(i) && t < sourceMap[i].durS)
                    stateArr[i] = (double)CellContaminationState.Source;
                else if (concentration[i] > toxThreshold)
                    stateArr[i] = (double)CellContaminationState.Contaminated;
                else if (peakConcentration[i] > toxThreshold && concentration[i] <= toxThreshold)
                    stateArr[i] = (double)CellContaminationState.Decayed;
                else
                    stateArr[i] = (double)CellContaminationState.Clean;
            }

            var gs = new GridSnapshot(grid, concArr, t, ti);
            gs.SetLayer("concentration", concArr);
            gs.SetLayer("contaminationState", stateArr);
            gs.SetLayer("velocity", velArr);
            gs.SetLayer("exposureTime", expArr);

            snapshots.Add(new SpreadSnapshot(gs));
        }

        return snapshots.AsReadOnly();
    }
}
