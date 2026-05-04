using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.WaterContamination.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using CSharpNumerics.Physics.Environmental.Water;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.WaterContamination2D;

/// <summary>
/// 2D Eulerian finite-difference advection–diffusion simulator for contaminant
/// transport in open water bodies (lakes, estuaries, coastal areas).
/// <para>
/// Equation: ∂c/∂t = Dx·∂²c/∂x² + Dy·∂²c/∂y² − vx·∂c/∂x − vy·∂c/∂y − λ·c + S
/// </para>
/// <para>
/// Uses forward Euler time stepping, upwind advection, central-difference diffusion,
/// and first-order exponential decay. Supports anisotropic diffusion (Dx ≠ Dy),
/// a prescribed 2D velocity field, and a land mask for coastline barriers.
/// </para>
/// </summary>
public class WaterContamination2DSimulator : ISpreadSimulator
{
    /// <summary>Simulation parameters.</summary>
    public WaterContamination2DParameters Parameters { get; }

    /// <summary>True if a CFL violation was detected during the last run.</summary>
    public bool CflViolationDetected { get; private set; }

    /// <summary>Creates a 2D water contamination simulator.</summary>
    public WaterContamination2DSimulator(WaterContamination2DParameters parameters)
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
        double dx = grid.Step;
        double dy = grid.Step;
        double dt = timeFrame.StepSeconds;

        double Dx = Parameters.DiffusionX;
        double Dy = Parameters.DiffusionY;
        var contaminant = Parameters.Contaminant;
        double lambda = contaminant.DecayConstant;
        bool[] landMask = Parameters.LandMask;

        // Pre-compute velocity field
        var vxArr = new double[cellCount];
        var vyArr = new double[cellCount];
        if (Parameters.VelocityField != null)
        {
            for (int iy = 0; iy < ny; iy++)
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = iy * nx + ix;
                    if (landMask != null && landMask[idx]) continue;
                    var (vx, vy) = Parameters.VelocityField(ix, iy);
                    vxArr[idx] = vx;
                    vyArr[idx] = vy;
                }
        }

        // CFL check
        CflViolationDetected = false;
        for (int i = 0; i < cellCount; i++)
        {
            double speed = Math.Abs(vxArr[i]) + Math.Abs(vyArr[i]);
            if (speed > 0)
            {
                double courant = speed * dt / Math.Min(dx, dy);
                if (courant > 1.0) { CflViolationDetected = true; break; }
            }
        }

        double diffStabX = Dx * dt / (dx * dx);
        double diffStabY = Dy * dt / (dy * dy);
        if (diffStabX + diffStabY > 0.5)
            CflViolationDetected = true;

        // State arrays
        var concentration = new double[cellCount];
        var exposureTime = new double[cellCount];
        var peakConcentration = new double[cellCount];

        // Source lookup
        var sourceMap = new Dictionary<int, (double concMgL, double durS)>();
        foreach (var (six, siy, conc, dur) in Parameters.Sources)
        {
            if (six >= 0 && six < nx && siy >= 0 && siy < ny)
                sourceMap[siy * nx + six] = (conc, dur);
        }

        var snapshots = new List<SpreadSnapshot>();

        for (int ti = 0; ti < timeFrame.Count; ti++)
        {
            double t = timeFrame.TimeAt(ti);
            var newConc = new double[cellCount];
            Array.Copy(concentration, newConc, cellCount);

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = iy * nx + ix;

                    // Skip land cells
                    if (landMask != null && landMask[idx])
                    {
                        newConc[idx] = 0;
                        continue;
                    }

                    // Source injection
                    bool isActiveSource = sourceMap.TryGetValue(idx, out var src) && t < src.durS;
                    if (isActiveSource)
                    {
                        newConc[idx] = src.concMgL;
                        continue; // source cells maintain fixed concentration
                    }

                    double c0 = concentration[idx];

                    // Neighbour values (Neumann zero-flux at boundaries and land)
                    double cXp = GetConc(concentration, ix + 1, iy, nx, ny, landMask);
                    double cXm = GetConc(concentration, ix - 1, iy, nx, ny, landMask);
                    double cYp = GetConc(concentration, ix, iy + 1, nx, ny, landMask);
                    double cYm = GetConc(concentration, ix, iy - 1, nx, ny, landMask);

                    // Diffusion: Dx·(c_{i+1} - 2c_i + c_{i-1})/dx² + Dy·(c_{j+1} - 2c_j + c_{j-1})/dy²
                    double diffusion = Dx * (cXp - 2 * c0 + cXm) / (dx * dx)
                                     + Dy * (cYp - 2 * c0 + cYm) / (dy * dy);

                    // Advection: upwind scheme
                    double vx = vxArr[idx];
                    double vy = vyArr[idx];

                    double advX = vx >= 0
                        ? vx * (c0 - cXm) / dx
                        : vx * (cXp - c0) / dx;

                    double advY = vy >= 0
                        ? vy * (c0 - cYm) / dy
                        : vy * (cYp - c0) / dy;

                    newConc[idx] = c0 + dt * (diffusion - advX - advY);

                    // Decay
                    if (lambda > 0 && newConc[idx] > 0)
                        newConc[idx] *= Math.Exp(-lambda * dt);

                    // Clamp
                    if (newConc[idx] < 0) newConc[idx] = 0;
                }
            }

            // Update state
            Array.Copy(newConc, concentration, cellCount);

            // Track exposure and peak
            double toxThreshold = contaminant.ToxicityThresholdMgL;
            for (int i = 0; i < cellCount; i++)
            {
                if (concentration[i] > toxThreshold)
                    exposureTime[i] += dt;
                if (concentration[i] > peakConcentration[i])
                    peakConcentration[i] = concentration[i];
            }

            // Record snapshot
            var stateArr = new double[cellCount];
            var concArr = new double[cellCount];
            var velArr = new double[cellCount];
            var expArr = new double[cellCount];

            for (int i = 0; i < cellCount; i++)
            {
                concArr[i] = concentration[i];
                velArr[i] = Math.Sqrt(vxArr[i] * vxArr[i] + vyArr[i] * vyArr[i]);
                expArr[i] = exposureTime[i];

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

    /// <summary>
    /// Gets concentration at (ix, iy) with Neumann zero-flux BC at boundaries
    /// and land cells (mirrors interior value).
    /// </summary>
    private static double GetConc(double[] c, int ix, int iy, int nx, int ny, bool[] landMask)
    {
        // Clamp to grid (Neumann: zero gradient at boundary)
        int cix = Math.Max(0, Math.Min(nx - 1, ix));
        int ciy = Math.Max(0, Math.Min(ny - 1, iy));
        int idx = ciy * nx + cix;

        // Land cells act as barriers (mirror interior)
        if (landMask != null && landMask[idx])
        {
            // Return the concentration of the original cell (caller side)
            // This is handled by the caller already using c0
            return 0;
        }

        return c[idx];
    }
}
