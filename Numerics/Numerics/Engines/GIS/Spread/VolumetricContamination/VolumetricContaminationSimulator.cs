using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Engines.GIS.Scenario;
using CSharpNumerics.Engines.GIS.Spread.VolumetricContamination.Enums;
using CSharpNumerics.Engines.GIS.Terrain;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.GIS.Spread.VolumetricContamination;

/// <summary>
/// 3D Eulerian finite-difference advection–diffusion simulator for contaminant
/// transport in volumetric water bodies (lakes, oceans, estuaries).
/// <para>
/// Equation: ∂c/∂t = Dh(∂²c/∂x² + ∂²c/∂y²) + Dv·∂²c/∂z²
///                  − vx·∂c/∂x − vy·∂c/∂y − vz·∂c/∂z − λ·c + S
/// </para>
/// <para>
/// Uses forward Euler time stepping, upwind advection, central-difference diffusion,
/// and first-order exponential decay. Supports anisotropic diffusion (Dh ≠ Dv),
/// a prescribed 3D velocity field, and a land/seabed mask for barriers.
/// </para>
/// </summary>
public class VolumetricContaminationSimulator : ISpreadSimulator
{
    /// <summary>Simulation parameters.</summary>
    public VolumetricContaminationParameters Parameters { get; }

    /// <summary>True if a CFL violation was detected during the last run.</summary>
    public bool CflViolationDetected { get; private set; }

    /// <summary>Creates a 3D volumetric contamination simulator.</summary>
    public VolumetricContaminationSimulator(VolumetricContaminationParameters parameters)
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
        int nz = grid.Nz;
        int nxy = nx * ny;
        int cellCount = nx * ny * nz;
        double dx = grid.Step;
        double dy = grid.Step;
        double dz = grid.Step;
        double dt = timeFrame.StepSeconds;

        double Dh = Parameters.HorizontalDiffusivity;
        double Dv = Parameters.VerticalDiffusivity;
        double lambda = Parameters.DecayRate;
        double toxThreshold = Parameters.ToxicityThresholdMgL;
        bool[] landMask = Parameters.LandMask;

        // Pre-compute velocity field
        var vxArr = new double[cellCount];
        var vyArr = new double[cellCount];
        var vzArr = new double[cellCount];
        if (Parameters.VelocityField != null)
        {
            for (int iz = 0; iz < nz; iz++)
                for (int iy = 0; iy < ny; iy++)
                    for (int ix = 0; ix < nx; ix++)
                    {
                        int idx = iz * nxy + iy * nx + ix;
                        if (landMask != null && landMask[idx]) continue;
                        var (vx, vy, vz) = Parameters.VelocityField(ix, iy, iz);
                        vxArr[idx] = vx;
                        vyArr[idx] = vy;
                        vzArr[idx] = vz;
                    }
        }

        // CFL check
        CflViolationDetected = false;
        for (int i = 0; i < cellCount; i++)
        {
            double speed = Math.Abs(vxArr[i]) + Math.Abs(vyArr[i]) + Math.Abs(vzArr[i]);
            if (speed > 0)
            {
                double minD = Math.Min(dx, Math.Min(dy, dz));
                double courant = speed * dt / minD;
                if (courant > 1.0) { CflViolationDetected = true; break; }
            }
        }

        double diffStabH = Dh * dt / (dx * dx) + Dh * dt / (dy * dy);
        double diffStabV = Dv * dt / (dz * dz);
        if (diffStabH + diffStabV > 0.5)
            CflViolationDetected = true;

        // State arrays
        var concentration = new double[cellCount];
        var exposureTime = new double[cellCount];
        var peakConcentration = new double[cellCount];

        // Source lookup
        var sourceMap = new Dictionary<int, (double concMgL, double durS)>();
        foreach (var (six, siy, siz, conc, dur) in Parameters.Sources)
        {
            if (six >= 0 && six < nx && siy >= 0 && siy < ny && siz >= 0 && siz < nz)
                sourceMap[siz * nxy + siy * nx + six] = (conc, dur);
        }

        var snapshots = new List<SpreadSnapshot>();

        for (int ti = 0; ti < timeFrame.Count; ti++)
        {
            double t = timeFrame.TimeAt(ti);
            var newConc = new double[cellCount];
            Array.Copy(concentration, newConc, cellCount);

            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        int idx = iz * nxy + iy * nx + ix;

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
                            continue;
                        }

                        double c0 = concentration[idx];

                        // Neighbour values (Neumann zero-flux at boundaries and land)
                        double cXp = GetConc(concentration, ix + 1, iy, iz, nx, ny, nz, landMask);
                        double cXm = GetConc(concentration, ix - 1, iy, iz, nx, ny, nz, landMask);
                        double cYp = GetConc(concentration, ix, iy + 1, iz, nx, ny, nz, landMask);
                        double cYm = GetConc(concentration, ix, iy - 1, iz, nx, ny, nz, landMask);
                        double cZp = GetConc(concentration, ix, iy, iz + 1, nx, ny, nz, landMask);
                        double cZm = GetConc(concentration, ix, iy, iz - 1, nx, ny, nz, landMask);

                        // Diffusion: Dh·(∂²c/∂x² + ∂²c/∂y²) + Dv·∂²c/∂z²
                        double diffusion = Dh * (cXp - 2 * c0 + cXm) / (dx * dx)
                                         + Dh * (cYp - 2 * c0 + cYm) / (dy * dy)
                                         + Dv * (cZp - 2 * c0 + cZm) / (dz * dz);

                        // Advection: upwind scheme
                        double vx = vxArr[idx];
                        double vy = vyArr[idx];
                        double vz = vzArr[idx];

                        double advX = vx >= 0
                            ? vx * (c0 - cXm) / dx
                            : vx * (cXp - c0) / dx;

                        double advY = vy >= 0
                            ? vy * (c0 - cYm) / dy
                            : vy * (cYp - c0) / dy;

                        double advZ = vz >= 0
                            ? vz * (c0 - cZm) / dz
                            : vz * (cZp - c0) / dz;

                        newConc[idx] = c0 + dt * (diffusion - advX - advY - advZ);

                        // Decay
                        if (lambda > 0 && newConc[idx] > 0)
                            newConc[idx] *= Math.Exp(-lambda * dt);

                        // Clamp
                        if (newConc[idx] < 0) newConc[idx] = 0;
                    }
                }
            }

            // Update state
            Array.Copy(newConc, concentration, cellCount);

            // Track exposure and peak
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
                velArr[i] = Math.Sqrt(vxArr[i] * vxArr[i] + vyArr[i] * vyArr[i] + vzArr[i] * vzArr[i]);
                expArr[i] = exposureTime[i];

                if (landMask != null && landMask[i])
                    stateArr[i] = (double)ContaminationCellState3D.Land;
                else if (sourceMap.ContainsKey(i) && t < sourceMap[i].durS)
                    stateArr[i] = (double)ContaminationCellState3D.Source;
                else if (concentration[i] > toxThreshold)
                    stateArr[i] = (double)ContaminationCellState3D.Contaminated;
                else
                    stateArr[i] = (double)ContaminationCellState3D.Clean;
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
    /// Gets concentration at (ix, iy, iz) with Neumann zero-flux BC at boundaries
    /// and land cells.
    /// </summary>
    private static double GetConc(double[] c, int ix, int iy, int iz,
        int nx, int ny, int nz, bool[] landMask)
    {
        int cix = Math.Max(0, Math.Min(nx - 1, ix));
        int ciy = Math.Max(0, Math.Min(ny - 1, iy));
        int ciz = Math.Max(0, Math.Min(nz - 1, iz));
        int idx = ciz * (nx * ny) + ciy * nx + cix;

        if (landMask != null && landMask[idx])
            return 0;

        return c[idx];
    }
}
