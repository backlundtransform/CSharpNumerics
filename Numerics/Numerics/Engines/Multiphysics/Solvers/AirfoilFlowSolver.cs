using CSharpNumerics.Engines.Multiphysics.Enums;
using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;
using System;

namespace CSharpNumerics.Engines.Multiphysics.Solvers;

/// <summary>
/// Solves 2D inviscid potential flow past an airfoil using the Hess-Smith panel method.
/// Computes pressure coefficient distribution, lift coefficient, moment coefficient,
/// and optionally the full velocity field around the airfoil.
///
/// Physics:
///   - Inviscid, incompressible, irrotational flow (potential flow).
///   - Kutta condition enforced at the trailing edge.
///   - Airfoil geometry from NACA 4-digit series.
///
/// Outputs:
///   - CpDistribution: surface pressure coefficient.
///   - LiftCoefficient: Cl from panel method.
///   - MomentCoefficient: Cm about quarter-chord.
///   - Vx, Vy: optional 2D velocity field on a grid around the airfoil.
/// </summary>
internal class AirfoilFlowSolver : IMultiphysicsSolver
{
    public SimulationResult Solve(SimulationBuilder cfg)
    {
        var mat = cfg.Material.Value;
        double alpha = cfg.AngleOfAttack;
        double uInf = cfg.Freestream;

        // ── Generate airfoil geometry ────────────────────────────
        var (xAirfoil, yAirfoil) = NACAGeometry.Generate(
            cfg.AirfoilNACA, cfg.AirfoilPanels, cfg.AirfoilChord);

        // Rotate for angle of attack
        var (xRot, yRot) = NACAGeometry.Rotate(xAirfoil, yAirfoil, alpha, cfg.AirfoilChord);

        // ── Solve panel method ───────────────────────────────────
        var panelResult = PanelMethod.Solve(xRot, yRot, alpha, uInf);

        // ── Build velocity field grid if geometry is specified ────
        double[,] vxField = null;
        double[,] vyField = null;
        double[,] pressureField = null;

        if (cfg.Nx > 0 && cfg.Ny > 0)
        {
            double width = cfg.GeomWidth > 0 ? cfg.GeomWidth : 4.0 * cfg.AirfoilChord;
            double height = cfg.GeomHeight > 0 ? cfg.GeomHeight : 4.0 * cfg.AirfoilChord;
            int nx = cfg.Nx;
            int ny = cfg.Ny;

            double xMin = -0.5 * cfg.AirfoilChord;
            double yMin = -height / 2.0;
            double dx = width / (nx - 1);
            double dy = height / (ny - 1);

            var queryX = new double[nx * ny];
            var queryY = new double[nx * ny];

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = iy * nx + ix;
                    queryX[idx] = xMin + ix * dx;
                    queryY[idx] = yMin + iy * dy;
                }
            }

            var (uArr, vArr) = PanelMethod.VelocityField(
                panelResult, xRot, yRot, queryX, queryY, uInf);

            vxField = new double[nx, ny];
            vyField = new double[nx, ny];
            pressureField = new double[nx, ny];

            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    int idx = iy * nx + ix;
                    vxField[ix, iy] = uArr[idx];
                    vyField[ix, iy] = vArr[idx];

                    // Cp = 1 - (V/U∞)²
                    double vMag2 = uArr[idx] * uArr[idx] + vArr[idx] * vArr[idx];
                    pressureField[ix, iy] = 0.5 * mat.Density * (uInf * uInf - vMag2);
                }
            }
        }

        // ── Compute dynamic pressure and dimensional lift ────────
        double qInf = 0.5 * mat.Density * uInf * uInf;
        double liftPerSpan = panelResult.Cl * qInf * cfg.AirfoilChord;

        return new SimulationResult
        {
            Type = MultiphysicsType.AirfoilFlow,
            CpDistribution = panelResult.Cp,
            SurfaceVelocity = panelResult.Vt,
            AirfoilX = panelResult.Xm,
            AirfoilY = panelResult.Ym,
            LiftCoefficient = panelResult.Cl,
            MomentCoefficient = panelResult.Cm,
            Circulation = panelResult.Gamma,
            DragCoefficient = 0.0, // inviscid → zero pressure drag for thin airfoil
            Vx = vxField,
            Vy = vyField,
            Pressure = pressureField,
            MaxValue = panelResult.Cl,
            MinValue = ArrayMin(panelResult.Cp),
            Iterations = 1
        };
    }

    private static double ArrayMin(double[] arr)
    {
        double min = double.MaxValue;
        for (int i = 0; i < arr.Length; i++)
            if (arr[i] < min) min = arr[i];
        return min;
    }
}
