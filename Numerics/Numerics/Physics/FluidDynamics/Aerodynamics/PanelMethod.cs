using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Result of a panel method solution for 2D potential flow around a body.
/// </summary>
public class PanelMethodResult
{
    /// <summary>Pressure coefficient at each panel midpoint.</summary>
    public double[] Cp { get; }

    /// <summary>Tangential velocity at each panel midpoint (normalised by freestream).</summary>
    public double[] Vt { get; }

    /// <summary>Panel midpoint x-coordinates.</summary>
    public double[] Xm { get; }

    /// <summary>Panel midpoint y-coordinates.</summary>
    public double[] Ym { get; }

    /// <summary>Lift coefficient (integrated from Cp).</summary>
    public double Cl { get; }

    /// <summary>Moment coefficient about the quarter-chord (positive nose-up).</summary>
    public double Cm { get; }

    /// <summary>Source strengths on each panel.</summary>
    public double[] Sigma { get; }

    /// <summary>Vortex strength (uniform circulation Γ per unit span).</summary>
    public double Gamma { get; }

    /// <summary>Angle of attack used in the solution (radians).</summary>
    public double Alpha { get; }

    /// <summary>Number of panels.</summary>
    public int NumPanels { get; }

    internal PanelMethodResult(
        double[] cp, double[] vt, double[] xm, double[] ym,
        double cl, double cm, double[] sigma, double gamma, double alpha, int numPanels)
    {
        Cp = cp;
        Vt = vt;
        Xm = xm;
        Ym = ym;
        Cl = cl;
        Cm = cm;
        Sigma = sigma;
        Gamma = gamma;
        Alpha = alpha;
        NumPanels = numPanels;
    }
}

/// <summary>
/// Hess-Smith panel method for 2D incompressible potential flow around arbitrary closed bodies.
/// Combines constant-strength source panels with a uniform vortex distribution to enforce
/// the Kutta condition at the trailing edge.
///
/// Algorithm:
///   1. Discretize body into N flat panels with constant source strength σᵢ and uniform vortex γ.
///   2. Apply N normal-flow boundary conditions (no penetration) + 1 Kutta condition.
///   3. Solve (N+1) × (N+1) linear system for σᵢ and γ.
///   4. Compute tangential velocity and pressure coefficient on each panel.
/// </summary>
public static class PanelMethod
{
    /// <summary>
    /// Solves potential flow around an airfoil defined by (x, y) coordinates.
    /// </summary>
    /// <param name="x">X-coordinates of panel endpoints (length N+1, closed contour).</param>
    /// <param name="y">Y-coordinates of panel endpoints (length N+1, closed contour).</param>
    /// <param name="alpha">Angle of attack in radians.</param>
    /// <param name="freestream">Freestream velocity magnitude (default 1.0).</param>
    /// <returns>Panel method result with Cp, Cl, velocities.</returns>
    public static PanelMethodResult Solve(double[] x, double[] y, double alpha, double freestream = 1.0)
    {
        if (x == null || y == null)
            throw new ArgumentNullException("Airfoil coordinates cannot be null.");
        if (x.Length != y.Length || x.Length < 4)
            throw new ArgumentException("Need at least 3 panels (4 endpoints).");

        int n = x.Length - 1; // number of panels

        // ── Panel geometry ───────────────────────────────────────
        var xm = new double[n];   // midpoints
        var ym = new double[n];
        var sLen = new double[n]; // panel lengths
        var sinT = new double[n]; // sin(theta_i) - panel angle
        var cosT = new double[n]; // cos(theta_i)

        for (int i = 0; i < n; i++)
        {
            xm[i] = 0.5 * (x[i] + x[i + 1]);
            ym[i] = 0.5 * (y[i] + y[i + 1]);
            double dx = x[i + 1] - x[i];
            double dy = y[i + 1] - y[i];
            sLen[i] = Math.Sqrt(dx * dx + dy * dy);
            cosT[i] = dx / sLen[i];
            sinT[i] = dy / sLen[i];
        }

        // ── Freestream components ────────────────────────────────
        double uInf = freestream * Math.Cos(alpha);
        double vInf = freestream * Math.Sin(alpha);

        // ── Influence coefficients ───────────────────────────────
        // A[i,j]: normal velocity at panel i due to unit source on panel j
        // B[i,j]: normal velocity at panel i due to unit vortex on panel j
        var A = new double[n + 1, n + 1];   // (N+1) system: N panels + Kutta
        var rhs = new double[n + 1];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    // Self-influence: source panel normal = 0.5, vortex = 0
                    A[i, j] = 0.5;
                    continue;
                }

                var (sn, st, vn, vt2) = PanelInfluence(
                    xm[i], ym[i], x[j], y[j], x[j + 1], y[j + 1],
                    sinT[i], cosT[i], sinT[j], cosT[j]);

                A[i, j] = sn;       // source normal influence
                A[i, n] += vn;      // vortex normal influence (accumulated in last column)
            }

            // RHS: -U∞ · n̂ᵢ
            rhs[i] = -(uInf * sinT[i] - vInf * cosT[i]);
        }

        // ── Kutta condition ──────────────────────────────────────
        // Tangential velocity on first + last panel = 0 (equal and opposite → zero net circulation at TE)
        // V_t(panel 0) + V_t(panel N-1) = 0
        for (int j = 0; j < n; j++)
        {
            if (j != 0)
            {
                var (_, st0, _, _) = PanelInfluence(
                    xm[0], ym[0], x[j], y[j], x[j + 1], y[j + 1],
                    sinT[0], cosT[0], sinT[j], cosT[j]);
                A[n, j] += st0;
            }
            if (j != n - 1)
            {
                var (_, stN, _, _) = PanelInfluence(
                    xm[n - 1], ym[n - 1], x[j], y[j], x[j + 1], y[j + 1],
                    sinT[n - 1], cosT[n - 1], sinT[j], cosT[j]);
                A[n, j] += stN;
            }
        }

        // Vortex tangential contribution to Kutta row
        double vtVortex0 = 0;
        double vtVortexN = 0;
        for (int j = 0; j < n; j++)
        {
            if (j != 0)
            {
                var (_, _, _, vt0) = PanelInfluence(
                    xm[0], ym[0], x[j], y[j], x[j + 1], y[j + 1],
                    sinT[0], cosT[0], sinT[j], cosT[j]);
                vtVortex0 += vt0;
            }
            if (j != n - 1)
            {
                var (_, _, _, vtN) = PanelInfluence(
                    xm[n - 1], ym[n - 1], x[j], y[j], x[j + 1], y[j + 1],
                    sinT[n - 1], cosT[n - 1], sinT[j], cosT[j]);
                vtVortexN += vtN;
            }
        }
        A[n, n] = vtVortex0 + vtVortexN;

        // Self tangential of panel 0 and panel N-1 (source self-tangential = 0 by symmetry)
        // RHS Kutta: negative of freestream tangential on both panels
        double utFS0 = uInf * cosT[0] + vInf * sinT[0];
        double utFSN = uInf * cosT[n - 1] + vInf * sinT[n - 1];
        rhs[n] = -(utFS0 + utFSN);

        // ── Solve linear system ──────────────────────────────────
        var solution = SolveLinearSystem(A, rhs, n + 1);
        var sigma = new double[n];
        Array.Copy(solution, sigma, n);
        double gamma = solution[n];

        // ── Compute tangential velocity & Cp ─────────────────────
        var vtArr = new double[n];
        var cp = new double[n];

        for (int i = 0; i < n; i++)
        {
            double vTang = uInf * cosT[i] + vInf * sinT[i];

            for (int j = 0; j < n; j++)
            {
                if (i == j) continue;

                var (_, stij, _, vtij) = PanelInfluence(
                    xm[i], ym[i], x[j], y[j], x[j + 1], y[j + 1],
                    sinT[i], cosT[i], sinT[j], cosT[j]);

                vTang += sigma[j] * stij + gamma * vtij;
            }

            vtArr[i] = vTang;
            cp[i] = 1.0 - (vTang * vTang) / (freestream * freestream);
        }

        // ── Integrate for Cl and Cm ─────────────────────────────
        double cl = 0;
        double cm = 0;
        double chord = MaxChord(x);

        for (int i = 0; i < n; i++)
        {
            // Force per panel: -Cp * panel_length * normal direction
            double fx = -cp[i] * sLen[i] * sinT[i];   // x-force (normal = (-sin, cos))
            double fy = cp[i] * sLen[i] * cosT[i];    // y-force

            // Lift = force perpendicular to freestream
            cl += -fx * Math.Sin(alpha) + fy * Math.Cos(alpha);

            // Moment about quarter-chord
            double xRef = 0.25 * chord;
            cm += cp[i] * sLen[i] * ((xm[i] - xRef) * cosT[i] + ym[i] * sinT[i]);
        }

        cl /= chord;
        cm /= (chord * chord);

        return new PanelMethodResult(cp, vtArr, xm, ym, cl, cm, sigma, gamma, alpha, n);
    }

    /// <summary>
    /// Computes the velocity field at arbitrary points in the flow using the panel method solution.
    /// </summary>
    /// <param name="result">A previously computed panel method result.</param>
    /// <param name="x">Panel endpoint x-coordinates used in the solution.</param>
    /// <param name="y">Panel endpoint y-coordinates used in the solution.</param>
    /// <param name="queryX">X-coordinates of query points.</param>
    /// <param name="queryY">Y-coordinates of query points.</param>
    /// <param name="freestream">Freestream velocity magnitude.</param>
    /// <returns>Velocity components (u, v) at each query point.</returns>
    public static (double[] u, double[] v) VelocityField(
        PanelMethodResult result, double[] x, double[] y,
        double[] queryX, double[] queryY, double freestream = 1.0)
    {
        int nq = queryX.Length;
        int n = result.NumPanels;
        var u = new double[nq];
        var v = new double[nq];

        double uInf = freestream * Math.Cos(result.Alpha);
        double vInf = freestream * Math.Sin(result.Alpha);

        var sinT = new double[n];
        var cosT = new double[n];
        for (int j = 0; j < n; j++)
        {
            double dx = x[j + 1] - x[j];
            double dy = y[j + 1] - y[j];
            double sLen = Math.Sqrt(dx * dx + dy * dy);
            cosT[j] = dx / sLen;
            sinT[j] = dy / sLen;
        }

        for (int q = 0; q < nq; q++)
        {
            double ux = uInf;
            double uy = vInf;

            for (int j = 0; j < n; j++)
            {
                var (uSrc, vSrc, uVort, vVort) = PanelVelocity(
                    queryX[q], queryY[q], x[j], y[j], x[j + 1], y[j + 1]);

                ux += result.Sigma[j] * uSrc + result.Gamma * uVort;
                uy += result.Sigma[j] * vSrc + result.Gamma * vVort;
            }

            u[q] = ux;
            v[q] = uy;
        }

        return (u, v);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Private helpers
    // ═══════════════════════════════════════════════════════════════

    /// <summary>
    /// Computes influence coefficients of panel j on control point i.
    /// Returns (source_normal, source_tangential, vortex_normal, vortex_tangential).
    /// </summary>
    private static (double sn, double st, double vn, double vt) PanelInfluence(
        double xci, double yci, double xj1, double yj1, double xj2, double yj2,
        double sinTi, double cosTi, double sinTj, double cosTj)
    {
        // Transform control point to panel j local coordinates
        double dx1 = xci - xj1;
        double dy1 = yci - yj1;
        double dx2 = xci - xj2;
        double dy2 = yci - yj2;

        // Panel length
        double sj = Math.Sqrt((xj2 - xj1) * (xj2 - xj1) + (yj2 - yj1) * (yj2 - yj1));

        // Local coordinates (ξ, η) in panel j frame
        double xi = dx1 * cosTj + dy1 * sinTj;
        double eta = -dx1 * sinTj + dy1 * cosTj;

        // Geometric quantities
        double r1sq = dx1 * dx1 + dy1 * dy1;
        double r2sq = dx2 * dx2 + dy2 * dy2;

        double r1 = Math.Sqrt(r1sq);
        double r2 = Math.Sqrt(r2sq);

        // Avoid log(0)
        double logTerm = (r1 < 1e-14 || r2 < 1e-14) ? 0.0 : Math.Log(r1sq / r2sq);
        double thetaTerm = Math.Atan2(eta * sj, eta * eta + xi * xi - xi * sj);

        // Velocity in panel j local frame due to unit source
        double uSource = logTerm / (4.0 * Math.PI);
        double vSource = thetaTerm / (2.0 * Math.PI);

        // Velocity due to unit vortex (rotated source)
        double uVortex = thetaTerm / (2.0 * Math.PI);
        double vVortex = -logTerm / (4.0 * Math.PI);

        // Transform to global frame
        double uSrcGlobal = uSource * cosTj - vSource * sinTj;
        double vSrcGlobal = uSource * sinTj + vSource * cosTj;
        double uVrtGlobal = uVortex * cosTj - vVortex * sinTj;
        double vVrtGlobal = uVortex * sinTj + vVortex * cosTj;

        // Project onto panel i normal and tangential
        double sn = uSrcGlobal * sinTi - vSrcGlobal * cosTi;   // -n̂ · v_source = sin·u - cos·v (outward normal)
        double st = uSrcGlobal * cosTi + vSrcGlobal * sinTi;
        double vn = uVrtGlobal * sinTi - vVrtGlobal * cosTi;
        double vt2 = uVrtGlobal * cosTi + vVrtGlobal * sinTi;

        return (sn, st, vn, vt2);
    }

    /// <summary>
    /// Computes velocity (u, v) at a field point due to a source panel and vortex panel.
    /// Returns (u_source, v_source, u_vortex, v_vortex) in global coordinates.
    /// </summary>
    private static (double uSrc, double vSrc, double uVrt, double vVrt) PanelVelocity(
        double xp, double yp, double xj1, double yj1, double xj2, double yj2)
    {
        double dxPanel = xj2 - xj1;
        double dyPanel = yj2 - yj1;
        double sj = Math.Sqrt(dxPanel * dxPanel + dyPanel * dyPanel);
        double cosTj = dxPanel / sj;
        double sinTj = dyPanel / sj;

        double dx1 = xp - xj1;
        double dy1 = yp - yj1;

        double xi = dx1 * cosTj + dy1 * sinTj;
        double eta = -dx1 * sinTj + dy1 * cosTj;

        double r1sq = dx1 * dx1 + dy1 * dy1;
        double dx2 = xp - xj2;
        double dy2 = yp - yj2;
        double r2sq = dx2 * dx2 + dy2 * dy2;

        double logTerm = (r1sq < 1e-28 || r2sq < 1e-28) ? 0.0 : Math.Log(r1sq / r2sq);
        double thetaTerm = Math.Atan2(eta * sj, eta * eta + xi * xi - xi * sj);

        double uLocal = logTerm / (4.0 * Math.PI);
        double vLocal = thetaTerm / (2.0 * Math.PI);

        double uSrc = uLocal * cosTj - vLocal * sinTj;
        double vSrc = uLocal * sinTj + vLocal * cosTj;
        double uVrt = thetaTerm / (2.0 * Math.PI) * cosTj - (-logTerm / (4.0 * Math.PI)) * sinTj;
        double vVrt = thetaTerm / (2.0 * Math.PI) * sinTj + (-logTerm / (4.0 * Math.PI)) * cosTj;

        return (uSrc, vSrc, uVrt, vVrt);
    }

    private static double[] SolveLinearSystem(double[,] A, double[] b, int n)
    {
        // Gaussian elimination with partial pivoting
        var aug = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = A[i, j];
            aug[i, n] = b[i];
        }

        for (int col = 0; col < n; col++)
        {
            // Partial pivoting
            int maxRow = col;
            double maxVal = Math.Abs(aug[col, col]);
            for (int row = col + 1; row < n; row++)
            {
                if (Math.Abs(aug[row, col]) > maxVal)
                {
                    maxVal = Math.Abs(aug[row, col]);
                    maxRow = row;
                }
            }
            if (maxRow != col)
            {
                for (int j = col; j <= n; j++)
                {
                    double tmp = aug[col, j];
                    aug[col, j] = aug[maxRow, j];
                    aug[maxRow, j] = tmp;
                }
            }

            double pivot = aug[col, col];
            if (Math.Abs(pivot) < 1e-14)
                continue; // singular row, skip

            for (int row = col + 1; row < n; row++)
            {
                double factor = aug[row, col] / pivot;
                for (int j = col; j <= n; j++)
                    aug[row, j] -= factor * aug[col, j];
            }
        }

        // Back substitution
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = aug[i, n];
            for (int j = i + 1; j < n; j++)
                sum -= aug[i, j] * x[j];
            x[i] = Math.Abs(aug[i, i]) > 1e-14 ? sum / aug[i, i] : 0.0;
        }

        return x;
    }

    private static double MaxChord(double[] x)
    {
        double min = double.MaxValue;
        double max = double.MinValue;
        for (int i = 0; i < x.Length; i++)
        {
            if (x[i] < min) min = x[i];
            if (x[i] > max) max = x[i];
        }
        return max - min;
    }
}
