using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Generates 2D airfoil surface coordinates for NACA 4-digit series airfoils.
/// Points are ordered clockwise from trailing edge (upper surface → leading edge → lower surface → trailing edge).
/// </summary>
public static class NACAGeometry
{
    /// <summary>
    /// Generates surface panel coordinates for a NACA 4-digit airfoil.
    /// </summary>
    /// <param name="naca">4-character NACA designation (e.g. "0012", "2412").</param>
    /// <param name="numPanels">Total number of panels (even number). Default 100.</param>
    /// <param name="chord">Chord length in metres. Default 1.0.</param>
    /// <returns>Arrays of (x, y) coordinates. Length = numPanels + 1 (closed contour).</returns>
    public static (double[] x, double[] y) Generate(string naca, int numPanels = 100, double chord = 1.0)
    {
        if (naca == null || naca.Length != 4)
            throw new ArgumentException("NACA designation must be a 4-character string (e.g. \"0012\").");
        if (numPanels < 4 || numPanels % 2 != 0)
            throw new ArgumentException("Number of panels must be an even integer ≥ 4.");

        double m = (naca[0] - '0') / 100.0;  // max camber
        double p = (naca[1] - '0') / 10.0;   // location of max camber
        double t = int.Parse(naca.Substring(2, 2)) / 100.0; // thickness ratio

        int halfPanels = numPanels / 2;
        int nPoints = numPanels + 1;
        var xCoords = new double[nPoints];
        var yCoords = new double[nPoints];

        // Cosine spacing for better leading-edge resolution
        var xc = new double[halfPanels + 1];
        for (int i = 0; i <= halfPanels; i++)
        {
            double beta = Math.PI * i / halfPanels;
            xc[i] = 0.5 * (1.0 - Math.Cos(beta));
        }

        // Upper surface: from trailing edge to leading edge
        for (int i = 0; i <= halfPanels; i++)
        {
            double x = xc[halfPanels - i];
            double yt = ThicknessDistribution(x, t);
            double yc = CamberLine(x, m, p);
            double theta = CamberSlope(x, m, p);

            xCoords[i] = (x - yt * Math.Sin(theta)) * chord;
            yCoords[i] = (yc + yt * Math.Cos(theta)) * chord;
        }

        // Lower surface: from leading edge to trailing edge
        for (int i = 1; i <= halfPanels; i++)
        {
            double x = xc[i];
            double yt = ThicknessDistribution(x, t);
            double yc = CamberLine(x, m, p);
            double theta = CamberSlope(x, m, p);

            xCoords[halfPanels + i] = (x + yt * Math.Sin(theta)) * chord;
            yCoords[halfPanels + i] = (yc - yt * Math.Cos(theta)) * chord;
        }

        return (xCoords, yCoords);
    }

    /// <summary>
    /// Generates surface coordinates for a symmetric NACA airfoil (e.g. NACA 0012).
    /// </summary>
    /// <param name="thicknessPercent">Thickness as percent of chord (e.g. 12 for NACA 0012).</param>
    /// <param name="numPanels">Total number of panels (even number).</param>
    /// <param name="chord">Chord length in metres.</param>
    public static (double[] x, double[] y) GenerateSymmetric(
        int thicknessPercent = 12, int numPanels = 100, double chord = 1.0)
    {
        string naca = $"00{thicknessPercent:D2}";
        return Generate(naca, numPanels, chord);
    }

    /// <summary>
    /// Rotates airfoil coordinates by a given angle of attack about the quarter-chord point.
    /// </summary>
    /// <param name="x">X coordinates.</param>
    /// <param name="y">Y coordinates.</param>
    /// <param name="alpha">Angle of attack in radians (positive = nose up).</param>
    /// <param name="chord">Chord length for locating the rotation centre (0.25c).</param>
    public static (double[] xRot, double[] yRot) Rotate(
        double[] x, double[] y, double alpha, double chord = 1.0)
    {
        double cx = 0.25 * chord;
        double cy = 0.0;
        double cosA = Math.Cos(-alpha);
        double sinA = Math.Sin(-alpha);

        var xRot = new double[x.Length];
        var yRot = new double[y.Length];

        for (int i = 0; i < x.Length; i++)
        {
            double dx = x[i] - cx;
            double dy = y[i] - cy;
            xRot[i] = cx + dx * cosA - dy * sinA;
            yRot[i] = cy + dx * sinA + dy * cosA;
        }

        return (xRot, yRot);
    }

    // ═══════════════════════════════════════════════════════════════
    //  Internal NACA formulas
    // ═══════════════════════════════════════════════════════════════

    private static double ThicknessDistribution(double x, double t)
    {
        // NACA 4-digit thickness formula (open trailing edge variant)
        return 5.0 * t * (
            0.2969 * Math.Sqrt(x)
            - 0.1260 * x
            - 0.3516 * x * x
            + 0.2843 * x * x * x
            - 0.1015 * x * x * x * x);
    }

    private static double CamberLine(double x, double m, double p)
    {
        if (m == 0 || p == 0) return 0.0;

        if (x < p)
            return m / (p * p) * (2.0 * p * x - x * x);
        else
            return m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x - x * x);
    }

    private static double CamberSlope(double x, double m, double p)
    {
        if (m == 0 || p == 0) return 0.0;

        if (x < p)
            return Math.Atan(2.0 * m / (p * p) * (p - x));
        else
            return Math.Atan(2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x));
    }
}
