using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Astro;

public static class KeplerOrbit
{
    /// <summary>
    /// Solves Kepler's equation M = E − e·sin(E) for the eccentric anomaly E
    /// using Newton-Raphson iteration, then converts to true anomaly ν.
    /// </summary>
    /// <param name="meanAnomaly">Mean anomaly M in radians.</param>
    /// <param name="eccentricity">Orbital eccentricity (0 ≤ e &lt; 1).</param>
    /// <param name="tolerance">Convergence tolerance (default 1e-12).</param>
    /// <param name="maxIterations">Maximum iterations (default 100).</param>
    /// <returns>True anomaly ν in radians.</returns>
    public static double TrueAnomaly(double meanAnomaly, double eccentricity, double tolerance = 1e-12, int maxIterations = 100)
    {
        if (eccentricity < 0 || eccentricity >= 1)
            throw new ArgumentOutOfRangeException(nameof(eccentricity), "Eccentricity must be in [0, 1).");

        // Solve Kepler's equation: M = E - e·sin(E)
        double E = meanAnomaly; // initial guess
        for (int i = 0; i < maxIterations; i++)
        {
            double dE = (E - eccentricity * Math.Sin(E) - meanAnomaly) / (1.0 - eccentricity * Math.Cos(E));
            E -= dE;
            if (Math.Abs(dE) < tolerance)
                break;
        }

        // Convert eccentric anomaly to true anomaly
        double sinNu = Math.Sqrt(1.0 - eccentricity * eccentricity) * Math.Sin(E) / (1.0 - eccentricity * Math.Cos(E));
        double cosNu = (Math.Cos(E) - eccentricity) / (1.0 - eccentricity * Math.Cos(E));

        return Math.Atan2(sinNu, cosNu);
    }

    /// <summary>
    /// Kepler's third law: a = (G·M★·P²/(4π²))^(1/3).
    /// </summary>
    /// <param name="period">Orbital period in seconds.</param>
    /// <param name="stellarMass">Stellar mass in kg.</param>
    /// <returns>Semi-major axis in meters.</returns>
    public static double SemiMajorAxis(double period, double stellarMass)
    {
        if (period <= 0) throw new ArgumentException("Period must be positive.", nameof(period));
        if (stellarMass <= 0) throw new ArgumentException("Stellar mass must be positive.", nameof(stellarMass));

        double G = PhysicsConstants.GravitationalConstant;
        double numerator = G * stellarMass * period * period;
        double denominator = 4.0 * Math.PI * Math.PI;

        return Math.Pow(numerator / denominator, 1.0 / 3.0);
    }

    /// <summary>
    /// Kepler's third law with period in days and mass in solar masses.
    /// Returns semi-major axis in AU.
    /// </summary>
    /// <param name="periodDays">Orbital period in days.</param>
    /// <param name="stellarMassSolar">Stellar mass in solar masses.</param>
    /// <returns>Semi-major axis in AU.</returns>
    public static double SemiMajorAxisAU(double periodDays, double stellarMassSolar)
    {
        double periodSeconds = periodDays * 86400.0;
        double massSI = stellarMassSolar * PhysicsConstants.SolarMass;
        double aSI = SemiMajorAxis(periodSeconds, massSI);
        return aSI / PhysicsConstants.AstronomicalUnit;
    }

    /// <summary>
    /// Mean orbital velocity for a circular orbit: v = 2πa/P.
    /// </summary>
    /// <param name="a">Semi-major axis (any consistent unit).</param>
    /// <param name="period">Orbital period (same time unit).</param>
    /// <returns>Orbital velocity (unit/time).</returns>
    public static double OrbitalVelocity(double a, double period)
    {
        if (period <= 0) throw new ArgumentException("Period must be positive.", nameof(period));
        return 2.0 * Math.PI * a / period;
    }
}
