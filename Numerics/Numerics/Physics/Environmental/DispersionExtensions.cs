using CSharpNumerics.Physics.Environmental.Enums;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.Environmental;

/// <summary>
/// Extension methods for atmospheric dispersion modelling:
/// Gaussian plume (steady-state), Gaussian puff (transient),
/// and Briggs dispersion parameters by Pasquill–Gifford stability class.
/// </summary>
public static class DispersionExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Gaussian plume
    // ═══════════════════════════════════════════════════════════════

    #region Gaussian Plume

    /// <summary>
    /// Creates a <see cref="ScalarField"/> representing the steady-state concentration
    /// from a Gaussian plume model with user-supplied dispersion functions.
    /// <para>
    /// C(x,y,z) = Q / (2π·u·σy·σz) · exp(−y²/(2σy²))
    ///            · [exp(−(z−H)²/(2σz²)) + exp(−(z+H)²/(2σz²))]
    /// </para>
    /// The second vertical term is the ground-reflection image source.
    /// Returns zero for locations upwind of the source (x ≤ 0).
    /// </summary>
    /// <param name="emissionRate">Source strength Q in kg/s (or consistent mass-flow unit).</param>
    /// <param name="windSpeed">Mean wind speed u at effective stack height in m/s.</param>
    /// <param name="stackHeight">Effective stack height H above ground in metres.</param>
    /// <param name="sourcePosition">Position of the emission source.</param>
    /// <param name="windDirection">Horizontal wind direction (z component is ignored).</param>
    /// <param name="sigmaY">Lateral dispersion σy(x) as a function of downwind distance x.</param>
    /// <param name="sigmaZ">Vertical dispersion σz(x) as a function of downwind distance x.</param>
    public static ScalarField GaussianPlume(
        this double emissionRate,
        double windSpeed,
        double stackHeight,
        Vector sourcePosition,
        Vector windDirection,
        Func<double, double> sigmaY,
        Func<double, double> sigmaZ)
    {
        if (windSpeed <= 0) throw new ArgumentException("Wind speed must be greater than zero.");

        var windUnit = new Vector(windDirection.x, windDirection.y, 0).GetUnitVector();
        var crossUnit = new Vector(-windUnit.y, windUnit.x, 0);
        var Q = emissionRate;
        var u = windSpeed;
        var H = stackHeight;
        var src = sourcePosition;

        return new ScalarField(r =>
        {
            var d = r - src;
            double x = d.x * windUnit.x + d.y * windUnit.y;

            if (x <= 0) return 0;

            double y = d.x * crossUnit.x + d.y * crossUnit.y;
            double z = r.z;

            double sy = sigmaY(x);
            double sz = sigmaZ(x);

            if (sy <= 0 || sz <= 0) return 0;

            double lateral = Math.Exp(-y * y / (2 * sy * sy));
            double vertDirect = Math.Exp(-(z - H) * (z - H) / (2 * sz * sz));
            double vertReflect = Math.Exp(-(z + H) * (z + H) / (2 * sz * sz));

            return Q / (2 * Math.PI * u * sy * sz) * lateral * (vertDirect + vertReflect);
        });
    }

    /// <summary>
    /// Creates a <see cref="ScalarField"/> representing the steady-state concentration
    /// from a Gaussian plume model using Briggs dispersion parameters for the
    /// given Pasquill–Gifford <see cref="StabilityClass"/>.
    /// </summary>
    /// <param name="emissionRate">Source strength Q in kg/s.</param>
    /// <param name="windSpeed">Mean wind speed u at effective stack height in m/s.</param>
    /// <param name="stackHeight">Effective stack height H above ground in metres.</param>
    /// <param name="sourcePosition">Position of the emission source.</param>
    /// <param name="windDirection">Horizontal wind direction (z component is ignored).</param>
    /// <param name="stability">Pasquill–Gifford stability class (default: D — neutral).</param>
    public static ScalarField GaussianPlume(
        this double emissionRate,
        double windSpeed,
        double stackHeight,
        Vector sourcePosition,
        Vector windDirection,
        StabilityClass stability = StabilityClass.D)
    {
        return emissionRate.GaussianPlume(
            windSpeed, stackHeight, sourcePosition, windDirection,
            x => BriggsSigmaY(x, stability),
            x => BriggsSigmaZ(x, stability));
    }

    /// <summary>
    /// Computes the ground-level concentration directly downwind at distance x
    /// from a Gaussian plume (y = 0, z = 0):
    /// C = Q / (π·u·σy·σz) · exp(−H²/(2σz²)).
    /// </summary>
    /// <param name="emissionRate">Source strength Q in kg/s.</param>
    /// <param name="windSpeed">Wind speed in m/s.</param>
    /// <param name="stackHeight">Effective stack height H in metres.</param>
    /// <param name="downwindDistance">Distance x downwind in metres.</param>
    /// <param name="stability">Pasquill–Gifford stability class.</param>
    public static double GaussianPlumeGroundLevel(
        this double emissionRate,
        double windSpeed,
        double stackHeight,
        double downwindDistance,
        StabilityClass stability = StabilityClass.D)
    {
        if (downwindDistance <= 0) return 0;
        double sy = BriggsSigmaY(downwindDistance, stability);
        double sz = BriggsSigmaZ(downwindDistance, stability);
        if (sy <= 0 || sz <= 0) return 0;

        return emissionRate / (Math.PI * windSpeed * sy * sz)
            * Math.Exp(-stackHeight * stackHeight / (2 * sz * sz));
    }

    /// <summary>
    /// Creates a <see cref="ScalarField"/> representing the transient concentration
    /// from a Gaussian puff model at a specific time after an instantaneous release
    /// of mass Q·Δt, advected downwind by a uniform wind field:
    /// <para>
    /// C(r,t) = (Q·Δt) / ((2π)^(3/2) · σx · σy · σz)
    ///          · exp(−(x−u·t)²/(2σx²)) · exp(−y²/(2σy²))
    ///          · [exp(−(z−H)²/(2σz²)) + exp(−(z+H)²/(2σz²))]
    /// </para>
    /// where σx = σy (isotropic lateral dispersion evaluated at the virtual
    /// downwind distance u·t) and σz is similarly evaluated.
    /// </summary>
    /// <param name="emissionRate">Source strength Q in kg/s.</param>
    /// <param name="releaseSeconds">Duration of release Δt in seconds (mass = Q·Δt).</param>
    /// <param name="windSpeed">Mean wind speed u in m/s.</param>
    /// <param name="stackHeight">Effective stack height H in metres.</param>
    /// <param name="sourcePosition">Position of the emission source.</param>
    /// <param name="windDirection">Horizontal wind direction (z ignored).</param>
    /// <param name="time">Time since start of release in seconds.</param>
    /// <param name="stability">Pasquill–Gifford stability class.</param>
    public static ScalarField GaussianPuff(
        this double emissionRate,
        double releaseSeconds,
        double windSpeed,
        double stackHeight,
        Vector sourcePosition,
        Vector windDirection,
        double time,
        StabilityClass stability = StabilityClass.D)
    {
        if (time <= 0) return new ScalarField(_ => 0);
        if (windSpeed <= 0) throw new ArgumentException("Wind speed must be greater than zero.");

        var windUnit = new Vector(windDirection.x, windDirection.y, 0).GetUnitVector();
        var crossUnit = new Vector(-windUnit.y, windUnit.x, 0);
        double mass = emissionRate * releaseSeconds;
        double u = windSpeed;
        double H = stackHeight;
        double t = time;

        // Virtual downwind distance for dispersion evaluation
        double xVirt = u * t;
        double sy = BriggsSigmaY(xVirt, stability);
        double sz = BriggsSigmaZ(xVirt, stability);
        double sx = sy; // isotropic lateral

        if (sx <= 0 || sy <= 0 || sz <= 0) return new ScalarField(_ => 0);

        double norm = mass / (Math.Pow(2 * Math.PI, 1.5) * sx * sy * sz);

        return new ScalarField(r =>
        {
            var d = r - sourcePosition;
            double xDownwind = d.x * windUnit.x + d.y * windUnit.y;
            double yCross = d.x * crossUnit.x + d.y * crossUnit.y;
            double z = r.z;

            double dx = xDownwind - u * t;
            double axial = Math.Exp(-dx * dx / (2 * sx * sx));
            double lateral = Math.Exp(-yCross * yCross / (2 * sy * sy));
            double vertDirect = Math.Exp(-(z - H) * (z - H) / (2 * sz * sz));
            double vertReflect = Math.Exp(-(z + H) * (z + H) / (2 * sz * sz));

            return norm * axial * lateral * (vertDirect + vertReflect);
        });
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Briggs dispersion parameters
    // ═══════════════════════════════════════════════════════════════

    #region Dispersion Parameters

    /// <summary>
    /// Computes the Briggs lateral dispersion parameter σy(x)
    /// for the given Pasquill–Gifford stability class.
    /// </summary>
    /// <param name="x">Downwind distance in metres.</param>
    /// <param name="stability">Stability class.</param>
    public static double BriggsSigmaY(double x, StabilityClass stability)
    {
        if (x <= 0) return 0;

        double a;
        switch (stability)
        {
            case StabilityClass.A: a = 0.22; break;
            case StabilityClass.B: a = 0.16; break;
            case StabilityClass.C: a = 0.11; break;
            case StabilityClass.D: a = 0.08; break;
            case StabilityClass.E: a = 0.06; break;
            case StabilityClass.F: a = 0.04; break;
            default: a = 0.08; break;
        }

        return a * x / Math.Sqrt(1 + 0.0001 * x);
    }

    /// <summary>
    /// Computes the Briggs vertical dispersion parameter σz(x)
    /// for the given Pasquill–Gifford stability class.
    /// </summary>
    /// <param name="x">Downwind distance in metres.</param>
    /// <param name="stability">Stability class.</param>
    public static double BriggsSigmaZ(double x, StabilityClass stability)
    {
        if (x <= 0) return 0;

        switch (stability)
        {
            case StabilityClass.A: return 0.20 * x;
            case StabilityClass.B: return 0.12 * x;
            case StabilityClass.C: return 0.08 * x / Math.Sqrt(1 + 0.0002 * x);
            case StabilityClass.D: return 0.06 * x / Math.Sqrt(1 + 0.0015 * x);
            case StabilityClass.E: return 0.03 * x / (1 + 0.0003 * x);
            case StabilityClass.F: return 0.016 * x / (1 + 0.0003 * x);
            default:               return 0.06 * x / Math.Sqrt(1 + 0.0015 * x);
        }
    }

    #endregion
}
