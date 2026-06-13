using System;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Physics.Astro;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Engines.Exoplanet;

/// <summary>
/// Analytical transit light curve model based on Mandel &amp; Agol (2002).
/// Computes the fractional flux drop as a planet transits a limb-darkened star.
/// Supports circular orbits with uniform, linear, and quadratic limb darkening.
/// </summary>
public class TransitModel
{
    private const double SolarRadiusMeters = 6.957e8;

    /// <summary>
    /// Evaluates the transit model flux at each given time.
    /// </summary>
    /// <param name="times">Array of observation times (days, BJD).</param>
    /// <param name="p">Transit parameters (Period, Epoch, RadiusRatio, ImpactParameter, Duration).</param>
    /// <param name="ldModel">Limb darkening model type.</param>
    /// <param name="ldCoefficients">Limb darkening coefficients.</param>
    /// <param name="star">Stellar properties (Radius in solar radii, Mass in solar masses).</param>
    /// <returns>Model flux array (1.0 = no transit, &lt;1.0 = in transit).</returns>
    public double[] Evaluate(double[] times, TransitParameters p, LimbDarkeningModel ldModel,
        double[] ldCoefficients, StellarProperties star)
    {
        if (times == null) throw new ArgumentNullException(nameof(times));
        if (p == null) throw new ArgumentNullException(nameof(p));

        double k = p.RadiusRatio; // Rp/R★
        double b = p.ImpactParameter;
        double period = p.Period;
        double epoch = p.Epoch;

        // Compute a/R★ from transit parameters
        double aOverRstar = ComputeAOverRstar(period, star);

        double inclination = Math.Acos(b / aOverRstar);

        double[] flux = new double[times.Length];

        for (int i = 0; i < times.Length; i++)
        {
            // Phase relative to mid-transit
            double phase = ((times[i] - epoch) % period + period) % period;
            if (phase > period / 2.0)
                phase -= period;

            // Projected separation z (in units of R★)
            double z = ProjectedSeparation(phase, period, aOverRstar, inclination);

            // Compute flux drop
            flux[i] = 1.0 - FluxDrop(z, k, ldModel, ldCoefficients);
        }

        return flux;
    }

    /// <summary>
    /// Evaluates the transit model with a/R★ specified directly (no stellar properties needed).
    /// </summary>
    public double[] Evaluate(double[] times, double period, double epoch, double radiusRatio,
        double aOverRstar, double inclination, LimbDarkeningModel ldModel, double[] ldCoefficients)
    {
        if (times == null) throw new ArgumentNullException(nameof(times));

        double k = radiusRatio;
        double[] flux = new double[times.Length];

        for (int i = 0; i < times.Length; i++)
        {
            double phase = ((times[i] - epoch) % period + period) % period;
            if (phase > period / 2.0)
                phase -= period;

            double z = ProjectedSeparation(phase, period, aOverRstar, inclination);
            flux[i] = 1.0 - FluxDrop(z, k, ldModel, ldCoefficients);
        }

        return flux;
    }

    /// <summary>
    /// Computes a/R★ from period (days) and stellar properties.
    /// </summary>
    private static double ComputeAOverRstar(double periodDays, StellarProperties star)
    {
        double periodSec = periodDays * 86400.0;
        double massSI = star.Mass * PhysicsConstants.SolarMass;
        double rStarSI = star.Radius * SolarRadiusMeters;

        double G = PhysicsConstants.GravitationalConstant;
        double a = Math.Pow(G * massSI * periodSec * periodSec / (4.0 * Math.PI * Math.PI), 1.0 / 3.0);

        return a / rStarSI;
    }

    /// <summary>
    /// Projected star-planet separation in units of R★ for circular orbit.
    /// z = (a/R★) · sqrt( sin²(ωt) + cos²(i)·cos²(ωt) )
    /// where ω = 2π/P, t = time from mid-transit.
    /// Returns double.MaxValue if the planet is behind the star (secondary eclipse).
    /// </summary>
    private static double ProjectedSeparation(double phase, double period, double aOverRstar, double inclination)
    {
        double omega = 2.0 * Math.PI / period;
        double t = phase;
        double sinOt = Math.Sin(omega * t);
        double cosOt = Math.Cos(omega * t);
        double cosI = Math.Cos(inclination);

        // Planet is behind the star when cos(ωt) < 0
        if (cosOt < 0)
            return double.MaxValue;

        return aOverRstar * Math.Sqrt(sinOt * sinOt + cosI * cosI * cosOt * cosOt);
    }

    /// <summary>
    /// Computes the fractional flux blocked by the planet at projected separation z.
    /// Based on Mandel &amp; Agol (2002) analytical expressions.
    /// </summary>
    private static double FluxDrop(double z, double k, LimbDarkeningModel ldModel, double[] coefficients)
    {
        if (z >= 1.0 + k)
            return 0.0; // No overlap

        switch (ldModel)
        {
            case LimbDarkeningModel.Uniform:
                return UniformFluxDrop(z, k);

            case LimbDarkeningModel.Linear:
                return LinearFluxDrop(z, k, coefficients[0]);

            case LimbDarkeningModel.Quadratic:
                return QuadraticFluxDrop(z, k, coefficients[0], coefficients[1]);

            default:
                return UniformFluxDrop(z, k);
        }
    }

    /// <summary>
    /// Uniform limb darkening: flux drop = area of overlap / π.
    /// Mandel &amp; Agol Eq. 1 (uniform source).
    /// </summary>
    private static double UniformFluxDrop(double z, double k)
    {
        if (z >= 1.0 + k)
            return 0.0;

        if (z <= 1.0 - k)
        {
            // Planet fully on disk
            return k * k;
        }

        if (z <= k - 1.0)
        {
            // Star fully inside planet (k > 1)
            return 1.0;
        }

        // Partial overlap — area of lens-shaped intersection
        double kappa1 = Math.Acos((1.0 - k * k + z * z) / (2.0 * z));
        double kappa0 = Math.Acos((k * k + z * z - 1.0) / (2.0 * k * z));

        double area = (1.0 / Math.PI) * (k * k * kappa0 + kappa1
            - 0.5 * Math.Sqrt(Math.Max(0, 4.0 * z * z - (1.0 + z * z - k * k) * (1.0 + z * z - k * k))));

        return area;
    }

    /// <summary>
    /// Linear limb darkening transit using numerical integration over the occulted stellar disk.
    /// I(μ) = 1 − u1·(1 − μ), integrated over the overlap.
    /// </summary>
    private static double LinearFluxDrop(double z, double k, double u1)
    {
        // Use the correction factor approach:
        // ΔF_LD = ΔF_uniform · (1 − u1) + u1 · I_LD_occult(z, k)
        // where I_LD_occult is the limb-darkened contribution from the occulted region.

        double uniformDrop = UniformFluxDrop(z, k);
        if (uniformDrop <= 0) return 0;

        // Normalization: total flux of limb-darkened star = 1 − u1/3
        double norm = 1.0 - u1 / 3.0;

        // For fully in-transit (z + k <= 1), compute the LD correction analytically
        if (z <= 1.0 - k)
        {
            // Occulted intensity at planet position
            double mu = Math.Sqrt(Math.Max(0, 1.0 - z * z));
            double intensity = LimbDarkening.Linear(mu, u1);
            return k * k * intensity / norm;
        }

        // For partial overlap, use numerical integration
        return NumericalLDFluxDrop(z, k, LimbDarkeningModel.Linear, new[] { u1 });
    }

    /// <summary>
    /// Quadratic limb darkening transit using numerical integration.
    /// I(μ) = 1 − u1·(1−μ) − u2·(1−μ)².
    /// </summary>
    private static double QuadraticFluxDrop(double z, double k, double u1, double u2)
    {
        double uniformDrop = UniformFluxDrop(z, k);
        if (uniformDrop <= 0) return 0;

        double norm = 1.0 - u1 / 3.0 - u2 / 6.0;

        if (z <= 1.0 - k)
        {
            double mu = Math.Sqrt(Math.Max(0, 1.0 - z * z));
            double intensity = LimbDarkening.Quadratic(mu, u1, u2);
            return k * k * intensity / norm;
        }

        return NumericalLDFluxDrop(z, k, LimbDarkeningModel.Quadratic, new[] { u1, u2 });
    }

    /// <summary>
    /// Numerical integration of limb-darkened flux drop for partial overlap cases.
    /// Uses a radial-angular integration scheme over the planet disk.
    /// </summary>
    private static double NumericalLDFluxDrop(double z, double k, LimbDarkeningModel model, double[] coeff)
    {
        double norm;
        switch (model)
        {
            case LimbDarkeningModel.Linear:
                norm = 1.0 - coeff[0] / 3.0;
                break;
            case LimbDarkeningModel.Quadratic:
                norm = 1.0 - coeff[0] / 3.0 - coeff[1] / 6.0;
                break;
            default:
                norm = 1.0;
                break;
        }

        // Integrate over the planet disk using polar coordinates centered on the planet
        int nR = 50;
        int nTheta = 100;
        double sum = 0;
        double dr = k / nR;

        for (int ir = 0; ir < nR; ir++)
        {
            double r = (ir + 0.5) * dr;
            double dTheta = 2.0 * Math.PI / nTheta;

            for (int it = 0; it < nTheta; it++)
            {
                double theta = (it + 0.5) * dTheta;

                // Position on stellar disk
                double x = z + r * Math.Cos(theta);
                double y = r * Math.Sin(theta);
                double rStar = Math.Sqrt(x * x + y * y);

                if (rStar > 1.0) continue; // Outside stellar disk

                double mu = Math.Sqrt(1.0 - rStar * rStar);
                double intensity;

                switch (model)
                {
                    case LimbDarkeningModel.Linear:
                        intensity = LimbDarkening.Linear(mu, coeff[0]);
                        break;
                    case LimbDarkeningModel.Quadratic:
                        intensity = LimbDarkening.Quadratic(mu, coeff[0], coeff[1]);
                        break;
                    default:
                        intensity = 1.0;
                        break;
                }

                sum += intensity * r * dr * dTheta;
            }
        }

        return sum / (Math.PI * norm);
    }
}
