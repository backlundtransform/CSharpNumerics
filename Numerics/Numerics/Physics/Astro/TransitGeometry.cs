using System;

namespace CSharpNumerics.Physics.Astro;

public static class TransitGeometry
{
    /// <summary>
    /// Computes the impact parameter b = (a / R★) · cos(i).
    /// </summary>
    /// <param name="aOverRstar">Semi-major axis in units of stellar radii (a/R★).</param>
    /// <param name="inclination">Orbital inclination in radians.</param>
    /// <returns>Impact parameter b (dimensionless).</returns>
    public static double ImpactParameter(double aOverRstar, double inclination)
    {
        return aOverRstar * Math.Cos(inclination);
    }

    /// <summary>
    /// Geometric transit probability P = (R★ + Rp) / a.
    /// </summary>
    /// <param name="a">Semi-major axis (any consistent unit).</param>
    /// <param name="rStar">Stellar radius (same unit as a).</param>
    /// <param name="rPlanet">Planet radius (same unit as a).</param>
    /// <returns>Transit probability (0–1).</returns>
    public static double TransitProbability(double a, double rStar, double rPlanet)
    {
        if (a <= 0) throw new ArgumentException("Semi-major axis must be positive.", nameof(a));
        return (rStar + rPlanet) / a;
    }

    /// <summary>
    /// Total transit duration (T14) for a circular orbit.
    /// T14 = (P/π) · arcsin[ (R★/a) · sqrt( (1 + k)² − b² ) / sin(i) ]
    /// where k = Rp/R★, b = impact parameter.
    /// </summary>
    /// <param name="period">Orbital period (days).</param>
    /// <param name="aOverRstar">Semi-major axis in units of stellar radii (a/R★).</param>
    /// <param name="radiusRatio">Planet-to-star radius ratio k = Rp/R★.</param>
    /// <param name="inclination">Orbital inclination in radians.</param>
    /// <returns>Transit duration T14 in same units as period (days).</returns>
    public static double TransitDuration(double period, double aOverRstar, double radiusRatio, double inclination)
    {
        double b = aOverRstar * Math.Cos(inclination);
        double sinI = Math.Sin(inclination);
        double onePlusK = 1.0 + radiusRatio;

        double arg = onePlusK * onePlusK - b * b;
        if (arg < 0) return 0;

        double inner = Math.Sqrt(arg) / (aOverRstar * sinI);
        if (inner > 1.0) inner = 1.0;

        return (period / Math.PI) * Math.Asin(inner);
    }

    /// <summary>
    /// Ingress/egress duration (T12 = T34) for a circular orbit.
    /// T12 = (P/π) · arcsin[ (R★/a) · sqrt( (1 − k)² − b² ) / sin(i) ]  (subtracted from full transit formula)
    /// Approximate: T_ingress ≈ T14 · [ sqrt((1+k)²−b²) − sqrt((1−k)²−b²) ] / [ 2·sqrt((1+k)²−b²) ]
    /// </summary>
    public static double IngressDuration(double period, double aOverRstar, double radiusRatio, double inclination)
    {
        double b = aOverRstar * Math.Cos(inclination);
        double sinI = Math.Sin(inclination);
        double oneMinusK = 1.0 - radiusRatio;

        double argFull = (1.0 + radiusRatio) * (1.0 + radiusRatio) - b * b;
        double argFlat = oneMinusK * oneMinusK - b * b;

        if (argFull < 0) return 0;
        if (argFlat < 0)
        {
            // Grazing transit — no flat bottom, ingress = half of total duration
            return TransitDuration(period, aOverRstar, radiusRatio, inclination) / 2.0;
        }

        double innerFull = Math.Sqrt(argFull) / (aOverRstar * sinI);
        double innerFlat = Math.Sqrt(argFlat) / (aOverRstar * sinI);

        if (innerFull > 1.0) innerFull = 1.0;
        if (innerFlat > 1.0) innerFlat = 1.0;

        double t14 = (period / Math.PI) * Math.Asin(innerFull);
        double t23 = (period / Math.PI) * Math.Asin(innerFlat);

        return (t14 - t23) / 2.0;
    }

    /// <summary>
    /// Contact times T1, T2, T3, T4 relative to mid-transit (epoch).
    /// T1 = first contact (planet touches stellar limb)
    /// T2 = second contact (planet fully on disk)
    /// T3 = third contact (planet starts to leave)
    /// T4 = fourth contact (planet fully off disk)
    /// </summary>
    /// <returns>(T1, T2, T3, T4) as absolute times.</returns>
    public static (double T1, double T2, double T3, double T4) ContactTimes(
        double epoch, double period, double aOverRstar, double radiusRatio, double inclination)
    {
        double t14 = TransitDuration(period, aOverRstar, radiusRatio, inclination);
        double tIngress = IngressDuration(period, aOverRstar, radiusRatio, inclination);

        double t1 = epoch - t14 / 2.0;
        double t2 = t1 + tIngress;
        double t3 = epoch + t14 / 2.0 - tIngress;
        double t4 = epoch + t14 / 2.0;

        return (t1, t2, t3, t4);
    }
}
