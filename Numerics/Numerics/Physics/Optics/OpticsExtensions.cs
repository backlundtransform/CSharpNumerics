using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// Extension methods for geometric and physical optics calculations:
/// Snell's law, Fresnel equations, total internal reflection, and
/// Beer–Lambert absorption.
/// </summary>
public static class OpticsExtensions
{
    // ───────────────────── Snell's Law ─────────────────────

    /// <summary>
    /// Computes the sine of the refracted angle using Snell's law:
    /// n₁ sin θ₁ = n₂ sin θ₂.
    /// Returns null when total internal reflection occurs.
    /// </summary>
    public static double? SnellSinRefracted(double n1, double n2, double sinTheta1)
    {
        double sinTheta2 = n1 / n2 * sinTheta1;
        if (Math.Abs(sinTheta2) > 1.0) return null; // TIR
        return sinTheta2;
    }

    /// <summary>
    /// Computes the refraction angle in radians given incidence angle and refractive indices.
    /// Returns null for total internal reflection.
    /// </summary>
    public static double? RefractionAngle(double n1, double n2, double thetaIncident)
    {
        var sin2 = SnellSinRefracted(n1, n2, Math.Sin(thetaIncident));
        return sin2.HasValue ? Math.Asin(sin2.Value) : null;
    }

    /// <summary>
    /// Returns the critical angle for total internal reflection (radians).
    /// Only valid when n1 &gt; n2.
    /// </summary>
    public static double CriticalAngle(double n1, double n2)
    {
        if (n1 <= n2) throw new ArgumentException("Critical angle only exists when n1 > n2.");
        return Math.Asin(n2 / n1);
    }

    /// <summary>
    /// Returns true when total internal reflection occurs.
    /// </summary>
    public static bool IsTotalInternalReflection(double n1, double n2, double thetaIncident) =>
        n1 > n2 && Math.Sin(thetaIncident) >= n2 / n1;

    // ───────────────────── Fresnel Equations ─────────────────────

    /// <summary>
    /// Fresnel reflectance for s-polarised light (TE).
    /// </summary>
    public static double FresnelReflectanceS(double n1, double n2, double thetaI, double thetaT) =>
        Sqr((n1 * Math.Cos(thetaI) - n2 * Math.Cos(thetaT)) /
            (n1 * Math.Cos(thetaI) + n2 * Math.Cos(thetaT)));

    /// <summary>
    /// Fresnel reflectance for p-polarised light (TM).
    /// </summary>
    public static double FresnelReflectanceP(double n1, double n2, double thetaI, double thetaT) =>
        Sqr((n2 * Math.Cos(thetaI) - n1 * Math.Cos(thetaT)) /
            (n2 * Math.Cos(thetaI) + n1 * Math.Cos(thetaT)));

    /// <summary>
    /// Average Fresnel reflectance for unpolarised light: R = (R_s + R_p) / 2.
    /// Returns 1.0 for total internal reflection.
    /// </summary>
    public static double FresnelReflectance(double n1, double n2, double thetaIncident)
    {
        var thetaT = RefractionAngle(n1, n2, thetaIncident);
        if (!thetaT.HasValue) return 1.0; // TIR

        double rs = FresnelReflectanceS(n1, n2, thetaIncident, thetaT.Value);
        double rp = FresnelReflectanceP(n1, n2, thetaIncident, thetaT.Value);
        return 0.5 * (rs + rp);
    }

    /// <summary>
    /// Schlick's fast approximation of Fresnel reflectance (used in real-time rendering).
    /// R(θ) ≈ R₀ + (1 − R₀)(1 − cos θ)⁵  where  R₀ = ((n₁−n₂)/(n₁+n₂))².
    /// </summary>
    public static double SchlickReflectance(double n1, double n2, double cosTheta)
    {
        double r0 = Sqr((n1 - n2) / (n1 + n2));
        double oneMinusCos = 1.0 - cosTheta;
        return r0 + (1.0 - r0) * oneMinusCos * oneMinusCos * oneMinusCos * oneMinusCos * oneMinusCos;
    }

    // ───────────────────── Beer–Lambert Absorption ─────────────────────

    /// <summary>
    /// Beer–Lambert transmittance: T = exp(−α · d).
    /// </summary>
    /// <param name="absorptionCoefficient">Absorption coefficient α in 1/m.</param>
    /// <param name="distance">Path length through the medium in metres.</param>
    public static double BeerLambertTransmittance(double absorptionCoefficient, double distance) =>
        Math.Exp(-absorptionCoefficient * distance);

    /// <summary>
    /// Returns the ray intensity after travelling through a medium using Beer–Lambert law.
    /// </summary>
    public static double AttenuateIntensity(double intensity, double absorptionCoefficient, double distance) =>
        intensity * BeerLambertTransmittance(absorptionCoefficient, distance);

    // ───────────────────── Ray-Level Helpers ─────────────────────

    /// <summary>
    /// Reflects a ray direction about a surface normal.
    /// Both <paramref name="incident"/> and <paramref name="normal"/> must be unit vectors.
    /// Returns the reflected direction (unit vector).
    /// </summary>
    public static Vector Reflect(Vector incident, Vector normal)
    {
        double dot = incident.Dot(normal);
        return (incident - 2.0 * dot * normal).GetUnitVector();
    }

    /// <summary>
    /// Refracts a ray direction through a surface using Snell's law.
    /// Returns null when total internal reflection occurs.
    /// Both <paramref name="incident"/> and <paramref name="normal"/> must be unit vectors.
    /// </summary>
    public static Vector? Refract(Vector incident, Vector normal, double n1, double n2)
    {
        double cosI = -incident.Dot(normal);
        double ratio = n1 / n2;
        double sin2T = ratio * ratio * (1.0 - cosI * cosI);
        if (sin2T > 1.0) return null; // TIR

        double cosT = Math.Sqrt(1.0 - sin2T);
        Vector refracted = ratio * incident + (ratio * cosI - cosT) * normal;
        return refracted.GetUnitVector();
    }

    // ───────────────────── Private ─────────────────────

    private static double Sqr(double x) => x * x;
}
