using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Optics;

/// <summary>
/// A triangular prism defined by its apex angle, material, and orientation.
/// Implements dispersion: the refractive index varies with wavelength via the
/// <see cref="OpticalMedium.RefractiveIndexAt"/> method.
/// </summary>
public class Prism
{
    /// <summary>Apex (top) angle of the prism in radians.</summary>
    public double ApexAngle { get; }

    /// <summary>Optical medium of the prism glass.</summary>
    public OpticalMedium Medium { get; }

    /// <summary>
    /// Creates a prism.
    /// </summary>
    /// <param name="apexAngleRadians">Apex angle A in radians.</param>
    /// <param name="medium">Glass material with dispersion data.</param>
    public Prism(double apexAngleRadians, OpticalMedium medium)
    {
        ApexAngle = apexAngleRadians;
        Medium = medium;
    }

    /// <summary>
    /// Computes the total deviation angle δ for a ray entering at angle θ₁
    /// (measured from the first surface normal).
    /// <para>
    /// δ = θ₁ + θ₄ − A, where θ₄ is the exit angle from the second surface.
    /// </para>
    /// Returns null if total internal reflection occurs inside the prism.
    /// </summary>
    /// <param name="thetaIncident">Incidence angle at the first surface (radians).</param>
    /// <param name="wavelengthNm">Wavelength for dispersion (nm).</param>
    /// <param name="nOutside">Refractive index of the surrounding medium (default air).</param>
    public double? DeviationAngle(double thetaIncident, double wavelengthNm = 550.0,
        double nOutside = 1.000293)
    {
        double nPrism = Medium.RefractiveIndexAt(wavelengthNm);

        // Refract at first surface
        var theta2 = OpticsExtensions.RefractionAngle(nOutside, nPrism, thetaIncident);
        if (!theta2.HasValue) return null;

        // Angle of incidence on second surface
        double theta3 = ApexAngle - theta2.Value;
        if (theta3 < 0) return null; // ray exits through base

        // Refract at second surface
        var theta4 = OpticsExtensions.RefractionAngle(nPrism, nOutside, theta3);
        if (!theta4.HasValue) return null; // TIR inside prism

        return thetaIncident + theta4.Value - ApexAngle;
    }

    /// <summary>
    /// Minimum deviation angle for this prism at a given wavelength.
    /// δ_min = 2 · arcsin(n · sin(A/2)) − A.
    /// </summary>
    public double MinimumDeviation(double wavelengthNm = 550.0)
    {
        double n = Medium.RefractiveIndexAt(wavelengthNm);
        return 2.0 * Math.Asin(n * Math.Sin(ApexAngle / 2.0)) - ApexAngle;
    }

    /// <summary>
    /// Angular dispersion: the difference in deviation between two wavelengths.
    /// </summary>
    public double AngularDispersion(double lambda1Nm, double lambda2Nm,
        double thetaIncident, double nOutside = 1.000293)
    {
        double? d1 = DeviationAngle(thetaIncident, lambda1Nm, nOutside);
        double? d2 = DeviationAngle(thetaIncident, lambda2Nm, nOutside);
        if (!d1.HasValue || !d2.HasValue)
            throw new InvalidOperationException("TIR occurred for one of the wavelengths.");
        return Math.Abs(d1.Value - d2.Value);
    }
}
