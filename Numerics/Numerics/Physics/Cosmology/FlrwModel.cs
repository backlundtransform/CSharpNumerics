using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Cosmology;

/// <summary>
/// A homogeneous, isotropic Friedmann–Lemaître–Robertson–Walker (FLRW) cosmology.
/// Holds the present-day density parameters and Hubble constant and derives the
/// expansion history H(z), the cosmological distance measures, lookback time and
/// the age of the universe by numerical integration of the Friedmann equation.
///
/// Conventions: H₀ is supplied in km·s⁻¹·Mpc⁻¹; distances are returned in metres,
/// times in seconds. The curvature density is fixed by flatness:
/// Ω_k = 1 − Ω_m − Ω_Λ − Ω_r.
/// </summary>
public class FlrwModel
{
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double G = PhysicsConstants.GravitationalConstant;
    private const double Mpc = PhysicsConstants.Megaparsec;

    private const int IntegrationSteps = 2000;
    private const double ScaleFactorFloor = 1e-8; // lower bound for age integrals

    /// <summary>Hubble constant H₀ in s⁻¹.</summary>
    public double HubbleConstant { get; }

    /// <summary>Hubble constant H₀ in km·s⁻¹·Mpc⁻¹.</summary>
    public double HubbleConstantKmSMpc { get; }

    /// <summary>Present-day matter density parameter Ω_m.</summary>
    public double OmegaMatter { get; }

    /// <summary>Present-day dark-energy density parameter Ω_Λ.</summary>
    public double OmegaLambda { get; }

    /// <summary>Present-day radiation density parameter Ω_r.</summary>
    public double OmegaRadiation { get; }

    /// <summary>Curvature density parameter Ω_k = 1 − Ω_m − Ω_Λ − Ω_r.</summary>
    public double OmegaCurvature { get; }

    /// <summary>Hubble distance D_H = c/H₀ (metres).</summary>
    public double HubbleDistance => C / HubbleConstant;

    /// <summary>Hubble time 1/H₀ (seconds).</summary>
    public double HubbleTime => 1.0 / HubbleConstant;

    /// <summary>
    /// Creates an FLRW cosmology.
    /// </summary>
    /// <param name="hubbleConstantKmSMpc">H₀ in km·s⁻¹·Mpc⁻¹.</param>
    /// <param name="omegaMatter">Present-day matter density parameter Ω_m.</param>
    /// <param name="omegaLambda">Present-day dark-energy density parameter Ω_Λ.</param>
    /// <param name="omegaRadiation">Present-day radiation density parameter Ω_r (default 0).</param>
    public FlrwModel(double hubbleConstantKmSMpc, double omegaMatter, double omegaLambda, double omegaRadiation = 0.0)
    {
        if (hubbleConstantKmSMpc <= 0)
            throw new ArgumentOutOfRangeException(nameof(hubbleConstantKmSMpc), "H₀ must be positive.");

        HubbleConstantKmSMpc = hubbleConstantKmSMpc;
        HubbleConstant = hubbleConstantKmSMpc * 1000.0 / Mpc; // (km/s)/Mpc → 1/s
        OmegaMatter = omegaMatter;
        OmegaLambda = omegaLambda;
        OmegaRadiation = omegaRadiation;
        OmegaCurvature = 1.0 - omegaMatter - omegaLambda - omegaRadiation;
    }

    /// <summary>
    /// The Planck 2018 flat ΛCDM concordance cosmology
    /// (H₀ = 67.4, Ω_m = 0.315, Ω_Λ = 0.685).
    /// </summary>
    public static FlrwModel Planck2018() => new FlrwModel(67.4, 0.315, 0.685);

    /// <summary>
    /// Dimensionless Hubble function E(z) = H(z)/H₀
    /// = √(Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_k(1+z)² + Ω_Λ).
    /// </summary>
    public double DimensionlessHubble(double z)
    {
        double zp1 = 1.0 + z;
        double zp1_2 = zp1 * zp1;
        double zp1_3 = zp1_2 * zp1;
        double zp1_4 = zp1_3 * zp1;

        return Math.Sqrt(
            OmegaMatter * zp1_3 +
            OmegaRadiation * zp1_4 +
            OmegaCurvature * zp1_2 +
            OmegaLambda);
    }

    /// <summary>Hubble parameter H(z) = H₀·E(z) in s⁻¹.</summary>
    public double HubbleParameter(double z) => HubbleConstant * DimensionlessHubble(z);

    /// <summary>Scale factor a = 1/(1+z) (a = 1 today).</summary>
    public double ScaleFactor(double z) => 1.0 / (1.0 + z);

    /// <summary>
    /// Critical density ρ_c(z) = 3H(z)²/(8πG) in kg·m⁻³.
    /// With no argument, the present-day value ρ_c,0.
    /// </summary>
    public double CriticalDensity(double z = 0.0)
    {
        double h = HubbleParameter(z);
        return 3.0 * h * h / (8.0 * Math.PI * G);
    }

    /// <summary>
    /// Line-of-sight comoving distance D_C = D_H·∫₀ᶻ dz'/E(z') in metres.
    /// </summary>
    public double ComovingDistance(double z)
    {
        if (z < 0) throw new ArgumentOutOfRangeException(nameof(z), "Redshift must be non-negative.");
        if (z == 0) return 0.0;

        double integral = Simpson(zp => 1.0 / DimensionlessHubble(zp), 0.0, z, IntegrationSteps);
        return HubbleDistance * integral;
    }

    /// <summary>
    /// Transverse comoving distance D_M (metres), accounting for spatial curvature.
    /// Equals D_C for a flat universe; uses the sinh/sin forms otherwise.
    /// </summary>
    public double TransverseComovingDistance(double z)
    {
        double dc = ComovingDistance(z);
        double dh = HubbleDistance;

        if (Math.Abs(OmegaCurvature) < 1e-12)
            return dc;

        double sqrtOk = Math.Sqrt(Math.Abs(OmegaCurvature));
        if (OmegaCurvature > 0) // open
            return dh / sqrtOk * Math.Sinh(sqrtOk * dc / dh);

        return dh / sqrtOk * Math.Sin(sqrtOk * dc / dh); // closed
    }

    /// <summary>
    /// Luminosity distance D_L = (1+z)·D_M in metres.
    /// </summary>
    public double LuminosityDistance(double z) => (1.0 + z) * TransverseComovingDistance(z);

    /// <summary>
    /// Angular-diameter distance D_A = D_M/(1+z) in metres.
    /// </summary>
    public double AngularDiameterDistance(double z) => TransverseComovingDistance(z) / (1.0 + z);

    /// <summary>
    /// Distance modulus μ = 5·log₁₀(D_L / 10 pc) in magnitudes.
    /// </summary>
    public double DistanceModulus(double z)
    {
        double dLparsecs = LuminosityDistance(z) / PhysicsConstants.Parsec;
        return 5.0 * Math.Log10(dLparsecs) - 5.0;
    }

    /// <summary>
    /// Age of the universe at redshift <paramref name="z"/> in seconds:
    /// t(a) = (1/H₀)·∫₀ᵃ da'/(a'·E(a')), with a = 1/(1+z).
    /// </summary>
    public double Age(double z)
    {
        if (z < 0) throw new ArgumentOutOfRangeException(nameof(z), "Redshift must be non-negative.");
        return AgeIntegral(ScaleFactorFloor, ScaleFactor(z));
    }

    /// <summary>Present-day age of the universe t₀ in seconds.</summary>
    public double AgeOfUniverse() => AgeIntegral(ScaleFactorFloor, 1.0);

    /// <summary>
    /// Lookback time to redshift <paramref name="z"/> in seconds:
    /// the difference between the present age and the age at that redshift.
    /// </summary>
    public double LookbackTime(double z)
    {
        if (z < 0) throw new ArgumentOutOfRangeException(nameof(z), "Redshift must be non-negative.");
        return AgeIntegral(ScaleFactor(z), 1.0);
    }

    // t = (1/H₀)·∫ da/(a·E(a)) with E(a) = E(z) at a = 1/(1+z).
    private double AgeIntegral(double aLow, double aHigh)
    {
        if (aHigh <= aLow) return 0.0;

        double integral = Simpson(a =>
        {
            double zp = 1.0 / a - 1.0;
            return 1.0 / (a * DimensionlessHubble(zp));
        }, aLow, aHigh, IntegrationSteps);

        return integral / HubbleConstant;
    }

    // Composite Simpson's rule over [a, b] with an even number of intervals.
    private static double Simpson(Func<double, double> f, double a, double b, int n)
    {
        if (n % 2 != 0) n++;
        double h = (b - a) / n;
        double sum = f(a) + f(b);

        for (int i = 1; i < n; i++)
        {
            double x = a + i * h;
            sum += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
        }

        return sum * h / 3.0;
    }
}
