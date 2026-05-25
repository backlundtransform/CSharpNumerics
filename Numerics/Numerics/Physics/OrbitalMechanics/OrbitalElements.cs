using System;

namespace CSharpNumerics.Physics.OrbitalMechanics;

/// <summary>
/// Classical Keplerian orbital elements defining an orbit in 3D space.
/// </summary>
public struct OrbitalElements
{
    /// <summary>Semi-major axis in meters. For hyperbolic orbits, a is negative.</summary>
    public double SemiMajorAxis;

    /// <summary>Eccentricity. 0=circle, 0&lt;e&lt;1=ellipse, e=1=parabola, e&gt;1=hyperbola.</summary>
    public double Eccentricity;

    /// <summary>Inclination in radians [0, π].</summary>
    public double Inclination;

    /// <summary>Right Ascension of Ascending Node (RAAN/Ω) in radians [0, 2π).</summary>
    public double RAAN;

    /// <summary>Argument of periapsis (ω) in radians [0, 2π).</summary>
    public double ArgumentOfPeriapsis;

    /// <summary>True anomaly (ν) in radians [0, 2π).</summary>
    public double TrueAnomaly;

    /// <summary>Gravitational parameter μ = GM used for this orbit (m³/s²).</summary>
    public double Mu;

    public OrbitalElements(double semiMajorAxis, double eccentricity, double inclination,
        double raan, double argumentOfPeriapsis, double trueAnomaly, double mu)
    {
        SemiMajorAxis = semiMajorAxis;
        Eccentricity = eccentricity;
        Inclination = inclination;
        RAAN = raan;
        ArgumentOfPeriapsis = argumentOfPeriapsis;
        TrueAnomaly = trueAnomaly;
        Mu = mu;
    }

    /// <summary>Orbital period in seconds (only valid for elliptical orbits, e &lt; 1).</summary>
    public double Period => 2.0 * Math.PI * Math.Sqrt(SemiMajorAxis * SemiMajorAxis * SemiMajorAxis / Mu);

    /// <summary>Mean motion n = 2π/T in rad/s.</summary>
    public double MeanMotion => Math.Sqrt(Mu / (SemiMajorAxis * SemiMajorAxis * SemiMajorAxis));

    /// <summary>Periapsis radius in meters.</summary>
    public double Periapsis => SemiMajorAxis * (1.0 - Eccentricity);

    /// <summary>Apoapsis radius in meters (only valid for elliptical orbits).</summary>
    public double Apoapsis => SemiMajorAxis * (1.0 + Eccentricity);

    /// <summary>Periapsis altitude above Earth's surface (using mean radius 6371 km).</summary>
    public double PeriapsisAltitude => Periapsis - 6371000.0;

    /// <summary>Apoapsis altitude above Earth's surface (using mean radius 6371 km).</summary>
    public double ApoapsisAltitude => Apoapsis - 6371000.0;

    /// <summary>Specific orbital energy ε = -μ/(2a).</summary>
    public double SpecificEnergy => -Mu / (2.0 * SemiMajorAxis);

    /// <summary>Specific angular momentum magnitude h = sqrt(μ·a·(1-e²)).</summary>
    public double SpecificAngularMomentum =>
        Math.Sqrt(Mu * SemiMajorAxis * (1.0 - Eccentricity * Eccentricity));

    /// <summary>Velocity at current true anomaly.</summary>
    public double VelocityAtTrueAnomaly
    {
        get
        {
            double r = SemiMajorAxis * (1.0 - Eccentricity * Eccentricity) /
                       (1.0 + Eccentricity * Math.Cos(TrueAnomaly));
            return Math.Sqrt(Mu * (2.0 / r - 1.0 / SemiMajorAxis));
        }
    }

    /// <summary>Radius at current true anomaly.</summary>
    public double RadiusAtTrueAnomaly =>
        SemiMajorAxis * (1.0 - Eccentricity * Eccentricity) /
        (1.0 + Eccentricity * Math.Cos(TrueAnomaly));

    /// <summary>Whether this orbit is bound (elliptical).</summary>
    public bool IsBound => Eccentricity < 1.0;

    /// <summary>Whether this is a circular orbit (e ≈ 0).</summary>
    public bool IsCircular => Eccentricity < 1e-6;
}
