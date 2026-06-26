using System;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Relativity;

/// <summary>
/// Weak-field gravitational-wave relations for a compact binary on a circular
/// orbit (the quadrupole approximation and the Peters inspiral formula):
/// chirp mass, gravitational-wave luminosity and frequency, time to coalescence,
/// and the strain amplitude. Masses are in kg, lengths in metres, time in seconds.
/// </summary>
public static class GravitationalWaves
{
    private const double G = PhysicsConstants.GravitationalConstant;
    private const double C = PhysicsConstants.SpeedOfLight;
    private const double C5 = C * C * C * C * C;

    /// <summary>
    /// Chirp mass M_c = (m₁·m₂)^(3/5) / (m₁ + m₂)^(1/5) — the mass combination
    /// that governs the inspiral waveform.
    /// </summary>
    public static double ChirpMass(double mass1, double mass2)
    {
        if (mass1 <= 0 || mass2 <= 0) throw new ArgumentOutOfRangeException(nameof(mass1), "Masses must be positive.");
        return Math.Pow(mass1 * mass2, 0.6) / Math.Pow(mass1 + mass2, 0.2);
    }

    /// <summary>
    /// Orbital frequency of the binary, f_orb = (1/2π)·√(G(m₁+m₂)/a³) (Hz).
    /// </summary>
    public static double OrbitalFrequency(double mass1, double mass2, double separation)
    {
        if (separation <= 0) throw new ArgumentOutOfRangeException(nameof(separation), "Separation must be positive.");
        return Math.Sqrt(G * (mass1 + mass2) / (separation * separation * separation)) / (2.0 * Math.PI);
    }

    /// <summary>
    /// Dominant gravitational-wave frequency f_gw = 2·f_orb (Hz).
    /// </summary>
    public static double GravitationalWaveFrequency(double mass1, double mass2, double separation)
        => 2.0 * OrbitalFrequency(mass1, mass2, separation);

    /// <summary>
    /// Gravitational-wave luminosity radiated by a circular binary
    /// (quadrupole formula): L = (32/5)·G⁴/c⁵·(m₁m₂)²(m₁+m₂)/a⁵ (watts).
    /// </summary>
    public static double Luminosity(double mass1, double mass2, double separation)
    {
        if (separation <= 0) throw new ArgumentOutOfRangeException(nameof(separation), "Separation must be positive.");

        double a5 = Math.Pow(separation, 5);
        double m1m2 = mass1 * mass2;
        return 32.0 / 5.0 * Math.Pow(G, 4) / C5 * m1m2 * m1m2 * (mass1 + mass2) / a5;
    }

    /// <summary>
    /// Time to coalescence from an initial circular separation (Peters 1964):
    /// t = (5/256)·c⁵·a₀⁴ / (G³·m₁·m₂·(m₁+m₂)) (seconds).
    /// </summary>
    public static double MergerTime(double mass1, double mass2, double initialSeparation)
    {
        if (initialSeparation <= 0) throw new ArgumentOutOfRangeException(nameof(initialSeparation), "Separation must be positive.");
        if (mass1 <= 0 || mass2 <= 0) throw new ArgumentOutOfRangeException(nameof(mass1), "Masses must be positive.");

        double a4 = Math.Pow(initialSeparation, 4);
        double g3 = G * G * G;
        return 5.0 / 256.0 * C5 * a4 / (g3 * mass1 * mass2 * (mass1 + mass2));
    }

    /// <summary>
    /// Order-of-magnitude gravitational-wave strain amplitude at a distance D:
    /// h ≈ (4/D)·(G·M_c/c²)^(5/3)·(π·f_gw/c)^(2/3) (dimensionless).
    /// </summary>
    /// <param name="chirpMass">Chirp mass M_c (kg).</param>
    /// <param name="distance">Distance to the source D (m).</param>
    /// <param name="gravitationalWaveFrequency">Gravitational-wave frequency f_gw (Hz).</param>
    public static double StrainAmplitude(double chirpMass, double distance, double gravitationalWaveFrequency)
    {
        if (distance <= 0) throw new ArgumentOutOfRangeException(nameof(distance), "Distance must be positive.");

        double geometricMass = G * chirpMass / (C * C);                    // metres
        double term1 = Math.Pow(geometricMass, 5.0 / 3.0);
        double term2 = Math.Pow(Math.PI * gravitationalWaveFrequency / C, 2.0 / 3.0);
        return 4.0 / distance * term1 * term2;
    }
}
