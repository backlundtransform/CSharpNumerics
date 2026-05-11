using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// International Standard Atmosphere (ISA) model.
/// Computes density, pressure, temperature, and speed of sound as functions of altitude
/// from sea level to 86 km using the standard lapse-rate layers.
/// 
/// Reference: ICAO Doc 7488, ISO 2533:1975.
/// 
/// Layer boundaries (geopotential altitude):
///   0–11 km    Troposphere      lapse −6.5 K/km
///  11–20 km    Tropopause       isothermal
///  20–32 km    Stratosphere 1   lapse +1.0 K/km
///  32–47 km    Stratosphere 2   lapse +2.8 K/km
///  47–51 km    Stratopause      isothermal
///  51–71 km    Mesosphere 1     lapse −2.8 K/km
///  71–86 km    Mesosphere 2     lapse −2.0 K/km
/// </summary>
public static class AtmosphereModel
{
    // Sea-level reference values
    private const double T0 = 288.15;       // K
    private const double P0 = 101325.0;     // Pa
    private const double Rho0 = 1.225;      // kg/m³
    private const double R = 287.05287;     // J/(kg·K) — specific gas constant for dry air
    private const double Gamma = 1.4;       // ratio of specific heats
    private const double G0 = 9.80665;      // m/s²

    // Layer base altitudes (m) and lapse rates (K/m)
    private static readonly double[] BaseAlt = { 0, 11000, 20000, 32000, 47000, 51000, 71000 };
    private static readonly double[] LapseRate = { -0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002 };

    // Precomputed base temperatures and pressures at each layer boundary
    private static readonly double[] BaseTemp;
    private static readonly double[] BasePressure;

    static AtmosphereModel()
    {
        BaseTemp = new double[BaseAlt.Length];
        BasePressure = new double[BaseAlt.Length];
        BaseTemp[0] = T0;
        BasePressure[0] = P0;

        for (int i = 1; i < BaseAlt.Length; i++)
        {
            double dh = BaseAlt[i] - BaseAlt[i - 1];
            double lr = LapseRate[i - 1];
            BaseTemp[i] = BaseTemp[i - 1] + lr * dh;

            if (Math.Abs(lr) < 1e-10)
            {
                // Isothermal layer
                BasePressure[i] = BasePressure[i - 1] * Math.Exp(-G0 * dh / (R * BaseTemp[i - 1]));
            }
            else
            {
                // Gradient layer
                BasePressure[i] = BasePressure[i - 1] *
                    Math.Pow(BaseTemp[i] / BaseTemp[i - 1], -G0 / (lr * R));
            }
        }
    }

    /// <summary>
    /// Returns the ISA layer index for the given geopotential altitude.
    /// </summary>
    private static int FindLayer(double altitude)
    {
        for (int i = BaseAlt.Length - 1; i >= 0; i--)
        {
            if (altitude >= BaseAlt[i])
                return i;
        }
        return 0;
    }

    /// <summary>
    /// Temperature at altitude in Kelvin.
    /// </summary>
    /// <param name="altitude">Geopotential altitude in metres (0–86 000 m).</param>
    public static double Temperature(double altitude)
    {
        altitude = Math.Clamp(altitude, 0, 86000);
        int i = FindLayer(altitude);
        return BaseTemp[i] + LapseRate[i] * (altitude - BaseAlt[i]);
    }

    /// <summary>
    /// Static pressure at altitude in Pascals.
    /// </summary>
    /// <param name="altitude">Geopotential altitude in metres.</param>
    public static double Pressure(double altitude)
    {
        altitude = Math.Clamp(altitude, 0, 86000);
        int i = FindLayer(altitude);
        double dh = altitude - BaseAlt[i];
        double lr = LapseRate[i];

        if (Math.Abs(lr) < 1e-10)
            return BasePressure[i] * Math.Exp(-G0 * dh / (R * BaseTemp[i]));

        double T = BaseTemp[i] + lr * dh;
        return BasePressure[i] * Math.Pow(T / BaseTemp[i], -G0 / (lr * R));
    }

    /// <summary>
    /// Air density at altitude in kg/m³ (from ideal gas law: ρ = p / (R·T)).
    /// </summary>
    /// <param name="altitude">Geopotential altitude in metres.</param>
    public static double Density(double altitude)
    {
        return Pressure(altitude) / (R * Temperature(altitude));
    }

    /// <summary>
    /// Speed of sound at altitude in m/s: a = √(γ·R·T).
    /// </summary>
    /// <param name="altitude">Geopotential altitude in metres.</param>
    public static double SpeedOfSound(double altitude)
    {
        return Math.Sqrt(Gamma * R * Temperature(altitude));
    }

    /// <summary>
    /// Dynamic viscosity at altitude in Pa·s using Sutherland's law.
    /// μ = μ_ref · (T/T_ref)^(3/2) · (T_ref + S) / (T + S)
    /// </summary>
    /// <param name="altitude">Geopotential altitude in metres.</param>
    public static double DynamicViscosity(double altitude)
    {
        const double muRef = 1.716e-5;  // Pa·s at T_ref
        const double tRef = 273.15;     // K
        const double S = 110.4;         // Sutherland constant K

        double T = Temperature(altitude);
        return muRef * Math.Pow(T / tRef, 1.5) * (tRef + S) / (T + S);
    }

    /// <summary>
    /// Kinematic viscosity at altitude in m²/s: ν = μ / ρ.
    /// </summary>
    /// <param name="altitude">Geopotential altitude in metres.</param>
    public static double KinematicViscosity(double altitude)
    {
        return DynamicViscosity(altitude) / Density(altitude);
    }

    /// <summary>
    /// Sea-level temperature in Kelvin (288.15 K).
    /// </summary>
    public static double SeaLevelTemperature => T0;

    /// <summary>
    /// Sea-level pressure in Pascals (101 325 Pa).
    /// </summary>
    public static double SeaLevelPressure => P0;

    /// <summary>
    /// Sea-level density in kg/m³ (1.225).
    /// </summary>
    public static double SeaLevelDensity => Rho0;
}
