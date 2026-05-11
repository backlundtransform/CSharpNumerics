using System;

namespace CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

/// <summary>
/// Models engine thrust as a function of throttle setting, altitude, and airspeed.
/// Supports both jet (turbojet/turbofan) and propeller (piston/turboprop) engines.
/// 
/// Jet model:  T = T_max · throttle · (ρ / ρ₀)^n
/// Prop model: T = (P_max · throttle · η) / max(V, V_min), capped at T_static
/// </summary>
public class PropulsionModel
{
    /// <summary>Engine type.</summary>
    public PropulsionType Type { get; }

    /// <summary>Maximum thrust at sea level, full throttle (N) — for jets.</summary>
    public double MaxThrustSeaLevel { get; }

    /// <summary>Maximum shaft power (W) — for propeller engines.</summary>
    public double MaxPower { get; }

    /// <summary>Propeller efficiency (0–1) — for propeller engines.</summary>
    public double PropellerEfficiency { get; }

    /// <summary>Static thrust limit (N) — for propeller engines at V ≈ 0.</summary>
    public double StaticThrust { get; }

    /// <summary>Density exponent for altitude lapse (default 0.7 for turbofan, 1.0 for piston).</summary>
    public double DensityExponent { get; }

    /// <summary>Number of engines.</summary>
    public int EngineCount { get; }

    /// <summary>
    /// Creates a jet propulsion model.
    /// </summary>
    /// <param name="maxThrustSeaLevel">Max thrust per engine at sea level in Newtons.</param>
    /// <param name="densityExponent">Altitude lapse exponent (default 0.7).</param>
    /// <param name="engineCount">Number of engines (default 1).</param>
    public static PropulsionModel Jet(
        double maxThrustSeaLevel,
        double densityExponent = 0.7,
        int engineCount = 1)
    {
        return new PropulsionModel(
            PropulsionType.Jet,
            maxThrustSeaLevel,
            maxPower: 0,
            propellerEfficiency: 0,
            staticThrust: 0,
            densityExponent: densityExponent,
            engineCount: engineCount);
    }

    /// <summary>
    /// Creates a propeller propulsion model.
    /// </summary>
    /// <param name="maxPower">Maximum shaft power per engine in Watts.</param>
    /// <param name="propellerEfficiency">Propeller efficiency (default 0.85).</param>
    /// <param name="staticThrust">Static thrust limit per engine in Newtons.</param>
    /// <param name="densityExponent">Altitude lapse exponent (default 1.0 for naturally aspirated).</param>
    /// <param name="engineCount">Number of engines (default 1).</param>
    public static PropulsionModel Propeller(
        double maxPower,
        double propellerEfficiency = 0.85,
        double staticThrust = double.MaxValue,
        double densityExponent = 1.0,
        int engineCount = 1)
    {
        return new PropulsionModel(
            PropulsionType.Propeller,
            maxThrustSeaLevel: 0,
            maxPower: maxPower,
            propellerEfficiency: propellerEfficiency,
            staticThrust: staticThrust,
            densityExponent: densityExponent,
            engineCount: engineCount);
    }

    private PropulsionModel(
        PropulsionType type,
        double maxThrustSeaLevel,
        double maxPower,
        double propellerEfficiency,
        double staticThrust,
        double densityExponent,
        int engineCount)
    {
        Type = type;
        MaxThrustSeaLevel = maxThrustSeaLevel;
        MaxPower = maxPower;
        PropellerEfficiency = propellerEfficiency;
        StaticThrust = staticThrust;
        DensityExponent = densityExponent;
        EngineCount = Math.Max(1, engineCount);
    }

    /// <summary>
    /// Computes total thrust (all engines) at given conditions.
    /// </summary>
    /// <param name="throttle">Throttle setting (0–1).</param>
    /// <param name="altitude">Geopotential altitude in metres.</param>
    /// <param name="airspeed">True airspeed in m/s.</param>
    /// <returns>Total thrust in Newtons (positive = forward along body X).</returns>
    public double Thrust(double throttle, double altitude, double airspeed)
    {
        throttle = Math.Clamp(throttle, 0, 1);

        double rho = AtmosphereModel.Density(altitude);
        double rho0 = AtmosphereModel.SeaLevelDensity;
        double densityRatio = Math.Pow(rho / rho0, DensityExponent);

        double thrustPerEngine;

        if (Type == PropulsionType.Jet)
        {
            thrustPerEngine = MaxThrustSeaLevel * throttle * densityRatio;
        }
        else
        {
            // Prop: T = P·η / V, capped at static thrust
            double vMin = 1.0; // avoid division by zero
            double v = Math.Max(airspeed, vMin);
            thrustPerEngine = (MaxPower * throttle * PropellerEfficiency * densityRatio) / v;
            thrustPerEngine = Math.Min(thrustPerEngine, StaticThrust * throttle * densityRatio);
        }

        return thrustPerEngine * EngineCount;
    }
}

/// <summary>
/// Engine type for propulsion models.
/// </summary>
public enum PropulsionType
{
    /// <summary>Turbojet or turbofan engine.</summary>
    Jet,
    /// <summary>Piston or turboprop with propeller.</summary>
    Propeller
}
