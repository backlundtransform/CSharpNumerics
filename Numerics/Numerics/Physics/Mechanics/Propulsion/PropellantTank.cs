using System;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Models a propellant tank with fuel mass, oxidizer mass, and mixture ratio.
/// Tracks consumption and computes remaining ΔV.
/// </summary>
public class PropellantTank
{
    private const double G0 = 9.80665;

    /// <summary>Initial fuel mass in kg.</summary>
    public double InitialFuelMass { get; }

    /// <summary>Initial oxidizer mass in kg.</summary>
    public double InitialOxidizerMass { get; }

    /// <summary>Oxidizer-to-fuel mixture ratio by mass (e.g. 2.56 for RP-1/LOX).</summary>
    public double MixtureRatio { get; }

    /// <summary>Current remaining fuel mass in kg.</summary>
    public double FuelMass { get; private set; }

    /// <summary>Current remaining oxidizer mass in kg.</summary>
    public double OxidizerMass { get; private set; }

    /// <summary>Total initial propellant mass (fuel + oxidizer) in kg.</summary>
    public double InitialPropellantMass => InitialFuelMass + InitialOxidizerMass;

    /// <summary>Total remaining propellant mass (fuel + oxidizer) in kg.</summary>
    public double PropellantMass => FuelMass + OxidizerMass;

    /// <summary>True when all propellant is exhausted.</summary>
    public bool IsEmpty => PropellantMass <= 1e-6;

    /// <summary>Propellant mass fraction consumed (0 = full, 1 = empty).</summary>
    public double ConsumedFraction => 1.0 - PropellantMass / InitialPropellantMass;

    /// <summary>
    /// Creates a propellant tank.
    /// </summary>
    /// <param name="fuelMass">Fuel mass in kg.</param>
    /// <param name="oxidizerMass">Oxidizer mass in kg.</param>
    public PropellantTank(double fuelMass, double oxidizerMass)
    {
        if (fuelMass < 0) throw new ArgumentOutOfRangeException(nameof(fuelMass));
        if (oxidizerMass < 0) throw new ArgumentOutOfRangeException(nameof(oxidizerMass));

        InitialFuelMass = fuelMass;
        InitialOxidizerMass = oxidizerMass;
        FuelMass = fuelMass;
        OxidizerMass = oxidizerMass;
        MixtureRatio = fuelMass > 0 ? oxidizerMass / fuelMass : 0;
    }

    /// <summary>
    /// Consumes propellant at the given total mass flow rate for a time step.
    /// Splits consumption between fuel and oxidizer by mixture ratio.
    /// Returns actual mass consumed (may be less than requested if tank runs dry).
    /// </summary>
    /// <param name="massFlowRate">Total mass flow rate in kg/s (fuel + oxidizer).</param>
    /// <param name="dt">Time step in seconds.</param>
    /// <returns>Actual propellant mass consumed in kg.</returns>
    public double Consume(double massFlowRate, double dt)
    {
        if (massFlowRate <= 0 || dt <= 0 || IsEmpty) return 0;

        double requested = massFlowRate * dt;
        double available = PropellantMass;
        double consumed = Math.Min(requested, available);

        // Split by mixture ratio: O/F = MixtureRatio
        // Total = fuel + oxidizer = fuel · (1 + MixtureRatio)
        double fuelFraction = 1.0 / (1.0 + MixtureRatio);
        double fuelConsumed = consumed * fuelFraction;
        double oxConsumed = consumed - fuelConsumed;

        FuelMass = Math.Max(0, FuelMass - fuelConsumed);
        OxidizerMass = Math.Max(0, OxidizerMass - oxConsumed);

        return consumed;
    }

    /// <summary>
    /// Computes remaining ΔV given the current propellant, dry mass, and engine Isp.
    /// ΔV = Isp · g₀ · ln((m_dry + m_prop) / m_dry)
    /// </summary>
    /// <param name="dryMass">Dry mass of the vehicle/stage in kg.</param>
    /// <param name="isp">Effective specific impulse in seconds.</param>
    public double RemainingDeltaV(double dryMass, double isp)
    {
        if (dryMass <= 0 || isp <= 0 || IsEmpty) return 0;
        double m0 = dryMass + PropellantMass;
        double mf = dryMass;
        return isp * G0 * Math.Log(m0 / mf);
    }

    /// <summary>
    /// Resets the tank to its initial propellant load.
    /// </summary>
    public void Reset()
    {
        FuelMass = InitialFuelMass;
        OxidizerMass = InitialOxidizerMass;
    }
}
