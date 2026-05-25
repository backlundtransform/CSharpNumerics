using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics.Propulsion;
using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Represents a single rocket stage with dry mass, engines, propellant tanks,
/// and separation trigger conditions.
/// </summary>
public class RocketStage
{
    /// <summary>Dry mass of the stage structure (no propellant) in kg.</summary>
    public double DryMass { get; }

    /// <summary>Engines mounted on this stage.</summary>
    public RocketEngine[] Engines { get; }

    /// <summary>Propellant tanks for this stage.</summary>
    public PropellantTank[] Tanks { get; }

    /// <summary>Center of mass tracker for this stage.</summary>
    public CenterOfMassTracker CgTracker { get; }

    /// <summary>Separation trigger condition.</summary>
    public StageSeparationTrigger SeparationTrigger { get; set; }

    /// <summary>Total current propellant mass across all tanks (kg).</summary>
    public double PropellantMass
    {
        get
        {
            double total = 0;
            for (int i = 0; i < Tanks.Length; i++)
                total += Tanks[i].PropellantMass;
            return total;
        }
    }

    /// <summary>Total initial propellant mass across all tanks (kg).</summary>
    public double InitialPropellantMass
    {
        get
        {
            double total = 0;
            for (int i = 0; i < Tanks.Length; i++)
                total += Tanks[i].InitialPropellantMass;
            return total;
        }
    }

    /// <summary>Total wet mass (dry + propellant) in kg.</summary>
    public double TotalMass => DryMass + PropellantMass;

    /// <summary>True when all tanks are empty.</summary>
    public bool IsDepleted
    {
        get
        {
            for (int i = 0; i < Tanks.Length; i++)
                if (!Tanks[i].IsEmpty) return false;
            return true;
        }
    }

    /// <summary>
    /// Creates a rocket stage.
    /// </summary>
    /// <param name="dryMass">Structural dry mass in kg.</param>
    /// <param name="engines">Array of rocket engines.</param>
    /// <param name="tanks">Array of propellant tanks.</param>
    /// <param name="cgTracker">CG tracker (optional, auto-created if null).</param>
    public RocketStage(double dryMass, RocketEngine[] engines, PropellantTank[] tanks, CenterOfMassTracker cgTracker = null)
    {
        if (dryMass <= 0) throw new ArgumentOutOfRangeException(nameof(dryMass));
        if (engines == null || engines.Length == 0) throw new ArgumentException("At least one engine required.", nameof(engines));
        if (tanks == null || tanks.Length == 0) throw new ArgumentException("At least one tank required.", nameof(tanks));

        DryMass = dryMass;
        Engines = engines;
        Tanks = tanks;

        double totalPropellant = 0;
        for (int i = 0; i < tanks.Length; i++)
            totalPropellant += tanks[i].InitialPropellantMass;

        CgTracker = cgTracker ?? new CenterOfMassTracker(dryMass, totalPropellant, cgFull: 0, cgEmpty: 0);
        SeparationTrigger = StageSeparationTrigger.FuelDepletion();
    }

    /// <summary>
    /// Computes total thrust from all active engines at the given pressure ratio.
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure / sea-level pressure (0–1).</param>
    public double TotalThrust(double pressureRatio)
    {
        double total = 0;
        for (int i = 0; i < Engines.Length; i++)
            total += Engines[i].Thrust(pressureRatio);
        return total;
    }

    /// <summary>
    /// Computes total mass flow rate from all active engines.
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure / sea-level pressure (0–1).</param>
    public double TotalMassFlowRate(double pressureRatio)
    {
        double total = 0;
        for (int i = 0; i < Engines.Length; i++)
            total += Engines[i].MassFlowRate(pressureRatio);
        return total;
    }

    /// <summary>
    /// Consumes propellant from tanks for the given time step.
    /// Distributes consumption evenly across tanks.
    /// </summary>
    /// <param name="dt">Time step in seconds.</param>
    /// <param name="pressureRatio">Ambient pressure ratio for mass flow calculation.</param>
    /// <returns>Total mass consumed in kg.</returns>
    public double ConsumePropellant(double dt, double pressureRatio)
    {
        double totalFlow = TotalMassFlowRate(pressureRatio);
        if (totalFlow <= 0) return 0;

        double flowPerTank = totalFlow / Tanks.Length;
        double totalConsumed = 0;
        for (int i = 0; i < Tanks.Length; i++)
            totalConsumed += Tanks[i].Consume(flowPerTank, dt);
        return totalConsumed;
    }

    /// <summary>
    /// Activates all engines on this stage.
    /// </summary>
    public void Ignite()
    {
        for (int i = 0; i < Engines.Length; i++)
            Engines[i].IsActive = true;
    }

    /// <summary>
    /// Shuts down all engines on this stage.
    /// </summary>
    public void Shutdown()
    {
        for (int i = 0; i < Engines.Length; i++)
            Engines[i].IsActive = false;
    }
}

/// <summary>
/// Defines the condition that triggers stage separation.
/// </summary>
public class StageSeparationTrigger
{
    /// <summary>Type of trigger condition.</summary>
    public SeparationType Type { get; }

    /// <summary>Value for the trigger (time in seconds, or altitude in metres).</summary>
    public double Value { get; }

    private StageSeparationTrigger(SeparationType type, double value = 0)
    {
        Type = type;
        Value = value;
    }

    /// <summary>Separate when fuel is depleted.</summary>
    public static StageSeparationTrigger FuelDepletion() => new(SeparationType.FuelDepletion);

    /// <summary>Separate at a specific time after launch.</summary>
    public static StageSeparationTrigger AtTime(double time) => new(SeparationType.Time, time);

    /// <summary>Separate at a specific altitude.</summary>
    public static StageSeparationTrigger AtAltitude(double altitude) => new(SeparationType.Altitude, altitude);
}

/// <summary>
/// Stage separation trigger type.
/// </summary>
public enum SeparationType
{
    /// <summary>Separate when propellant is exhausted.</summary>
    FuelDepletion,
    /// <summary>Separate at a specific time.</summary>
    Time,
    /// <summary>Separate at a specific altitude.</summary>
    Altitude
}
