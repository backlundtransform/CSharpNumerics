using CSharpNumerics.Physics.Mechanics.Propulsion;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Multi-stage rocket vehicle: a stack of <see cref="RocketStage"/> objects
/// with active stage tracking, total mass computation, and composite CG.
/// Supports parallel boosters (strap-ons) that fire alongside the core stage.
/// </summary>
public class RocketVehicle
{
    /// <summary>All stages from bottom (index 0) to top (payload).</summary>
    public RocketStage[] Stages { get; }

    /// <summary>Index of the currently active (firing) stage. Increments on separation.</summary>
    public int ActiveStageIndex { get; private set; }

    /// <summary>Payload mass at the top of the stack (kg).</summary>
    public double PayloadMass { get; }

    /// <summary>Strap-on boosters that fire in parallel with the core first stage.</summary>
    public List<RocketStage> Boosters { get; } = new List<RocketStage>();

    /// <summary>True if boosters have been separated.</summary>
    public bool BoostersSeparated { get; private set; }

    /// <summary>The currently active stage, or null if all stages are spent.</summary>
    public RocketStage ActiveStage => ActiveStageIndex < Stages.Length ? Stages[ActiveStageIndex] : null;

    /// <summary>True if all stages have been separated/consumed.</summary>
    public bool AllStagesSpent => ActiveStageIndex >= Stages.Length;

    /// <summary>
    /// Total mass of the vehicle: sum of all remaining stages (from active upward) + payload + boosters.
    /// </summary>
    public double TotalMass
    {
        get
        {
            double mass = PayloadMass;
            for (int i = ActiveStageIndex; i < Stages.Length; i++)
                mass += Stages[i].TotalMass;
            if (!BoostersSeparated)
            {
                for (int i = 0; i < Boosters.Count; i++)
                    mass += Boosters[i].TotalMass;
            }
            return mass;
        }
    }

    /// <summary>
    /// Total dry mass (no propellant) of all remaining stages + payload.
    /// </summary>
    public double DryMass
    {
        get
        {
            double mass = PayloadMass;
            for (int i = ActiveStageIndex; i < Stages.Length; i++)
                mass += Stages[i].DryMass;
            return mass;
        }
    }

    /// <summary>
    /// Creates a rocket vehicle.
    /// Stages are ordered bottom-first: index 0 is the first stage to fire.
    /// </summary>
    /// <param name="stages">Array of stages (bottom-first order).</param>
    /// <param name="payloadMass">Payload mass in kg (default 0).</param>
    public RocketVehicle(RocketStage[] stages, double payloadMass = 0)
    {
        if (stages == null || stages.Length == 0)
            throw new ArgumentException("At least one stage required.", nameof(stages));

        Stages = stages;
        PayloadMass = payloadMass;
        ActiveStageIndex = 0;
    }

    /// <summary>
    /// Separates the current active stage and advances to the next.
    /// Returns the separated stage, or null if no stages remain.
    /// </summary>
    public RocketStage Separate()
    {
        if (AllStagesSpent) return null;

        var separated = Stages[ActiveStageIndex];
        separated.Shutdown();
        ActiveStageIndex++;

        // Ignite next stage if available
        if (!AllStagesSpent)
            Stages[ActiveStageIndex].Ignite();

        return separated;
    }

    /// <summary>
    /// Checks if the current active stage should be separated based on its trigger.
    /// </summary>
    /// <param name="time">Current simulation time in seconds.</param>
    /// <param name="altitude">Current altitude in metres.</param>
    public bool ShouldSeparate(double time, double altitude)
    {
        if (AllStagesSpent) return false;
        var stage = ActiveStage;
        var trigger = stage.SeparationTrigger;

        return trigger.Type switch
        {
            SeparationType.FuelDepletion => stage.IsDepleted,
            SeparationType.Time => time >= trigger.Value,
            SeparationType.Altitude => altitude >= trigger.Value,
            _ => false
        };
    }

    /// <summary>
    /// Computes composite CG position of the remaining stack.
    /// Returns the mass-weighted average CG along the longitudinal axis.
    /// </summary>
    /// <param name="stageOffsets">Longitudinal offset of each stage's reference point (m).</param>
    public double CompositeCg(double[] stageOffsets)
    {
        double totalMoment = 0;
        double totalMass = PayloadMass;

        for (int i = ActiveStageIndex; i < Stages.Length; i++)
        {
            var stage = Stages[i];
            double stageCg = stage.CgTracker.ComputeCg(stage.PropellantMass);
            double absoluteCg = stageOffsets[i] + stageCg;
            double stageMass = stage.TotalMass;
            totalMoment += absoluteCg * stageMass;
            totalMass += stageMass;
        }

        return totalMass > 0 ? totalMoment / totalMass : 0;
    }

    /// <summary>
    /// Ignites the first stage to begin the launch sequence.
    /// Also ignites any attached boosters.
    /// </summary>
    public void Ignite()
    {
        if (!AllStagesSpent)
            ActiveStage.Ignite();
        for (int i = 0; i < Boosters.Count; i++)
            Boosters[i].Ignite();
    }

    /// <summary>
    /// Checks if boosters should be separated (all depleted).
    /// </summary>
    public bool ShouldSeparateBoosters()
    {
        if (BoostersSeparated || Boosters.Count == 0) return false;
        for (int i = 0; i < Boosters.Count; i++)
            if (!Boosters[i].IsDepleted) return false;
        return true;
    }

    /// <summary>
    /// Separates all boosters simultaneously.
    /// Returns the total mass of separated boosters.
    /// </summary>
    public double SeparateBoosters()
    {
        if (BoostersSeparated || Boosters.Count == 0) return 0;

        double mass = 0;
        for (int i = 0; i < Boosters.Count; i++)
        {
            Boosters[i].Shutdown();
            mass += Boosters[i].TotalMass;
        }
        BoostersSeparated = true;
        return mass;
    }

    /// <summary>
    /// Gets total thrust from boosters (if still attached and active).
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure ratio.</param>
    public double BoosterThrust(double pressureRatio)
    {
        if (BoostersSeparated) return 0;
        double total = 0;
        for (int i = 0; i < Boosters.Count; i++)
            total += Boosters[i].TotalThrust(pressureRatio);
        return total;
    }

    /// <summary>
    /// Gets total mass flow rate from boosters.
    /// </summary>
    /// <param name="pressureRatio">Ambient pressure ratio.</param>
    public double BoosterMassFlowRate(double pressureRatio)
    {
        if (BoostersSeparated) return 0;
        double total = 0;
        for (int i = 0; i < Boosters.Count; i++)
            total += Boosters[i].TotalMassFlowRate(pressureRatio);
        return total;
    }

    /// <summary>
    /// Consumes propellant from all active boosters.
    /// </summary>
    public void ConsumeBoosterPropellant(double dt, double pressureRatio)
    {
        if (BoostersSeparated) return;
        for (int i = 0; i < Boosters.Count; i++)
            Boosters[i].ConsumePropellant(dt, pressureRatio);
    }

    /// <summary>
    /// Resets vehicle to initial state (all stages active, tanks full).
    /// </summary>
    public void Reset()
    {
        ActiveStageIndex = 0;
        BoostersSeparated = false;
        for (int i = 0; i < Stages.Length; i++)
        {
            Stages[i].Shutdown();
            for (int j = 0; j < Stages[i].Tanks.Length; j++)
                Stages[i].Tanks[j].Reset();
        }
        for (int i = 0; i < Boosters.Count; i++)
        {
            Boosters[i].Shutdown();
            for (int j = 0; j < Boosters[i].Tanks.Length; j++)
                Boosters[i].Tanks[j].Reset();
        }
    }
}
