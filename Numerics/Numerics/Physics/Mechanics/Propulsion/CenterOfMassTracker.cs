using System;

namespace CSharpNumerics.Physics.Mechanics.Propulsion;

/// <summary>
/// Tracks center of mass position as propellant is consumed.
/// Models CG shift along the longitudinal axis of a rocket stage.
/// </summary>
public class CenterOfMassTracker
{
    /// <summary>CG position when tanks are full (distance from reference point along body X-axis, metres).</summary>
    public double CgFull { get; }

    /// <summary>CG position when tanks are empty (distance from reference point along body X-axis, metres).</summary>
    public double CgEmpty { get; }

    /// <summary>Dry mass of the structure (kg).</summary>
    public double DryMass { get; }

    /// <summary>Maximum propellant mass (kg).</summary>
    public double MaxPropellantMass { get; }

    /// <summary>
    /// Creates a CG tracker for a rocket stage.
    /// </summary>
    /// <param name="dryMass">Structural dry mass in kg.</param>
    /// <param name="maxPropellantMass">Maximum propellant mass in kg.</param>
    /// <param name="cgFull">CG position at full load (m from reference, along body X).</param>
    /// <param name="cgEmpty">CG position at empty tanks (m from reference, along body X).</param>
    public CenterOfMassTracker(double dryMass, double maxPropellantMass, double cgFull, double cgEmpty)
    {
        if (dryMass <= 0) throw new ArgumentOutOfRangeException(nameof(dryMass));
        if (maxPropellantMass < 0) throw new ArgumentOutOfRangeException(nameof(maxPropellantMass));

        DryMass = dryMass;
        MaxPropellantMass = maxPropellantMass;
        CgFull = cgFull;
        CgEmpty = cgEmpty;
    }

    /// <summary>
    /// Computes current CG position based on remaining propellant mass.
    /// Linearly interpolates between full and empty CG positions.
    /// </summary>
    /// <param name="currentPropellantMass">Current propellant mass in kg.</param>
    /// <returns>CG position in metres from reference point along body X-axis.</returns>
    public double ComputeCg(double currentPropellantMass)
    {
        double fraction = MaxPropellantMass > 0
            ? Math.Clamp(currentPropellantMass / MaxPropellantMass, 0, 1)
            : 0;

        return CgEmpty + fraction * (CgFull - CgEmpty);
    }

    /// <summary>
    /// Computes the vector from CG to the engine mount point (aft of vehicle).
    /// This is the moment arm for thrust vectoring torque calculations.
    /// </summary>
    /// <param name="currentPropellantMass">Current propellant mass in kg.</param>
    /// <param name="enginePosition">Engine mount position (m from reference along body X).</param>
    /// <returns>Distance from CG to engine (positive = engine is aft of CG).</returns>
    public double CgToEngineDistance(double currentPropellantMass, double enginePosition)
    {
        double cg = ComputeCg(currentPropellantMass);
        return enginePosition - cg;
    }
}
