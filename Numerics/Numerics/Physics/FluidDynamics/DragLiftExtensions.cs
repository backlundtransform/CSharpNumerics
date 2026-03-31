using CSharpNumerics.Physics.Constants;
using System;

namespace CSharpNumerics.Physics.FluidDynamics;

/// <summary>
/// Extension methods for drag, lift, and terminal velocity calculations.
/// </summary>
public static class DragLiftExtensions
{
    // ═══════════════════════════════════════════════════════════════
    //  Drag & lift
    // ═══════════════════════════════════════════════════════════════

    #region Drag & Lift

    /// <summary>
    /// Computes the drag force on a body: F_D = ½ρv²C_D A.
    /// </summary>
    /// <param name="dragCoefficient">Drag coefficient C_D (dimensionless).</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="speed">Flow speed relative to the body in m/s.</param>
    /// <param name="referenceArea">Reference (frontal) area A in m².</param>
    public static double DragForce(
        this double dragCoefficient, double density, double speed, double referenceArea)
    {
        return 0.5 * density * speed * speed * dragCoefficient * referenceArea;
    }

    /// <summary>
    /// Computes the lift force on a body: F_L = ½ρv²C_L A.
    /// </summary>
    /// <param name="liftCoefficient">Lift coefficient C_L (dimensionless).</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="speed">Flow speed relative to the body in m/s.</param>
    /// <param name="referenceArea">Reference (wing) area A in m².</param>
    public static double LiftForce(
        this double liftCoefficient, double density, double speed, double referenceArea)
    {
        return 0.5 * density * speed * speed * liftCoefficient * referenceArea;
    }

    /// <summary>
    /// Computes the terminal velocity of a falling object:
    /// v_t = √(2mg / (ρC_D A)).
    /// </summary>
    /// <param name="mass">Mass of the object in kg.</param>
    /// <param name="dragCoefficient">Drag coefficient C_D.</param>
    /// <param name="density">Fluid density ρ in kg/m³.</param>
    /// <param name="referenceArea">Reference (frontal) area A in m².</param>
    public static double TerminalVelocity(
        this double mass, double dragCoefficient, double density, double referenceArea)
    {
        if (dragCoefficient <= 0 || density <= 0 || referenceArea <= 0)
            throw new ArgumentException("Drag coefficient, density, and area must be positive.");
        return Math.Sqrt(2.0 * mass * PhysicsConstants.GravitationalAcceleration
            / (density * dragCoefficient * referenceArea));
    }

    #endregion
}
