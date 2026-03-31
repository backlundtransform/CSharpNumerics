using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics.Electromagnetism;

/// <summary>
/// Extension methods for electrostatics: Coulomb's law, point-charge fields,
/// electric potentials, superposition, and electric dipoles.
/// </summary>
public static class ElectrostaticsExtensions
{
    /// <summary>
    /// Coulomb constant: k = 1 / (4πε₀) ≈ 8.988 × 10⁹ N·m²/C².
    /// </summary>
    public static readonly double CoulombConstant =
        1.0 / (4.0 * Math.PI * PhysicsConstants.VacuumPermittivity);

    // ═══════════════════════════════════════════════════════════════
    //  Point-charge fields
    // ═══════════════════════════════════════════════════════════════

    #region Point Charges

    /// <summary>
    /// Computes the electric field from a point charge at a field point:
    /// E = kq / |r|² · r̂, where r = fieldPoint − chargePosition.
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="chargePosition">Position of the charge.</param>
    /// <param name="fieldPoint">Position where the field is evaluated.</param>
    public static Vector ElectricField(this double charge, Vector chargePosition, Vector fieldPoint)
    {
        var r = fieldPoint - chargePosition;
        var mag = r.GetMagnitude();
        if (mag == 0) throw new ArgumentException("Field point cannot coincide with charge position.");
        return CoulombConstant * charge / (mag * mag * mag) * r;
    }

    /// <summary>
    /// Computes the electric potential from a point charge: V = kq / |r|.
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="chargePosition">Position of the charge.</param>
    /// <param name="fieldPoint">Position where the potential is evaluated.</param>
    public static double ElectricPotential(this double charge, Vector chargePosition, Vector fieldPoint)
    {
        var r = (fieldPoint - chargePosition).GetMagnitude();
        if (r == 0) throw new ArgumentException("Field point cannot coincide with charge position.");
        return CoulombConstant * charge / r;
    }

    /// <summary>
    /// Computes the Coulomb force on charge1 due to charge2:
    /// F₁₂ = kq₁q₂ / |r|² · r̂, where r = position1 − position2.
    /// Positive → repulsion, negative → attraction.
    /// </summary>
    /// <param name="charge1">First charge in coulombs.</param>
    /// <param name="charge2">Second charge in coulombs.</param>
    /// <param name="position1">Position of the first charge.</param>
    /// <param name="position2">Position of the second charge.</param>
    public static Vector CoulombForce(this double charge1, double charge2, Vector position1, Vector position2)
    {
        var r = position1 - position2;
        var mag = r.GetMagnitude();
        if (mag == 0) throw new ArgumentException("Charges cannot occupy the same position.");
        return CoulombConstant * charge1 * charge2 / (mag * mag * mag) * r;
    }

    /// <summary>
    /// Computes the superposed electric field from multiple point charges at a field point.
    /// E(r) = Σ kqᵢ / |r − rᵢ|² · (r − rᵢ) / |r − rᵢ|.
    /// </summary>
    /// <param name="charges">Array of (charge, position) tuples.</param>
    /// <param name="fieldPoint">Position where the field is evaluated.</param>
    public static Vector ElectricFieldSuperposition(
        this (double charge, Vector position)[] charges, Vector fieldPoint)
    {
        var sum = new Vector(0, 0, 0);
        for (int i = 0; i < charges.Length; i++)
            sum = sum + charges[i].charge.ElectricField(charges[i].position, fieldPoint);
        return sum;
    }

    /// <summary>
    /// Computes the superposed electric potential from multiple point charges.
    /// V(r) = Σ kqᵢ / |r − rᵢ|.
    /// </summary>
    /// <param name="charges">Array of (charge, position) tuples.</param>
    /// <param name="fieldPoint">Position where the potential is evaluated.</param>
    public static double ElectricPotentialSuperposition(
        this (double charge, Vector position)[] charges, Vector fieldPoint)
    {
        double sum = 0;
        for (int i = 0; i < charges.Length; i++)
            sum += charges[i].charge.ElectricPotential(charges[i].position, fieldPoint);
        return sum;
    }

    /// <summary>
    /// Creates a <see cref="ScalarField"/> representing the electric potential from a point charge:
    /// V(r) = kq / |r − r_q|.
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="chargePosition">Position of the charge.</param>
    public static ScalarField ElectricScalarPotential(this double charge, Vector chargePosition)
    {
        var k = CoulombConstant;
        return new ScalarField(r => k * charge / (r - chargePosition).GetMagnitude());
    }

    /// <summary>
    /// Creates a <see cref="ScalarField"/> representing the superposed electric potential
    /// from multiple point charges: V(r) = Σ kqᵢ / |r − rᵢ|.
    /// </summary>
    /// <param name="charges">Array of (charge, position) tuples.</param>
    public static ScalarField ElectricScalarPotential(this (double charge, Vector position)[] charges)
    {
        var k = CoulombConstant;
        return new ScalarField(r =>
        {
            double sum = 0;
            for (int i = 0; i < charges.Length; i++)
                sum += k * charges[i].charge / (r - charges[i].position).GetMagnitude();
            return sum;
        });
    }

    /// <summary>
    /// Creates a <see cref="VectorField"/> representing the electric field from a point charge.
    /// Computed as E = −∇V via <see cref="ElectricScalarPotential"/>.
    /// The returned field can be used with ∇· and ∇× operators.
    /// </summary>
    /// <param name="charge">Charge in coulombs.</param>
    /// <param name="chargePosition">Position of the charge.</param>
    public static VectorField ElectricVectorField(this double charge, Vector chargePosition)
    {
        return (-charge.ElectricScalarPotential(chargePosition)).GradientField();
    }

    /// <summary>
    /// Creates a <see cref="VectorField"/> representing the superposed electric field
    /// from multiple point charges. Computed as E = −∇V via <see cref="ElectricScalarPotential"/>.
    /// </summary>
    /// <param name="charges">Array of (charge, position) tuples.</param>
    public static VectorField ElectricVectorField(this (double charge, Vector position)[] charges)
    {
        return (-charges.ElectricScalarPotential()).GradientField();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════
    //  Dipole fields
    // ═══════════════════════════════════════════════════════════════

    #region Dipoles

    /// <summary>
    /// Computes the electric potential from a dipole at a field point:
    /// V = (1/4πε₀) · (p · r̂) / r², where p is the dipole moment and r = fieldPoint − dipolePosition.
    /// </summary>
    /// <param name="dipoleMoment">Electric dipole moment vector (p = qd) in C·m.</param>
    /// <param name="dipolePosition">Position of the dipole.</param>
    /// <param name="fieldPoint">Position where the potential is evaluated.</param>
    public static double DipolePotential(this Vector dipoleMoment, Vector dipolePosition, Vector fieldPoint)
    {
        var r = fieldPoint - dipolePosition;
        var mag = r.GetMagnitude();
        if (mag == 0) throw new ArgumentException("Field point cannot coincide with dipole position.");
        return CoulombConstant * dipoleMoment.Dot(r) / (mag * mag * mag);
    }

    /// <summary>
    /// Computes the electric field from a dipole at a field point:
    /// E = (1/4πε₀) · [3(p·r̂)r̂ − p] / r³.
    /// </summary>
    /// <param name="dipoleMoment">Electric dipole moment vector in C·m.</param>
    /// <param name="dipolePosition">Position of the dipole.</param>
    /// <param name="fieldPoint">Position where the field is evaluated.</param>
    public static Vector DipoleElectricField(this Vector dipoleMoment, Vector dipolePosition, Vector fieldPoint)
    {
        var r = fieldPoint - dipolePosition;
        var mag = r.GetMagnitude();
        if (mag == 0) throw new ArgumentException("Field point cannot coincide with dipole position.");
        var rHat = r.GetUnitVector();
        var pDotRhat = dipoleMoment.Dot(rHat);
        return CoulombConstant / (mag * mag * mag) * (3.0 * pDotRhat * rHat - dipoleMoment);
    }

    #endregion
}
