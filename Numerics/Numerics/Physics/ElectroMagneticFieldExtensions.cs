using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Enums;
using CSharpNumerics.Numerics;
using System;

namespace CSharpNumerics.Physics
{
    /// <summary>
    /// Provides extension methods for electromagnetic field calculations:
    /// Coulomb's law, Lorentz force, Maxwell's equations (differential form),
    /// Poynting vector, energy density, and electromagnetic potentials.
    /// <para>
    /// Bridges <see cref="VectorField"/> (∇·, ∇×) and
    /// <see cref="VectorFieldExtensions"/> (∇, ∇²) to classical electrodynamics.
    /// </para>
    /// </summary>
    public static class ElectroMagneticFieldExtensions
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
        //  Lorentz force
        // ═══════════════════════════════════════════════════════════════

        #region Lorentz Force

        /// <summary>
        /// Computes the Lorentz force on a charged particle: F = q(E + v × B).
        /// </summary>
        /// <param name="charge">Charge in coulombs.</param>
        /// <param name="velocity">Velocity of the particle in m/s.</param>
        /// <param name="electricField">Electric field at the particle position.</param>
        /// <param name="magneticField">Magnetic field at the particle position.</param>
        public static Vector LorentzForce(this double charge, Vector velocity, Vector electricField, Vector magneticField)
        {
            return charge * (electricField + velocity.Cross(magneticField));
        }

        /// <summary>
        /// Computes only the electric part of the Lorentz force: F = qE.
        /// </summary>
        /// <param name="charge">Charge in coulombs.</param>
        /// <param name="electricField">Electric field at the particle position.</param>
        public static Vector ElectricForce(this double charge, Vector electricField)
        {
            return charge * electricField;
        }

        /// <summary>
        /// Computes only the magnetic part of the Lorentz force: F = q(v × B).
        /// </summary>
        /// <param name="charge">Charge in coulombs.</param>
        /// <param name="velocity">Velocity of the particle in m/s.</param>
        /// <param name="magneticField">Magnetic field at the particle position.</param>
        public static Vector MagneticForce(this double charge, Vector velocity, Vector magneticField)
        {
            return charge * velocity.Cross(magneticField);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Maxwell's equations — differential form
        //
        //  (1) ∇·E = ρ/ε₀          Gauss (electric)
        //  (2) ∇·B = 0              Gauss (magnetic)
        //  (3) ∇×E = −∂B/∂t        Faraday
        //  (4) ∇×B = μ₀J + μ₀ε₀∂E/∂t   Ampère–Maxwell
        // ═══════════════════════════════════════════════════════════════

        #region Maxwell's Equations

        /// <summary>
        /// Gauss's law for electricity: computes the charge density at a point
        /// from the electric field divergence.
        /// ρ = ε₀ ∇·E.
        /// </summary>
        /// <param name="electricField">The electric vector field.</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static double ChargeDensity(this VectorField electricField, (double, double, double) point)
        {
            return PhysicsConstants.VacuumPermittivity * electricField.Divergence(point);
        }

        /// <summary>
        /// Gauss's law for magnetism: computes ∇·B which should be zero
        /// for any physical magnetic field. Returns the numerical divergence
        /// for verification.
        /// </summary>
        /// <param name="magneticField">The magnetic vector field.</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static double GaussMagnetic(this VectorField magneticField, (double, double, double) point)
        {
            return magneticField.Divergence(point);
        }

        /// <summary>
        /// Faraday's law of induction: computes −∇×E, which equals ∂B/∂t.
        /// For electrostatic fields, the result should be approximately zero.
        /// </summary>
        /// <param name="electricField">The electric vector field.</param>
        /// <param name="point">The point at which to evaluate.</param>
        /// <returns>The time rate of change of the magnetic field: ∂B/∂t = −∇×E.</returns>
        public static Vector FaradayLaw(this VectorField electricField, (double, double, double) point)
        {
            var curl = electricField.Curl(point);
            return new Vector(-curl.x, -curl.y, -curl.z);
        }

        /// <summary>
        /// Ampère–Maxwell law: computes (1/μ₀)∇×B, which equals J + ε₀ ∂E/∂t.
        /// For magnetostatic fields with no changing E, the result gives
        /// the current density J.
        /// </summary>
        /// <param name="magneticField">The magnetic vector field.</param>
        /// <param name="point">The point at which to evaluate.</param>
        /// <returns>J + ε₀ ∂E/∂t = (1/μ₀)∇×B.</returns>
        public static Vector AmpereLaw(this VectorField magneticField, (double, double, double) point)
        {
            var curl = magneticField.Curl(point);
            var invMu = 1.0 / PhysicsConstants.VacuumPermeability;
            return invMu * curl;
        }

        /// <summary>
        /// Verifies the wave equation for a scalar component of an electromagnetic field:
        /// returns ∇²f − μ₀ε₀ ∂²f/∂t², which should be zero in free space.
        /// </summary>
        /// <param name="field">A scalar field component f(r, t) parameterised as f(t)(r).</param>
        /// <param name="t">Time at which to evaluate.</param>
        /// <param name="point">Spatial point at which to evaluate.</param>
        /// <param name="dt">Time step for numerical second derivative (default 1 × 10⁻⁶).</param>
        public static double WaveEquationResidual(
            this Func<double, Func<Vector, double>> field,
            double t, (double, double, double) point, double dt = 1e-6)
        {
            var laplacian = field(t).Laplacian(point);

            var fPlus = field(t + dt)(new Vector(point));
            var fMid = field(t)(new Vector(point));
            var fMinus = field(t - dt)(new Vector(point));
            var d2fdt2 = (fPlus - 2 * fMid + fMinus) / (dt * dt);

            var mu0eps0 = PhysicsConstants.VacuumPermeability * PhysicsConstants.VacuumPermittivity;
            return laplacian - mu0eps0 * d2fdt2;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Energy & momentum
        // ═══════════════════════════════════════════════════════════════

        #region Energy & Momentum

        /// <summary>
        /// Computes the electromagnetic energy density:
        /// u = ½(ε₀|E|² + |B|²/μ₀).
        /// </summary>
        /// <param name="electricField">Electric field vector at the point.</param>
        /// <param name="magneticField">Magnetic field vector at the point.</param>
        public static double EnergyDensity(this Vector electricField, Vector magneticField)
        {
            var eMag2 = electricField.Dot(electricField);
            var bMag2 = magneticField.Dot(magneticField);
            return 0.5 * (PhysicsConstants.VacuumPermittivity * eMag2
                        + bMag2 / PhysicsConstants.VacuumPermeability);
        }

        /// <summary>
        /// Computes the Poynting vector: S = (1/μ₀)(E × B).
        /// Represents the directional energy flux (power per unit area).
        /// </summary>
        /// <param name="electricField">Electric field vector at the point.</param>
        /// <param name="magneticField">Magnetic field vector at the point.</param>
        public static Vector PoyntingVector(this Vector electricField, Vector magneticField)
        {
            return (1.0 / PhysicsConstants.VacuumPermeability) * electricField.Cross(magneticField);
        }

        /// <summary>
        /// Computes the radiation pressure from the Poynting vector magnitude:
        /// P = |S| / c (for full absorption) or P = 2|S| / c (for full reflection).
        /// </summary>
        /// <param name="electricField">Electric field vector.</param>
        /// <param name="magneticField">Magnetic field vector.</param>
        /// <param name="reflected">If true, assumes total reflection (pressure doubled).</param>
        public static double RadiationPressure(this Vector electricField, Vector magneticField, bool reflected = false)
        {
            var sMag = electricField.PoyntingVector(magneticField).GetMagnitude();
            return (reflected ? 2.0 : 1.0) * sMag / PhysicsConstants.SpeedOfLight;
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Potentials
        // ═══════════════════════════════════════════════════════════════

        #region Potentials

        /// <summary>
        /// Creates an electric vector field from a scalar potential: E = −∇V.
        /// The returned <see cref="VectorField"/> can be used with ∇· and ∇× operators.
        /// </summary>
        /// <param name="potential">The electric potential V(r).</param>
        public static VectorField ElectricFieldFromPotential(this Func<Vector, double> potential)
        {
            return (-new ScalarField(potential)).GradientField();
        }

        /// <summary>
        /// Creates an electric vector field from a <see cref="ScalarField"/> potential: E = −∇V.
        /// </summary>
        /// <param name="potential">The electric potential V(r) as a ScalarField.</param>
        public static VectorField ElectricFieldFromPotential(this ScalarField potential)
        {
            return (-potential).GradientField();
        }

        /// <summary>
        /// Computes the magnetic field from a vector potential at a point: B = ∇ × A.
        /// </summary>
        /// <param name="vectorPotential">The magnetic vector potential A(r).</param>
        /// <param name="point">The point at which to evaluate.</param>
        public static Vector MagneticFieldFromVectorPotential(
            this VectorField vectorPotential, (double, double, double) point)
        {
            return vectorPotential.Curl(point);
        }

        #endregion

        // ═══════════════════════════════════════════════════════════════
        //  Magnetic field sources (Biot–Savart)
        // ═══════════════════════════════════════════════════════════════

        #region Biot–Savart

        /// <summary>
        /// Computes the magnetic field magnitude at perpendicular distance d
        /// from an infinite straight wire: B = μ₀I / (2πd).
        /// </summary>
        /// <param name="current">Current in amperes.</param>
        /// <param name="distance">Perpendicular distance from the wire in meters.</param>
        public static double MagneticFieldFromWire(this double current, double distance)
        {
            if (distance <= 0) throw new ArgumentException("Distance must be greater than zero.");
            return PhysicsConstants.VacuumPermeability * Math.Abs(current) / (2.0 * Math.PI * distance);
        }

        /// <summary>
        /// Computes the magnetic field vector from an infinite straight wire
        /// at a given field point: B = (μ₀I / 2πd) · (ŵ × ρ̂),
        /// where ŵ is the wire direction and ρ̂ is the perpendicular direction
        /// from the wire to the field point.
        /// </summary>
        /// <param name="current">Current in amperes (sign determines direction along wire).</param>
        /// <param name="wireDirection">Direction of current flow (need not be unit length).</param>
        /// <param name="wirePoint">Any point on the wire.</param>
        /// <param name="fieldPoint">Position where the field is evaluated.</param>
        public static Vector MagneticFieldFromWire(this double current, Vector wireDirection, Vector wirePoint, Vector fieldPoint)
        {
            var wireUnit = wireDirection.GetUnitVector();
            var r = fieldPoint - wirePoint;
            var rParallel = r.Dot(wireUnit) * wireUnit;
            var rPerp = r - rParallel;
            var d = rPerp.GetMagnitude();
            if (d == 0) throw new ArgumentException("Field point lies on the wire.");

            var bDirection = wireUnit.Cross(rPerp.GetUnitVector());
            var bMagnitude = PhysicsConstants.VacuumPermeability * current / (2.0 * Math.PI * d);

            return bMagnitude * bDirection;
        }

        /// <summary>
        /// Computes the magnetic field at the center of a circular current loop:
        /// B = μ₀I / (2R), directed along the loop normal.
        /// </summary>
        /// <param name="current">Current in amperes.</param>
        /// <param name="radius">Radius of the loop in meters.</param>
        /// <param name="normal">Normal direction of the loop (right-hand rule).</param>
        public static Vector MagneticFieldFromLoop(this double current, double radius, Vector normal)
        {
            if (radius <= 0) throw new ArgumentException("Radius must be greater than zero.");
            var bMagnitude = PhysicsConstants.VacuumPermeability * current / (2.0 * radius);
            return bMagnitude * normal.GetUnitVector();
        }

        /// <summary>
        /// Computes the magnetic field inside an ideal solenoid: B = μ₀nI,
        /// where n = turnsPerUnitLength.
        /// </summary>
        /// <param name="current">Current in amperes.</param>
        /// <param name="turnsPerUnitLength">Number of turns per meter.</param>
        /// <param name="direction">Direction along the solenoid axis.</param>
        public static Vector MagneticFieldFromSolenoid(this double current, double turnsPerUnitLength, Vector direction)
        {
            var bMagnitude = PhysicsConstants.VacuumPermeability * turnsPerUnitLength * current;
            return bMagnitude * direction.GetUnitVector();
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
}
