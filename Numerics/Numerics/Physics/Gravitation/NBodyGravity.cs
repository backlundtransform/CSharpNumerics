using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Constants;

namespace CSharpNumerics.Physics.Gravitation;

/// <summary>
/// Newtonian gravitation for systems of point masses by direct summation.
/// Works with arbitrary bodies (not Earth-specific) and supports Plummer
/// softening to tame the singularity during close encounters. Positions are in
/// metres, masses in kg, accelerations in m/s².
/// </summary>
public static class NBodyGravity
{
    private const double G = PhysicsConstants.GravitationalConstant;

    /// <summary>
    /// Gravitational acceleration on a body at <paramref name="position"/> due to
    /// a single point mass <paramref name="otherMass"/> at <paramref name="otherPosition"/>:
    /// a = G·m·(r_other − r)/|r_other − r|³.
    /// </summary>
    public static Vector PairwiseAcceleration(Vector position, Vector otherPosition, double otherMass, double softening = 0.0)
    {
        double dx = otherPosition.x - position.x;
        double dy = otherPosition.y - position.y;
        double dz = otherPosition.z - position.z;

        double r2 = dx * dx + dy * dy + dz * dz + softening * softening;
        if (r2 <= 0) return new Vector(0, 0, 0);

        double r = Math.Sqrt(r2);
        double factor = G * otherMass / (r2 * r);

        return new Vector(factor * dx, factor * dy, factor * dz);
    }

    /// <summary>
    /// Acceleration on body <paramref name="index"/> from every other body in the
    /// system: a_i = Σ_{j≠i} G·m_j·(r_j − r_i)/|r_j − r_i|³.
    /// </summary>
    public static Vector Acceleration(Vector[] positions, double[] masses, int index, double softening = 0.0)
    {
        ValidateInputs(positions, masses);
        if (index < 0 || index >= positions.Length)
            throw new ArgumentOutOfRangeException(nameof(index));

        double ax = 0, ay = 0, az = 0;
        Vector ri = positions[index];

        for (int j = 0; j < positions.Length; j++)
        {
            if (j == index) continue;
            Vector a = PairwiseAcceleration(ri, positions[j], masses[j], softening);
            ax += a.x; ay += a.y; az += a.z;
        }

        return new Vector(ax, ay, az);
    }

    /// <summary>
    /// Accelerations on all bodies in the system (one per body).
    /// </summary>
    public static Vector[] Accelerations(Vector[] positions, double[] masses, double softening = 0.0)
    {
        ValidateInputs(positions, masses);

        var result = new Vector[positions.Length];
        for (int i = 0; i < positions.Length; i++)
            result[i] = Acceleration(positions, masses, i, softening);

        return result;
    }

    /// <summary>
    /// Total gravitational potential energy of the system:
    /// U = −Σ_{i&lt;j} G·m_i·m_j/r_ij (joules).
    /// </summary>
    public static double PotentialEnergy(Vector[] positions, double[] masses)
    {
        ValidateInputs(positions, masses);

        double u = 0.0;
        for (int i = 0; i < positions.Length; i++)
        {
            for (int j = i + 1; j < positions.Length; j++)
            {
                double dx = positions[j].x - positions[i].x;
                double dy = positions[j].y - positions[i].y;
                double dz = positions[j].z - positions[i].z;
                double r = Math.Sqrt(dx * dx + dy * dy + dz * dz);
                if (r > 0)
                    u -= G * masses[i] * masses[j] / r;
            }
        }

        return u;
    }

    /// <summary>
    /// Centre of mass (barycentre) of the system:
    /// R = Σ m_i·r_i / Σ m_i.
    /// </summary>
    public static Vector CenterOfMass(Vector[] positions, double[] masses)
    {
        ValidateInputs(positions, masses);

        double totalMass = 0, x = 0, y = 0, z = 0;
        for (int i = 0; i < positions.Length; i++)
        {
            totalMass += masses[i];
            x += masses[i] * positions[i].x;
            y += masses[i] * positions[i].y;
            z += masses[i] * positions[i].z;
        }

        if (totalMass <= 0) throw new ArgumentException("Total mass must be positive.", nameof(masses));
        return new Vector(x / totalMass, y / totalMass, z / totalMass);
    }

    private static void ValidateInputs(Vector[] positions, double[] masses)
    {
        if (positions == null) throw new ArgumentNullException(nameof(positions));
        if (masses == null) throw new ArgumentNullException(nameof(masses));
        if (positions.Length != masses.Length)
            throw new ArgumentException("positions and masses must have the same length.");
    }
}
