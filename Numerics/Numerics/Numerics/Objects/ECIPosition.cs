using System;

namespace CSharpNumerics.Numerics.Objects;

/// <summary>
/// Represents a position in the Earth-Centered Inertial (ECI) coordinate system (J2000).
/// X points to the vernal equinox.
/// Z points to the celestial north pole (aligned with Earth's rotation axis at J2000).
/// Units: meters.
/// </summary>
public struct ECIPosition
{
    public double X;
    public double Y;
    public double Z;

    public ECIPosition(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    /// <summary>Distance from Earth's center in meters.</summary>
    public double Radius => Math.Sqrt(X * X + Y * Y + Z * Z);

    public static ECIPosition operator +(ECIPosition a, ECIPosition b)
        => new ECIPosition(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

    public static ECIPosition operator -(ECIPosition a, ECIPosition b)
        => new ECIPosition(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

    public static ECIPosition operator *(double s, ECIPosition a)
        => new ECIPosition(s * a.X, s * a.Y, s * a.Z);

    /// <summary>Converts to a Vector (x, y, z).</summary>
    public Vector ToVector() => new Vector(X, Y, Z);

    /// <summary>Creates from a Vector.</summary>
    public static ECIPosition FromVector(Vector v) => new ECIPosition(v.x, v.y, v.z);
}
