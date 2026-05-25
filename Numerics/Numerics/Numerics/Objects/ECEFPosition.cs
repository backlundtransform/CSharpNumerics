using System;

namespace CSharpNumerics.Numerics.Objects;

/// <summary>
/// Represents a position in the Earth-Centered Earth-Fixed (ECEF) coordinate system.
/// X points to 0° latitude, 0° longitude (intersection of equator and prime meridian).
/// Y points to 0° latitude, 90°E longitude.
/// Z points to the North Pole.
/// Units: meters.
/// </summary>
public struct ECEFPosition
{
    public double X;
    public double Y;
    public double Z;

    public ECEFPosition(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    /// <summary>Distance from Earth's center in meters.</summary>
    public double Radius => Math.Sqrt(X * X + Y * Y + Z * Z);

    public static ECEFPosition operator +(ECEFPosition a, ECEFPosition b)
        => new ECEFPosition(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

    public static ECEFPosition operator -(ECEFPosition a, ECEFPosition b)
        => new ECEFPosition(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

    public static ECEFPosition operator *(double s, ECEFPosition a)
        => new ECEFPosition(s * a.X, s * a.Y, s * a.Z);

    /// <summary>Converts to a Vector (x, y, z).</summary>
    public Vector ToVector() => new Vector(X, Y, Z);

    /// <summary>Creates from a Vector.</summary>
    public static ECEFPosition FromVector(Vector v) => new ECEFPosition(v.x, v.y, v.z);
}
