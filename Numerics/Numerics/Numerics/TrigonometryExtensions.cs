using System;
namespace CSharpNumerics.Numerics;

public static class TrigonometryExtensions
{
    /// <summary>
    /// Converts an angle from degrees to radians.
    /// radians = degrees * π / 180.
    /// </summary>
    /// <param name="angle">Angle in degrees.</param>
    /// <returns>Angle in radians.</returns>
    public static double DegreeToRadians(this double angle)
    {
        return (angle * Math.PI) / 180;
    }
}
