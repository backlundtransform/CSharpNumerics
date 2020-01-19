using System;

namespace Numerics.Objects
{
    public static class VectorExtensions
    {
        public static double Dot(this Vector a, Vector b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;
        public static Vector Cross(this Vector a, Vector b) => new Vector(a.Y * b.Z - b.Y * a.Z, -(a.X * b.Z - b.X * a.Z), a.X * b.Y - b.X * a.Y);
        public static double GetAngle(this Vector a, Vector b) => Math.Acos(a.Dot(b) / (a.GetMagnitude() * b.GetMagnitude()));
        public static Vector Projection(this Vector a, Vector b)=> b.Dot(a) / Math.Pow(a.GetMagnitude(), 2) * a;
        public static Vector Reflection(this Vector a, Vector b) => 2.0*a.Projection(b)-b;
        public static Vector ToSphericalCoordinates(this Vector a) =>new Vector(a.GetMagnitude(),a.GetInclination(),a.GetAzimuth());
    }
}
