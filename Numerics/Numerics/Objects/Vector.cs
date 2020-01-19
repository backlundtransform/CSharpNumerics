using System;

namespace Numerics.Objects
{
    public struct Vector
    {
        public double X;

        public double Y;

        public double Z;


        public Vector(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public Vector((double, double, double) p1, (double, double, double) p2)
        {
            X = p2.Item1 - p1.Item1;
            Y = p2.Item2 - p1.Item2;
            Z = p2.Item3 - p1.Item3;
        }
        public double GetMagnitude() => Math.Sqrt(Math.Pow(X, 2) + Math.Pow(Y, 2) + Math.Pow(Z, 2));

        public double GetInclination() => Math.Acos(Z/GetMagnitude());

        public double GetAzimuth() => Math.Atan2(Y, X);


        public Vector GetUnitVector() => (1 / GetMagnitude()) * new Vector(X, Y, Z);

        public static Vector FromSphericalCoordinates(double radius, double inclination, double azimuth) =>
            new Vector(radius * Math.Sin(inclination) * Math.Cos(azimuth),
                radius * Math.Sin(inclination) * Math.Sin(azimuth), radius * Math.Cos(inclination));

        public static Vector operator *(double a, Vector b) => new Vector(a * b.X, a * b.Y, a * b.Z);

        public static Vector operator -(Vector a, Vector b) => new Vector(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

        public static Vector operator +(Vector a, Vector b) => new Vector(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

        public static double operator *(Vector a, Vector b) => a.Dot(b);

        public override string ToString() => $"{X}*î̤+{Y}*ĵ+{Z}*k̂";

    }
}
