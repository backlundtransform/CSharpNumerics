using System;

namespace Numerics.Objects
{
    public struct Vector
    {
        public double x;

        public double y;

        public double z;


        public Vector(double a, double b, double c)
        {
            x = a;
            y = b;
            z =c;
        }

        public Vector((double, double, double) p1, (double, double, double) p2)
        {
            x = p2.Item1 - p1.Item1;
            y = p2.Item2 - p1.Item2;
            z = p2.Item3 - p1.Item3;
        }
        public double GetMagnitude() => Math.Sqrt(Math.Pow(x, 2) + Math.Pow(y, 2) + Math.Pow(z, 2));

        public double GetInclination() => Math.Acos(z/GetMagnitude());

        public double GetAzimuth() => Math.Atan2(y, x);


        public Vector GetUnitVector() => (1 / GetMagnitude()) * new Vector(x, y, z);

        public static Vector FromSphericalCoordinates(double radius, double inclination, double azimuth) =>
            new Vector(radius * Math.Sin(inclination) * Math.Cos(azimuth),
                radius * Math.Sin(inclination) * Math.Sin(azimuth), radius * Math.Cos(inclination));

        public static Vector operator *(double a, Vector b) => new Vector(a * b.x, a * b.y, a * b.z);

        public static Vector operator -(Vector a, Vector b) => new Vector(a.x - b.x, a.y - b.y, a.z - b.z);

        public static Vector operator +(Vector a, Vector b) => new Vector(a.x + b.x, a.y + b.y, a.z + b.z);

        public static double operator *(Vector a, Vector b) => a.Dot(b);

        public override string ToString() => $"{x}*î̤+{y}*ĵ+{z}*k̂";

        public double Dot(Vector b) => x * b.x + y * b.y + z * b.z;
        public Vector Cross(Vector b) => new Vector(y * b.z - b.y * z, -(x * b.z - b.x * z), x * b.y - b.x * y);
        public double GetAngle(Vector b) => Math.Acos(Dot(b) / (GetMagnitude() * b.GetMagnitude()));
        public  Vector Projection(Vector b) => b.Dot(this) / Math.Pow(GetMagnitude(), 2) * this;
        public  Vector Reflection( Vector b) => 2.0 * Projection(b) - b;
        public Vector ToSphericalCoordinates() => new Vector(GetMagnitude(), GetInclination(), GetAzimuth());

    }
}
