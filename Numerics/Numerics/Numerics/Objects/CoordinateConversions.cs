using System;
using CSharpNumerics.Physics.Mechanics;

namespace CSharpNumerics.Numerics.Objects;

/// <summary>
/// Provides conversions between geodetic (WGS84), ECEF, and ECI coordinate systems.
/// </summary>
public static class CoordinateConversions
{
    /// <summary>
    /// Converts geodetic coordinates to ECEF.
    /// </summary>
    public static ECEFPosition GeodeticToECEF(GeoCoordinate geo)
    {
        double sinLat = Math.Sin(geo.Latitude);
        double cosLat = Math.Cos(geo.Latitude);
        double sinLon = Math.Sin(geo.Longitude);
        double cosLon = Math.Cos(geo.Longitude);

        double N = EarthModel.PrimeVerticalRadius(geo.Latitude);

        double x = (N + geo.Altitude) * cosLat * cosLon;
        double y = (N + geo.Altitude) * cosLat * sinLon;
        double z = (N * (1.0 - EarthModel.EccentricitySquared) + geo.Altitude) * sinLat;

        return new ECEFPosition(x, y, z);
    }

    /// <summary>
    /// Converts ECEF to geodetic coordinates using the iterative Bowring method.
    /// Converges to sub-millimeter accuracy in 2-3 iterations.
    /// </summary>
    public static GeoCoordinate ECEFToGeodetic(ECEFPosition ecef)
    {
        double x = ecef.X;
        double y = ecef.Y;
        double z = ecef.Z;

        double lon = Math.Atan2(y, x);

        double p = Math.Sqrt(x * x + y * y);
        double e2 = EarthModel.EccentricitySquared;
        double a = EarthModel.SemiMajorAxis;
        double b = EarthModel.SemiMinorAxis;

        // Initial estimate using Bowring's method
        double theta = Math.Atan2(z * a, p * b);
        double sinTheta = Math.Sin(theta);
        double cosTheta = Math.Cos(theta);

        double lat = Math.Atan2(
            z + EarthModel.SecondEccentricitySquared * b * sinTheta * sinTheta * sinTheta,
            p - e2 * a * cosTheta * cosTheta * cosTheta);

        // Iterate for convergence
        for (int i = 0; i < 5; i++)
        {
            double sinLat = Math.Sin(lat);
            double cosLat = Math.Cos(lat);
            double N = a / Math.Sqrt(1.0 - e2 * sinLat * sinLat);

            lat = Math.Atan2(z + e2 * N * sinLat, p);
        }

        double sinLatFinal = Math.Sin(lat);
        double Nfinal = a / Math.Sqrt(1.0 - e2 * sinLatFinal * sinLatFinal);
        double cosLatFinal = Math.Cos(lat);

        double alt;
        if (Math.Abs(cosLatFinal) > 1e-10)
            alt = p / cosLatFinal - Nfinal;
        else
            alt = Math.Abs(z) / Math.Abs(sinLatFinal) - Nfinal * (1.0 - e2);

        return new GeoCoordinate(lat, lon, alt);
    }

    /// <summary>
    /// Converts ECEF position to ECI position given Greenwich Mean Sidereal Time angle in radians.
    /// </summary>
    public static ECIPosition ECEFToECI(ECEFPosition ecef, double gmstRadians)
    {
        double cosG = Math.Cos(gmstRadians);
        double sinG = Math.Sin(gmstRadians);

        double xEci = ecef.X * cosG - ecef.Y * sinG;
        double yEci = ecef.X * sinG + ecef.Y * cosG;
        double zEci = ecef.Z;

        return new ECIPosition(xEci, yEci, zEci);
    }

    /// <summary>
    /// Converts ECI position to ECEF position given Greenwich Mean Sidereal Time angle in radians.
    /// </summary>
    public static ECEFPosition ECIToECEF(ECIPosition eci, double gmstRadians)
    {
        double cosG = Math.Cos(gmstRadians);
        double sinG = Math.Sin(gmstRadians);

        // Inverse rotation (rotate by -gmst)
        double xEcef = eci.X * cosG + eci.Y * sinG;
        double yEcef = -eci.X * sinG + eci.Y * cosG;
        double zEcef = eci.Z;

        return new ECEFPosition(xEcef, yEcef, zEcef);
    }

    /// <summary>
    /// Converts ECEF velocity to ECI velocity given GMST angle and ECEF position.
    /// Accounts for Earth's rotation: v_eci = R·v_ecef + ω×r_eci
    /// </summary>
    public static Vector ECEFVelocityToECI(Vector velocityECEF, ECEFPosition positionECEF, double gmstRadians)
    {
        double cosG = Math.Cos(gmstRadians);
        double sinG = Math.Sin(gmstRadians);

        // Rotate velocity
        double vxEci = velocityECEF.x * cosG - velocityECEF.y * sinG;
        double vyEci = velocityECEF.x * sinG + velocityECEF.y * cosG;
        double vzEci = velocityECEF.z;

        // Add ω × r_eci (ω = [0, 0, ω_earth])
        ECIPosition posEci = ECEFToECI(positionECEF, gmstRadians);
        double omega = EarthModel.RotationRate;

        vxEci += -omega * posEci.Y;
        vyEci += omega * posEci.X;

        return new Vector(vxEci, vyEci, vzEci);
    }

    /// <summary>
    /// Converts ECI velocity to ECEF velocity given GMST angle and ECI position.
    /// v_ecef = R^T·(v_eci - ω×r_eci)
    /// </summary>
    public static Vector ECIVelocityToECEF(Vector velocityECI, ECIPosition positionECI, double gmstRadians)
    {
        double omega = EarthModel.RotationRate;

        // Subtract ω × r_eci
        double vxInertial = velocityECI.x - (-omega * positionECI.Y);
        double vyInertial = velocityECI.y - (omega * positionECI.X);
        double vzInertial = velocityECI.z;

        // Inverse rotation
        double cosG = Math.Cos(gmstRadians);
        double sinG = Math.Sin(gmstRadians);

        double vxEcef = vxInertial * cosG + vyInertial * sinG;
        double vyEcef = -vxInertial * sinG + vyInertial * cosG;
        double vzEcef = vzInertial;

        return new Vector(vxEcef, vyEcef, vzEcef);
    }

    /// <summary>
    /// Convenience: Geodetic directly to ECI.
    /// </summary>
    public static ECIPosition GeodeticToECI(GeoCoordinate geo, double gmstRadians)
    {
        return ECEFToECI(GeodeticToECEF(geo), gmstRadians);
    }

    /// <summary>
    /// Convenience: ECI directly to Geodetic.
    /// </summary>
    public static GeoCoordinate ECIToGeodetic(ECIPosition eci, double gmstRadians)
    {
        return ECEFToGeodetic(ECIToECEF(eci, gmstRadians));
    }
}
