using System;

namespace CSharpNumerics.Numerics.Objects;

/// <summary>
/// Represents a geodetic coordinate (latitude, longitude, altitude) on the WGS84 ellipsoid.
/// </summary>
public struct GeoCoordinate
{
    /// <summary>Geodetic latitude in radians. Positive north.</summary>
    public double Latitude;

    /// <summary>Longitude in radians. Positive east.</summary>
    public double Longitude;

    /// <summary>Altitude above the WGS84 ellipsoid in meters.</summary>
    public double Altitude;

    public GeoCoordinate(double latitude, double longitude, double altitude)
    {
        Latitude = latitude;
        Longitude = longitude;
        Altitude = altitude;
    }

    /// <summary>Creates a GeoCoordinate from degrees.</summary>
    public static GeoCoordinate FromDegrees(double latitudeDeg, double longitudeDeg, double altitudeMeters)
    {
        return new GeoCoordinate(
            latitudeDeg * Math.PI / 180.0,
            longitudeDeg * Math.PI / 180.0,
            altitudeMeters);
    }

    /// <summary>Latitude in degrees.</summary>
    public double LatitudeDegrees => Latitude * 180.0 / Math.PI;

    /// <summary>Longitude in degrees.</summary>
    public double LongitudeDegrees => Longitude * 180.0 / Math.PI;
}
