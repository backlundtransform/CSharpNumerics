using System;
using System.Globalization;

namespace CSharpNumerics.Engines.GIS.Coordinates
{
    /// <summary>
    /// An immutable WGS-84 geographic coordinate (latitude, longitude, altitude).
    /// Latitude is in degrees north (−90 to +90), longitude in degrees east (−180 to +180),
    /// altitude in metres above the WGS-84 ellipsoid.
    /// </summary>
    public readonly struct GeoCoordinate : IEquatable<GeoCoordinate>
    {
        /// <summary>Latitude in decimal degrees (positive = north).</summary>
        public double Latitude { get; }

        /// <summary>Longitude in decimal degrees (positive = east).</summary>
        public double Longitude { get; }

        /// <summary>Altitude in metres above the WGS-84 ellipsoid.</summary>
        public double Altitude { get; }

        /// <summary>
        /// Creates a new geographic coordinate.
        /// </summary>
        /// <param name="latitude">Latitude in degrees (−90 to +90).</param>
        /// <param name="longitude">Longitude in degrees (−180 to +180).</param>
        /// <param name="altitude">Altitude in metres (default 0).</param>
        public GeoCoordinate(double latitude, double longitude, double altitude = 0)
        {
            if (latitude < -90 || latitude > 90)
                throw new ArgumentOutOfRangeException(nameof(latitude), "Latitude must be between −90 and +90.");
            if (longitude < -180 || longitude > 180)
                throw new ArgumentOutOfRangeException(nameof(longitude), "Longitude must be between −180 and +180.");

            Latitude = latitude;
            Longitude = longitude;
            Altitude = altitude;
        }

        /// <summary>
        /// Haversine distance in metres to another coordinate (ignoring altitude).
        /// </summary>
        public double DistanceTo(GeoCoordinate other)
        {
            const double R = 6371000.0; // mean Earth radius in metres
            double lat1 = Latitude * Math.PI / 180.0;
            double lat2 = other.Latitude * Math.PI / 180.0;
            double dLat = (other.Latitude - Latitude) * Math.PI / 180.0;
            double dLon = (other.Longitude - Longitude) * Math.PI / 180.0;

            double a = Math.Sin(dLat / 2) * Math.Sin(dLat / 2) +
                       Math.Cos(lat1) * Math.Cos(lat2) *
                       Math.Sin(dLon / 2) * Math.Sin(dLon / 2);
            double c = 2 * Math.Atan2(Math.Sqrt(a), Math.Sqrt(1 - a));
            return R * c;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Equality
        // ═══════════════════════════════════════════════════════════════

        public bool Equals(GeoCoordinate other) =>
            Latitude == other.Latitude && Longitude == other.Longitude && Altitude == other.Altitude;

        public override bool Equals(object obj) => obj is GeoCoordinate g && Equals(g);

        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 31 + Latitude.GetHashCode();
                hash = hash * 31 + Longitude.GetHashCode();
                hash = hash * 31 + Altitude.GetHashCode();
                return hash;
            }
        }

        public static bool operator ==(GeoCoordinate a, GeoCoordinate b) => a.Equals(b);
        public static bool operator !=(GeoCoordinate a, GeoCoordinate b) => !a.Equals(b);

        public override string ToString() =>
            string.Format(CultureInfo.InvariantCulture,
                "({0:F6}°N, {1:F6}°E, {2:F1}m)", Latitude, Longitude, Altitude);
    }
}
