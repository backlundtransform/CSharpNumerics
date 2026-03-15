using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Engines.GIS.Coordinates
{
    /// <summary>
    /// Projects between WGS-84 geographic coordinates and a local Cartesian
    /// coordinate system (metres). Two modes are supported:
    /// <list type="bullet">
    /// <item><description>
    /// <b>Local Tangent Plane (ENU)</b> — east/north/up centred on an origin.
    /// Good for areas up to ~50 km.
    /// </description></item>
    /// <item><description>
    /// <b>UTM</b> — Universal Transverse Mercator. Automatically determines
    /// the zone from the origin longitude.
    /// </description></item>
    /// </list>
    /// </summary>
    public class Projection
    {
        private const double Deg2Rad = Math.PI / 180.0;
        private const double Rad2Deg = 180.0 / Math.PI;

        // WGS-84 ellipsoid parameters
        private const double A = 6378137.0;            // semi-major axis (m)
        private const double F = 1.0 / 298.257223563;  // flattening
        private const double E2 = 2 * F - F * F;       // eccentricity squared

        /// <summary>The projection origin in geographic coordinates.</summary>
        public GeoCoordinate Origin { get; }

        /// <summary>The projection type.</summary>
        public ProjectionType Type { get; }

        // UTM zone fields
        private readonly int _utmZone;
        private readonly bool _utmNorth;

        // Precomputed values for local tangent plane
        private readonly double _cosLat0;
        private readonly double _sinLat0;
        private readonly double _mPerDegLat;
        private readonly double _mPerDegLon;

        /// <summary>
        /// Creates a projection centred on <paramref name="origin"/>.
        /// </summary>
        /// <param name="origin">The geographic origin (becomes (0,0) in local coordinates).</param>
        /// <param name="type">Projection method to use.</param>
        public Projection(GeoCoordinate origin, ProjectionType type = ProjectionType.LocalTangentPlane)
        {
            Origin = origin;
            Type = type;

            double latRad = origin.Latitude * Deg2Rad;
            _cosLat0 = Math.Cos(latRad);
            _sinLat0 = Math.Sin(latRad);

            // Radii of curvature at origin latitude
            double N0 = A / Math.Sqrt(1 - E2 * _sinLat0 * _sinLat0);
            double M0 = A * (1 - E2) / Math.Pow(1 - E2 * _sinLat0 * _sinLat0, 1.5);

            _mPerDegLat = M0 * Deg2Rad;
            _mPerDegLon = N0 * _cosLat0 * Deg2Rad;

            // UTM zone from origin longitude
            _utmZone = (int)Math.Floor((origin.Longitude + 180.0) / 6.0) + 1;
            _utmNorth = origin.Latitude >= 0;
        }

        /// <summary>UTM zone number (1–60). Only meaningful when <see cref="Type"/> is <see cref="ProjectionType.UTM"/>.</summary>
        public int UtmZone => _utmZone;

        /// <summary>Whether the UTM zone is in the northern hemisphere.</summary>
        public bool UtmNorth => _utmNorth;

        // ═══════════════════════════════════════════════════════════════
        //  Forward: geographic → local metres
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Projects a geographic coordinate to local (x, y, z) in metres.
        /// x = east, y = north, z = altitude above origin.
        /// </summary>
        public Vector ToLocal(GeoCoordinate coord)
        {
            if (Type == ProjectionType.UTM)
                return ToLocalUtm(coord);

            return ToLocalTangentPlane(coord);
        }

        /// <summary>
        /// Projects geographic coordinates to a local (x, y, z) tuple.
        /// Convenience overload accepting raw lat/lon/alt.
        /// </summary>
        public Vector ToLocal(double latitude, double longitude, double altitude = 0)
            => ToLocal(new GeoCoordinate(latitude, longitude, altitude));

        // ═══════════════════════════════════════════════════════════════
        //  Inverse: local metres → geographic
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Converts local (x, y, z) metres back to a geographic coordinate.
        /// </summary>
        public GeoCoordinate ToGeo(Vector local)
        {
            if (Type == ProjectionType.UTM)
                return ToGeoUtm(local);

            return ToGeoTangentPlane(local);
        }

        /// <summary>Convenience overload accepting (x, y, z) as doubles.</summary>
        public GeoCoordinate ToGeo(double x, double y, double z = 0)
            => ToGeo(new Vector(x, y, z));

        // ═══════════════════════════════════════════════════════════════
        //  Local Tangent Plane (ENU)
        // ═══════════════════════════════════════════════════════════════

        private Vector ToLocalTangentPlane(GeoCoordinate coord)
        {
            double dLat = coord.Latitude - Origin.Latitude;
            double dLon = coord.Longitude - Origin.Longitude;
            double x = dLon * _mPerDegLon;
            double y = dLat * _mPerDegLat;
            double z = coord.Altitude - Origin.Altitude;
            return new Vector(x, y, z);
        }

        private GeoCoordinate ToGeoTangentPlane(Vector local)
        {
            double lat = Origin.Latitude + local.y / _mPerDegLat;
            double lon = Origin.Longitude + local.x / _mPerDegLon;
            double alt = Origin.Altitude + local.z;
            return new GeoCoordinate(lat, lon, alt);
        }

        // ═══════════════════════════════════════════════════════════════
        //  UTM (Transverse Mercator)
        // ═══════════════════════════════════════════════════════════════

        private Vector ToLocalUtm(GeoCoordinate coord)
        {
            var (easting, northing) = LatLonToUtm(coord.Latitude, coord.Longitude);
            var (e0, n0) = LatLonToUtm(Origin.Latitude, Origin.Longitude);
            return new Vector(easting - e0, northing - n0, coord.Altitude - Origin.Altitude);
        }

        private GeoCoordinate ToGeoUtm(Vector local)
        {
            var (e0, n0) = LatLonToUtm(Origin.Latitude, Origin.Longitude);
            double easting = e0 + local.x;
            double northing = n0 + local.y;
            var (lat, lon) = UtmToLatLon(easting, northing);
            return new GeoCoordinate(lat, lon, Origin.Altitude + local.z);
        }

        // ─── UTM Forward ────────────────────────────────────────────

        private (double easting, double northing) LatLonToUtm(double lat, double lon)
        {
            double latRad = lat * Deg2Rad;
            double lonRad = lon * Deg2Rad;

            double lonOrigin = (_utmZone - 1) * 6 - 180 + 3; // centre meridian
            double lonOriginRad = lonOrigin * Deg2Rad;

            double ep2 = E2 / (1 - E2); // e'^2
            double N = A / Math.Sqrt(1 - E2 * Math.Sin(latRad) * Math.Sin(latRad));
            double T = Math.Tan(latRad) * Math.Tan(latRad);
            double C = ep2 * Math.Cos(latRad) * Math.Cos(latRad);
            double AA = Math.Cos(latRad) * (lonRad - lonOriginRad);

            double M = MeridianArc(latRad);

            const double k0 = 0.9996;

            double easting = k0 * N * (AA
                + (1 - T + C) * AA * AA * AA / 6
                + (5 - 18 * T + T * T + 72 * C - 58 * ep2) * AA * AA * AA * AA * AA / 120)
                + 500000.0;

            double northing = k0 * (M
                + N * Math.Tan(latRad) * (
                    AA * AA / 2
                    + (5 - T + 9 * C + 4 * C * C) * AA * AA * AA * AA / 24
                    + (61 - 58 * T + T * T + 600 * C - 330 * ep2) * Math.Pow(AA, 6) / 720));

            if (lat < 0)
                northing += 10000000.0;

            return (easting, northing);
        }

        // ─── UTM Inverse ────────────────────────────────────────────

        private (double lat, double lon) UtmToLatLon(double easting, double northing)
        {
            const double k0 = 0.9996;
            double ep2 = E2 / (1 - E2);

            double x = easting - 500000.0;
            double y = _utmNorth ? northing : northing - 10000000.0;

            double M = y / k0;
            double mu = M / (A * (1 - E2 / 4 - 3 * E2 * E2 / 64 - 5 * E2 * E2 * E2 / 256));

            double e1 = (1 - Math.Sqrt(1 - E2)) / (1 + Math.Sqrt(1 - E2));
            double phi1 = mu
                + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * Math.Sin(2 * mu)
                + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * Math.Sin(4 * mu)
                + (151 * e1 * e1 * e1 / 96) * Math.Sin(6 * mu);

            double N1 = A / Math.Sqrt(1 - E2 * Math.Sin(phi1) * Math.Sin(phi1));
            double T1 = Math.Tan(phi1) * Math.Tan(phi1);
            double C1 = ep2 * Math.Cos(phi1) * Math.Cos(phi1);
            double R1 = A * (1 - E2) / Math.Pow(1 - E2 * Math.Sin(phi1) * Math.Sin(phi1), 1.5);
            double D = x / (N1 * k0);

            double lat = phi1 - (N1 * Math.Tan(phi1) / R1) * (
                D * D / 2
                - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * ep2) * D * D * D * D / 24
                + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * ep2 - 3 * C1 * C1) * Math.Pow(D, 6) / 720);

            double lonOrigin = (_utmZone - 1) * 6 - 180 + 3;
            double lon = lonOrigin + (
                D
                - (1 + 2 * T1 + C1) * D * D * D / 6
                + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * ep2 + 24 * T1 * T1) * Math.Pow(D, 5) / 120)
                / Math.Cos(phi1) * Rad2Deg;

            return (lat * Rad2Deg, lon);
        }

        // ─── Meridian arc length ────────────────────────────────────

        private static double MeridianArc(double latRad)
        {
            return A * (
                (1 - E2 / 4 - 3 * E2 * E2 / 64 - 5 * E2 * E2 * E2 / 256) * latRad
                - (3 * E2 / 8 + 3 * E2 * E2 / 32 + 45 * E2 * E2 * E2 / 1024) * Math.Sin(2 * latRad)
                + (15 * E2 * E2 / 256 + 45 * E2 * E2 * E2 / 1024) * Math.Sin(4 * latRad)
                - (35 * E2 * E2 * E2 / 3072) * Math.Sin(6 * latRad));
        }
    }

    /// <summary>
    /// Specifies the projection method used to convert geographic to local coordinates.
    /// </summary>
    public enum ProjectionType
    {
        /// <summary>
        /// Local tangent plane (East-North-Up). Accurate for small areas (~50 km).
        /// </summary>
        LocalTangentPlane,

        /// <summary>
        /// Universal Transverse Mercator. More accurate over larger areas.
        /// </summary>
        UTM
    }
}
