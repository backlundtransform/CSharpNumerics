using CSharpNumerics.Physics.Astro.Enums;
using CSharpNumerics.Physics.Constants;
using System;

namespace CSharpNumerics.Physics.Astro
{
    /// <summary>
    /// Provides extension methods for astronomical calculations:
    /// distance conversions, sidereal time, and equatorial coordinate transforms.
    /// All methods are self-contained — no external libraries are used.
    /// </summary>
    public static class AstronomyExtensions
    {
        #region Distance Conversions

        /// <summary>
        /// Converts a distance in light-years to parsecs: pc = ly / 3.26156.
        /// </summary>
        /// <param name="lightYears">Distance in light-years.</param>
        public static double LightYearsToParsecs(this double lightYears)
        {
            return lightYears * PhysicsConstants.LightYear / PhysicsConstants.Parsec;
        }

        /// <summary>
        /// Converts a distance in parsecs to light-years: ly = pc * 3.26156.
        /// </summary>
        /// <param name="parsecs">Distance in parsecs.</param>
        public static double ParsecsToLightYears(this double parsecs)
        {
            return parsecs * PhysicsConstants.Parsec / PhysicsConstants.LightYear;
        }

        /// <summary>
        /// Converts a distance in light-years to astronomical units: AU = ly * (LightYear / AU).
        /// </summary>
        /// <param name="lightYears">Distance in light-years.</param>
        public static double LightYearsToAU(this double lightYears)
        {
            return lightYears * PhysicsConstants.LightYear / PhysicsConstants.AstronomicalUnit;
        }

        /// <summary>
        /// Converts a distance in parsecs to astronomical units: AU = pc * (Parsec / AU).
        /// </summary>
        /// <param name="parsecs">Distance in parsecs.</param>
        public static double ParsecsToAU(this double parsecs)
        {
            return parsecs * PhysicsConstants.Parsec / PhysicsConstants.AstronomicalUnit;
        }

        #endregion

        #region Julian Date

        /// <summary>
        /// Computes the Julian Date for a given UTC DateTime.
        /// Uses the standard algorithm valid for dates after the Gregorian reform.
        /// </summary>
        /// <param name="utc">The UTC date and time.</param>
        public static double ToJulianDate(this DateTime utc)
        {
            int y = utc.Year;
            int m = utc.Month;
            int d = utc.Day;

            if (m <= 2)
            {
                y -= 1;
                m += 12;
            }

            int A = y / 100;
            int B = 2 - A + A / 4;

            double dayFraction = (utc.Hour + utc.Minute / 60.0 + utc.Second / 3600.0) / 24.0;

            return Math.Floor(365.25 * (y + 4716))
                 + Math.Floor(30.6001 * (m + 1))
                 + d + dayFraction + B - 1524.5;
        }

        /// <summary>
        /// Computes the number of Julian centuries since J2000.0 (2000-01-01 12:00 TT).
        /// </summary>
        /// <param name="utc">The UTC date and time.</param>
        public static double JulianCenturiesSinceJ2000(this DateTime utc)
        {
            double jd = utc.ToJulianDate();
            return (jd - 2451545.0) / 36525.0;
        }

        #endregion

        #region Sidereal Time

        /// <summary>
        /// Computes Greenwich Mean Sidereal Time (GMST) in degrees for a given UTC time.
        /// Based on the IAU formula using Julian centuries since J2000.0.
        /// </summary>
        /// <param name="utc">The UTC date and time.</param>
        public static double GreenwichMeanSiderealTimeDegrees(this DateTime utc)
        {
            double T = utc.JulianCenturiesSinceJ2000();

            // GMST at 0h UT1 in seconds of time, then add rotation for UT1 fraction
            double gmstDeg = 280.46061837
                           + 360.98564736629 * (utc.ToJulianDate() - 2451545.0)
                           + 0.000387933 * T * T
                           - T * T * T / 38710000.0;

            return NormalizeDegrees(gmstDeg);
        }

        /// <summary>
        /// Computes Greenwich Mean Sidereal Time (GMST) in hours for a given UTC time.
        /// </summary>
        /// <param name="utc">The UTC date and time.</param>
        public static double GreenwichMeanSiderealTimeHours(this DateTime utc)
        {
            return utc.GreenwichMeanSiderealTimeDegrees() / 15.0;
        }

        /// <summary>
        /// Computes Local Mean Sidereal Time (LMST) in degrees for a given UTC time and longitude.
        /// </summary>
        /// <param name="utc">The UTC date and time.</param>
        /// <param name="longitudeDegrees">Observer's longitude in degrees (east positive).</param>
        public static double LocalMeanSiderealTimeDegrees(this DateTime utc, double longitudeDegrees)
        {
            return NormalizeDegrees(utc.GreenwichMeanSiderealTimeDegrees() + longitudeDegrees);
        }

        /// <summary>
        /// Computes Local Mean Sidereal Time (LMST) in hours for a given UTC time and longitude.
        /// </summary>
        /// <param name="utc">The UTC date and time.</param>
        /// <param name="longitudeDegrees">Observer's longitude in degrees (east positive).</param>
        public static double LocalMeanSiderealTimeHours(this DateTime utc, double longitudeDegrees)
        {
            return utc.LocalMeanSiderealTimeDegrees(longitudeDegrees) / 15.0;
        }

        /// <summary>
        /// Computes Local Mean Sidereal Time from a local time, time zone offset, and longitude.
        /// The local time is first converted to UTC, then LMST is computed.
        /// </summary>
        /// <param name="localTime">The observer's local date and time.</param>
        /// <param name="utcOffsetHours">UTC offset in hours (e.g. +1 for CET, -5 for EST).</param>
        /// <param name="longitudeDegrees">Observer's longitude in degrees (east positive).</param>
        public static double LocalSiderealTimeFromLocal(
            this DateTime localTime,
            double utcOffsetHours,
            double longitudeDegrees)
        {
            DateTime utc = localTime.AddHours(-utcOffsetHours);
            return utc.LocalMeanSiderealTimeHours(longitudeDegrees);
        }

        #endregion

        #region Horizontal to Equatorial Coordinate Transform

        /// <summary>
        /// Converts horizontal coordinates (altitude, azimuth) to equatorial coordinates
        /// (right ascension, declination) given the observer's latitude and local sidereal time.
        /// Returns (rightAscensionDegrees, declinationDegrees).
        /// </summary>
        /// <param name="altitudeDegrees">Altitude above the horizon in degrees (-90 to +90).</param>
        /// <param name="azimuthDegrees">Azimuth in degrees measured from north through east (0 to 360).</param>
        /// <param name="latitudeDegrees">Observer's geographic latitude in degrees (-90 to +90).</param>
        /// <param name="localSiderealTimeDegrees">Local sidereal time in degrees (0 to 360).</param>
        public static (double RightAscensionDegrees, double DeclinationDegrees) HorizontalToEquatorial(
            double altitudeDegrees,
            double azimuthDegrees,
            double latitudeDegrees,
            double localSiderealTimeDegrees)
        {
            double alt = DegreesToRadians(altitudeDegrees);
            double az = DegreesToRadians(azimuthDegrees);
            double lat = DegreesToRadians(latitudeDegrees);

            // Declination: sin(dec) = sin(alt)*sin(lat) + cos(alt)*cos(lat)*cos(az)
            double sinDec = Math.Sin(alt) * Math.Sin(lat)
                          + Math.Cos(alt) * Math.Cos(lat) * Math.Cos(az);
            double dec = Math.Asin(Clamp(sinDec, -1.0, 1.0));

            // Hour angle: sin(H) = -sin(az)*cos(alt) / cos(dec)
            //             cos(H) = (sin(alt) - sin(lat)*sin(dec)) / (cos(lat)*cos(dec))
            double cosDecCosLat = Math.Cos(dec) * Math.Cos(lat);
            double cosH = Math.Abs(cosDecCosLat) < 1e-15
                ? 0.0
                : (Math.Sin(alt) - Math.Sin(lat) * Math.Sin(dec)) / cosDecCosLat;
            double sinH = Math.Abs(Math.Cos(dec)) < 1e-15
                ? 0.0
                : -Math.Sin(az) * Math.Cos(alt) / Math.Cos(dec);

            double H = RadiansToDegrees(Math.Atan2(sinH, cosH));

            // Right ascension = LST - H
            double ra = NormalizeDegrees(localSiderealTimeDegrees - H);

            return (ra, RadiansToDegrees(dec));
        }

        /// <summary>
        /// Converts horizontal coordinates to equatorial coordinates using UTC time and observer position.
        /// Returns (rightAscensionDegrees, declinationDegrees).
        /// </summary>
        /// <param name="altitudeDegrees">Altitude above the horizon in degrees.</param>
        /// <param name="azimuthDegrees">Azimuth in degrees measured from north through east.</param>
        /// <param name="latitudeDegrees">Observer's geographic latitude in degrees.</param>
        /// <param name="longitudeDegrees">Observer's geographic longitude in degrees (east positive).</param>
        /// <param name="utc">The UTC date and time of observation.</param>
        public static (double RightAscensionDegrees, double DeclinationDegrees) HorizontalToEquatorial(
            double altitudeDegrees,
            double azimuthDegrees,
            double latitudeDegrees,
            double longitudeDegrees,
            DateTime utc)
        {
            double lst = utc.LocalMeanSiderealTimeDegrees(longitudeDegrees);
            return HorizontalToEquatorial(altitudeDegrees, azimuthDegrees, latitudeDegrees, lst);
        }

        #endregion

        #region Equatorial to Horizontal Coordinate Transform

        /// <summary>
        /// Converts equatorial coordinates (right ascension, declination) to horizontal coordinates
        /// (altitude, azimuth) given the observer's latitude and local sidereal time.
        /// Returns (altitudeDegrees, azimuthDegrees).
        /// </summary>
        /// <param name="rightAscensionDegrees">Right ascension in degrees (0 to 360).</param>
        /// <param name="declinationDegrees">Declination in degrees (-90 to +90).</param>
        /// <param name="latitudeDegrees">Observer's geographic latitude in degrees (-90 to +90).</param>
        /// <param name="localSiderealTimeDegrees">Local sidereal time in degrees (0 to 360).</param>
        public static (double AltitudeDegrees, double AzimuthDegrees) EquatorialToHorizontal(
            double rightAscensionDegrees,
            double declinationDegrees,
            double latitudeDegrees,
            double localSiderealTimeDegrees)
        {
            double H = DegreesToRadians(NormalizeDegrees(localSiderealTimeDegrees - rightAscensionDegrees));
            double dec = DegreesToRadians(declinationDegrees);
            double lat = DegreesToRadians(latitudeDegrees);

            // Altitude: sin(alt) = sin(dec)*sin(lat) + cos(dec)*cos(lat)*cos(H)
            double sinAlt = Math.Sin(dec) * Math.Sin(lat)
                          + Math.Cos(dec) * Math.Cos(lat) * Math.Cos(H);
            double alt = Math.Asin(Clamp(sinAlt, -1.0, 1.0));

            // Azimuth: sin(A) = -sin(H)*cos(dec) / cos(alt)
            //          cos(A) = (sin(dec) - sin(lat)*sin(alt)) / (cos(lat)*cos(alt))
            double cosLatCosAlt = Math.Cos(lat) * Math.Cos(alt);
            double cosA = Math.Abs(cosLatCosAlt) < 1e-15
                ? 0.0
                : (Math.Sin(dec) - Math.Sin(lat) * Math.Sin(alt)) / cosLatCosAlt;
            double sinA = Math.Abs(Math.Cos(alt)) < 1e-15
                ? 0.0
                : -Math.Sin(H) * Math.Cos(dec) / Math.Cos(alt);

            double az = NormalizeDegrees(RadiansToDegrees(Math.Atan2(sinA, cosA)));

            return (RadiansToDegrees(alt), az);
        }

        /// <summary>
        /// Converts equatorial coordinates to horizontal coordinates using UTC time and observer position.
        /// Returns (altitudeDegrees, azimuthDegrees).
        /// </summary>
        /// <param name="rightAscensionDegrees">Right ascension in degrees.</param>
        /// <param name="declinationDegrees">Declination in degrees.</param>
        /// <param name="latitudeDegrees">Observer's geographic latitude in degrees.</param>
        /// <param name="longitudeDegrees">Observer's geographic longitude in degrees (east positive).</param>
        /// <param name="utc">The UTC date and time of observation.</param>
        public static (double AltitudeDegrees, double AzimuthDegrees) EquatorialToHorizontal(
            double rightAscensionDegrees,
            double declinationDegrees,
            double latitudeDegrees,
            double longitudeDegrees,
            DateTime utc)
        {
            double lst = utc.LocalMeanSiderealTimeDegrees(longitudeDegrees);
            return EquatorialToHorizontal(rightAscensionDegrees, declinationDegrees, latitudeDegrees, lst);
        }

        #endregion

        #region Angle Helpers

        /// <summary>
        /// Converts right ascension from hours, minutes, seconds to degrees.
        /// </summary>
        public static double RightAscensionToDegrees(double hours, double minutes = 0, double seconds = 0)
        {
            return (hours + minutes / 60.0 + seconds / 3600.0) * 15.0;
        }

        /// <summary>
        /// Converts right ascension in degrees to (hours, minutes, seconds).
        /// </summary>
        public static (int Hours, int Minutes, double Seconds) DegreesToRightAscension(double degrees)
        {
            double totalHours = NormalizeDegrees(degrees) / 15.0;
            int h = (int)totalHours;
            double remainMinutes = (totalHours - h) * 60.0;
            int m = (int)remainMinutes;
            double s = (remainMinutes - m) * 60.0;
            return (h, m, s);
        }

        /// <summary>
        /// Converts declination from degrees, arcminutes, arcseconds to decimal degrees.
        /// </summary>
        public static double DeclinationToDegrees(double degrees, double arcminutes = 0, double arcseconds = 0)
        {
            double sign = degrees < 0 ? -1.0 : 1.0;
            return sign * (Math.Abs(degrees) + arcminutes / 60.0 + arcseconds / 3600.0);
        }

        private static double DegreesToRadians(double deg) => deg * Math.PI / 180.0;

        private static double RadiansToDegrees(double rad) => rad * 180.0 / Math.PI;

        private static double NormalizeDegrees(double deg)
        {
            deg %= 360.0;
            if (deg < 0) deg += 360.0;
            return deg;
        }

        private static double Clamp(double value, double min, double max)
        {
            if (value < min) return min;
            if (value > max) return max;
            return value;
        }

        #endregion

        #region Exoplanet Classification

        /// <summary>
        /// Returns the Harvard spectral type for a given stellar effective temperature (Kelvin).
        /// Boundaries follow the standard Morgan–Keenan classification.
        /// </summary>
        /// <param name="effectiveTemperatureK">Stellar effective temperature in Kelvin.</param>
        public static SpectralType GetSpectralFromTemp(double effectiveTemperatureK)
        {
            if (effectiveTemperatureK <= 0) return SpectralType.Unknown;

            if (effectiveTemperatureK >= 30000) return SpectralType.O;
            if (effectiveTemperatureK >= 10000) return SpectralType.B;
            if (effectiveTemperatureK >= 7500)  return SpectralType.A;
            if (effectiveTemperatureK >= 6000)  return SpectralType.F;
            if (effectiveTemperatureK >= 5200)  return SpectralType.G;
            if (effectiveTemperatureK >= 3700)  return SpectralType.K;
            if (effectiveTemperatureK >= 2400)  return SpectralType.M;
            if (effectiveTemperatureK >= 1300)  return SpectralType.L;
            if (effectiveTemperatureK >= 550)   return SpectralType.T;

            return SpectralType.Y;
        }

        /// <summary>
        /// Calculates the conservative habitable zone (Goldilocks zone) boundaries in AU
        /// using the Kopparapu et al. (2013) parameterisation.
        /// Returns (innerEdgeAU, outerEdgeAU).
        /// </summary>
        /// <param name="stellarLuminositySolar">Stellar luminosity in solar luminosities.</param>
        /// <param name="effectiveTemperatureK">Stellar effective temperature in Kelvin.</param>
        public static (double InnerEdgeAU, double OuterEdgeAU) CalculateGoldilocksZone(
            double stellarLuminositySolar,
            double effectiveTemperatureK)
        {
            double tStar = effectiveTemperatureK - 5780.0;
            double tStar2 = tStar * tStar;
            double tStar3 = tStar2 * tStar;
            double tStar4 = tStar3 * tStar;

            // Runaway greenhouse (inner edge, conservative)
            double sInner = 1.0146
                          + 8.1884e-5 * tStar
                          + 1.9394e-9 * tStar2
                          - 4.3618e-12 * tStar3
                          - 6.8260e-16 * tStar4;

            // Maximum greenhouse (outer edge, conservative)
            double sOuter = 0.3507
                          + 5.9578e-5 * tStar
                          + 1.6707e-9 * tStar2
                          - 3.0058e-12 * tStar3
                          - 5.1925e-16 * tStar4;

            double inner = Math.Sqrt(stellarLuminositySolar / sInner);
            double outer = Math.Sqrt(stellarLuminositySolar / sOuter);

            return (inner, outer);
        }

        /// <summary>
        /// Calculates the conservative habitable zone boundaries in AU from stellar radius and temperature.
        /// Luminosity is derived via the Stefan–Boltzmann law: L = (R/R☉)² × (T/T☉)⁴.
        /// Returns (innerEdgeAU, outerEdgeAU).
        /// </summary>
        /// <param name="stellarRadiusSolar">Stellar radius in solar radii.</param>
        /// <param name="effectiveTemperatureK">Stellar effective temperature in Kelvin.</param>
        public static (double InnerEdgeAU, double OuterEdgeAU) CalculateGoldilocksZoneFromRadius(
            double stellarRadiusSolar,
            double effectiveTemperatureK)
        {
            const double solarTeff = 5778.0;
            double luminosity = stellarRadiusSolar * stellarRadiusSolar
                              * Math.Pow(effectiveTemperatureK / solarTeff, 4);
            return CalculateGoldilocksZone(luminosity, effectiveTemperatureK);
        }

        /// <summary>
        /// Computes the Earth Similarity Index (ESI), a value between 0 and 1
        /// that measures how physically similar a planet is to Earth.
        /// Uses the standard four-parameter formulation (Schulze-Makuch et al. 2011).
        /// </summary>
        /// <param name="radiusEarth">Planet radius in Earth radii (Earth = 1.0).</param>
        /// <param name="densityEarth">Planet bulk density in Earth densities (Earth = 1.0).</param>
        /// <param name="escapeVelocityEarth">Planet escape velocity in Earth escape velocities (Earth = 1.0).</param>
        /// <param name="surfaceTemperatureK">Planet surface temperature in Kelvin (Earth ≈ 288 K).</param>
        public static double CalculateEsi(
            double radiusEarth,
            double densityEarth,
            double escapeVelocityEarth,
            double surfaceTemperatureK)
        {
            const int n = 4;
            const double refTemp = 288.0;
            const double wRadius = 0.57;
            const double wDensity = 1.07;
            const double wEscapeVelocity = 0.70;
            const double wTemperature = 5.58;

            double EsiTerm(double value, double reference, double weight)
            {
                double ratio = Math.Abs(value - reference) / (value + reference);
                return Math.Pow(1.0 - ratio, weight / n);
            }

            return EsiTerm(radiusEarth, 1.0, wRadius)
                 * EsiTerm(densityEarth, 1.0, wDensity)
                 * EsiTerm(escapeVelocityEarth, 1.0, wEscapeVelocity)
                 * EsiTerm(surfaceTemperatureK, refTemp, wTemperature);
        }

        #endregion
    }
}
