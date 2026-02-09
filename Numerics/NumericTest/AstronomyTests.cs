using CSharpNumerics.Physics;

namespace NumericTest
{
    [TestClass]
    public class AstronomyTests
    {
        private const double Tolerance = 1e-4;

        #region Distance Conversions

        [TestMethod]
        public void LightYearsToParsecs_OneLightYear()
        {
            double pc = 1.0.LightYearsToParsecs();
            // 1 ly ≈ 0.3066 pc
            Assert.AreEqual(0.3066, pc, 0.001);
        }

        [TestMethod]
        public void ParsecsToLightYears_OneParsec()
        {
            double ly = 1.0.ParsecsToLightYears();
            // 1 pc ≈ 3.2616 ly
            Assert.AreEqual(3.2616, ly, 0.001);
        }

        [TestMethod]
        public void LightYearsToParsecs_RoundTrip()
        {
            double original = 4.37; // Proxima Centauri distance
            double roundTrip = original.LightYearsToParsecs().ParsecsToLightYears();
            Assert.AreEqual(original, roundTrip, 1e-10);
        }

        [TestMethod]
        public void LightYearsToAU_OneLightYear()
        {
            double au = 1.0.LightYearsToAU();
            // 1 ly ≈ 63241 AU
            Assert.AreEqual(63241, au, 10);
        }

        [TestMethod]
        public void ParsecsToAU_OneParsec()
        {
            double au = 1.0.ParsecsToAU();
            // 1 pc ≈ 206265 AU
            Assert.AreEqual(206265, au, 10);
        }

        #endregion

        #region Julian Date

        [TestMethod]
        public void JulianDate_J2000_Epoch()
        {
            // J2000.0 = 2000-01-01 12:00:00 UTC → JD 2451545.0
            var j2000 = new DateTime(2000, 1, 1, 12, 0, 0, DateTimeKind.Utc);
            double jd = j2000.ToJulianDate();
            Assert.AreEqual(2451545.0, jd, 0.001);
        }

        [TestMethod]
        public void JulianDate_KnownDate()
        {
            // 2024-03-20 00:00 UTC → JD 2460389.5
            var date = new DateTime(2024, 3, 20, 0, 0, 0, DateTimeKind.Utc);
            double jd = date.ToJulianDate();
            Assert.AreEqual(2460389.5, jd, 0.001);
        }

        [TestMethod]
        public void JulianCenturies_J2000_IsZero()
        {
            var j2000 = new DateTime(2000, 1, 1, 12, 0, 0, DateTimeKind.Utc);
            double T = j2000.JulianCenturiesSinceJ2000();
            Assert.AreEqual(0.0, T, 1e-10);
        }

        #endregion

        #region Sidereal Time

        [TestMethod]
        public void GMST_J2000_Epoch()
        {
            // At J2000.0, GMST ≈ 280.46° ≈ 18.697 hours
            var j2000 = new DateTime(2000, 1, 1, 12, 0, 0, DateTimeKind.Utc);
            double gmstDeg = j2000.GreenwichMeanSiderealTimeDegrees();
            Assert.AreEqual(280.46, gmstDeg, 0.1);
        }

        [TestMethod]
        public void GMST_Hours_MatchesDegrees()
        {
            var utc = new DateTime(2024, 6, 15, 21, 0, 0, DateTimeKind.Utc);
            double deg = utc.GreenwichMeanSiderealTimeDegrees();
            double hours = utc.GreenwichMeanSiderealTimeHours();
            Assert.AreEqual(deg / 15.0, hours, 1e-10);
        }

        [TestMethod]
        public void LMST_GreenwichEqualsGMST()
        {
            // At longitude 0° (Greenwich), LMST = GMST
            var utc = new DateTime(2024, 6, 15, 21, 0, 0, DateTimeKind.Utc);
            double gmst = utc.GreenwichMeanSiderealTimeDegrees();
            double lmst = utc.LocalMeanSiderealTimeDegrees(0.0);
            Assert.AreEqual(gmst, lmst, 1e-10);
        }

        [TestMethod]
        public void LMST_EastLongitude_IsGreaterThanGMST()
        {
            var utc = new DateTime(2024, 6, 15, 21, 0, 0, DateTimeKind.Utc);
            double gmst = utc.GreenwichMeanSiderealTimeDegrees();
            // Stockholm ≈ 18.07° E
            double lmst = utc.LocalMeanSiderealTimeDegrees(18.07);
            double expected = (gmst + 18.07) % 360.0;
            Assert.AreEqual(expected, lmst, 1e-10);
        }

        [TestMethod]
        public void LocalSiderealTime_FromLocalTime()
        {
            // Stockholm: 2024-06-15 23:00 local, UTC+2, longitude 18.07°E
            var local = new DateTime(2024, 6, 15, 23, 0, 0);
            double utcOffset = 2.0;
            double longitude = 18.07;

            double lstHours = local.LocalSiderealTimeFromLocal(utcOffset, longitude);

            // Verify same as converting to UTC and calling LMST directly
            var utc = new DateTime(2024, 6, 15, 21, 0, 0, DateTimeKind.Utc);
            double expected = utc.LocalMeanSiderealTimeHours(longitude);
            Assert.AreEqual(expected, lstHours, 1e-10);
        }

        #endregion

        #region Horizontal ↔ Equatorial Transforms

        [TestMethod]
        public void HorizontalToEquatorial_Zenith_GivesDeclEqualLatitude()
        {
            // A star at zenith (alt=90°) must have dec = lat
            double lat = 52.0; // Berlin
            double lst = 90.0;
            var (ra, dec) = AstronomyExtensions.HorizontalToEquatorial(90.0, 0.0, lat, lst);
            Assert.AreEqual(lat, dec, 0.01);
        }

        [TestMethod]
        public void EquatorialToHorizontal_RoundTrip()
        {
            double lat = 59.33; // Stockholm
            double lst = 120.0;
            double raIn = 200.0;
            double decIn = 30.0;

            var (alt, az) = AstronomyExtensions.EquatorialToHorizontal(raIn, decIn, lat, lst);
            var (raOut, decOut) = AstronomyExtensions.HorizontalToEquatorial(alt, az, lat, lst);

            Assert.AreEqual(raIn, raOut, 0.01);
            Assert.AreEqual(decIn, decOut, 0.01);
        }

        [TestMethod]
        public void HorizontalToEquatorial_WithUtc()
        {
            double alt = 45.0;
            double az = 180.0; // due south
            double lat = 51.48; // London
            double lon = -0.0;  // Greenwich
            var utc = new DateTime(2024, 3, 20, 22, 0, 0, DateTimeKind.Utc);

            var (ra1, dec1) = AstronomyExtensions.HorizontalToEquatorial(alt, az, lat, lon, utc);

            // Same via manual LST
            double lst = utc.LocalMeanSiderealTimeDegrees(lon);
            var (ra2, dec2) = AstronomyExtensions.HorizontalToEquatorial(alt, az, lat, lst);

            Assert.AreEqual(ra1, ra2, 1e-10);
            Assert.AreEqual(dec1, dec2, 1e-10);
        }

        [TestMethod]
        public void EquatorialToHorizontal_WithUtc()
        {
            double ra = 100.0;
            double dec = 20.0;
            double lat = 51.48;
            double lon = -0.0;
            var utc = new DateTime(2024, 3, 20, 22, 0, 0, DateTimeKind.Utc);

            var (alt1, az1) = AstronomyExtensions.EquatorialToHorizontal(ra, dec, lat, lon, utc);

            double lst = utc.LocalMeanSiderealTimeDegrees(lon);
            var (alt2, az2) = AstronomyExtensions.EquatorialToHorizontal(ra, dec, lat, lst);

            Assert.AreEqual(alt1, alt2, 1e-10);
            Assert.AreEqual(az1, az2, 1e-10);
        }

        [TestMethod]
        public void EquatorialToHorizontal_RoundTrip_MultiplePositions()
        {
            double lat = 45.0;
            double lst = 200.0;

            double[] ras = { 0, 45, 90, 180, 270, 359.9 };
            double[] decs = { -60, -30, 0, 30, 60, 89 };

            foreach (double ra in ras)
            {
                foreach (double dec in decs)
                {
                    var (alt, az) = AstronomyExtensions.EquatorialToHorizontal(ra, dec, lat, lst);
                    var (raBack, decBack) = AstronomyExtensions.HorizontalToEquatorial(alt, az, lat, lst);

                    Assert.AreEqual(ra, raBack, 0.1, $"RA failed for ra={ra}, dec={dec}");
                    Assert.AreEqual(dec, decBack, 0.1, $"Dec failed for ra={ra}, dec={dec}");
                }
            }
        }

        #endregion

        #region Angle Helpers

        [TestMethod]
        public void RightAscensionToDegrees_6h30m()
        {
            // 6h 30m 0s = 97.5°
            double deg = AstronomyExtensions.RightAscensionToDegrees(6, 30, 0);
            Assert.AreEqual(97.5, deg, 1e-10);
        }

        [TestMethod]
        public void DegreesToRightAscension_97_5_Degrees()
        {
            var (h, m, s) = AstronomyExtensions.DegreesToRightAscension(97.5);
            Assert.AreEqual(6, h);
            Assert.AreEqual(30, m);
            Assert.AreEqual(0.0, s, 0.01);
        }

        [TestMethod]
        public void DeclinationToDegrees_Negative()
        {
            // -16° 42' 58" ≈ -16.7161°
            double dec = AstronomyExtensions.DeclinationToDegrees(-16, 42, 58);
            Assert.AreEqual(-16.7161, dec, 0.001);
        }

        [TestMethod]
        public void RightAscension_RoundTrip()
        {
            double original = 213.75;
            var (h, m, s) = AstronomyExtensions.DegreesToRightAscension(original);
            double roundTrip = AstronomyExtensions.RightAscensionToDegrees(h, m, s);
            Assert.AreEqual(original, roundTrip, 0.001);
        }

        #endregion
    }
}
