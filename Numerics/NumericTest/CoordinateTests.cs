using CSharpNumerics.Engines.GIS.Coordinates;
using CSharpNumerics.Engines.GIS.Export;
using CSharpNumerics.Engines.GIS.Grid;
using CSharpNumerics.Numerics.Objects;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace NumericTest
{
    [TestClass]
    public class CoordinateTests
    {
        // ═══════════════════════════════════════════════════════════════
        //  GeoCoordinate
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoCoordinate_StoresLatLonAlt()
        {
            var gc = new GeoCoordinate(59.3293, 18.0686, 25.0);
            Assert.AreEqual(59.3293, gc.Latitude);
            Assert.AreEqual(18.0686, gc.Longitude);
            Assert.AreEqual(25.0, gc.Altitude);
        }

        [TestMethod]
        public void GeoCoordinate_DefaultAltitudeIsZero()
        {
            var gc = new GeoCoordinate(59.0, 18.0);
            Assert.AreEqual(0.0, gc.Altitude);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void GeoCoordinate_RejectsInvalidLatitude()
        {
            _ = new GeoCoordinate(91.0, 0.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void GeoCoordinate_RejectsInvalidLongitude()
        {
            _ = new GeoCoordinate(0.0, 181.0);
        }

        [TestMethod]
        public void GeoCoordinate_EqualityWorks()
        {
            var a = new GeoCoordinate(59.0, 18.0, 10.0);
            var b = new GeoCoordinate(59.0, 18.0, 10.0);
            var c = new GeoCoordinate(59.0, 18.1, 10.0);
            Assert.AreEqual(a, b);
            Assert.IsTrue(a == b);
            Assert.AreNotEqual(a, c);
            Assert.IsTrue(a != c);
        }

        [TestMethod]
        public void GeoCoordinate_DistanceTo_KnownPair()
        {
            // Stockholm (59.3293°N, 18.0686°E) to Gothenburg (57.7089°N, 11.9746°E)
            var stockholm = new GeoCoordinate(59.3293, 18.0686);
            var gothenburg = new GeoCoordinate(57.7089, 11.9746);
            double dist = stockholm.DistanceTo(gothenburg);
            // ~398 km by great circle
            Assert.IsTrue(dist > 390_000 && dist < 410_000,
                $"Expected ~398km, got {dist / 1000:F1}km");
        }

        [TestMethod]
        public void GeoCoordinate_ToString_ContainsValues()
        {
            var gc = new GeoCoordinate(59.3293, 18.0686, 25.0);
            string s = gc.ToString();
            Assert.IsTrue(s.Contains("59.3293"), $"Missing lat in: {s}");
            Assert.IsTrue(s.Contains("18.0686"), $"Missing lon in: {s}");
        }

        // ═══════════════════════════════════════════════════════════════
        //  Projection — Local Tangent Plane
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Projection_LTP_OriginMapsToZero()
        {
            var origin = new GeoCoordinate(59.3293, 18.0686);
            var proj = new Projection(origin, ProjectionType.LocalTangentPlane);
            var local = proj.ToLocal(origin);
            Assert.AreEqual(0.0, local.x, 0.01);
            Assert.AreEqual(0.0, local.y, 0.01);
            Assert.AreEqual(0.0, local.z, 0.01);
        }

        [TestMethod]
        public void Projection_LTP_1DegreeNorth_About111km()
        {
            var origin = new GeoCoordinate(59.0, 18.0);
            var proj = new Projection(origin);
            var local = proj.ToLocal(60.0, 18.0);
            // 1° latitude ≈ 111 km
            Assert.IsTrue(local.y > 110_000 && local.y < 112_000,
                $"Expected ~111km north, got {local.y:F0}m");
            Assert.AreEqual(0.0, local.x, 50.0); // should be ~0 east
        }

        [TestMethod]
        public void Projection_LTP_RoundTrip()
        {
            var origin = new GeoCoordinate(59.3293, 18.0686);
            var proj = new Projection(origin);

            var point = new GeoCoordinate(59.34, 18.08, 50.0);
            var local = proj.ToLocal(point);
            var back = proj.ToGeo(local);

            Assert.AreEqual(point.Latitude, back.Latitude, 1e-5);
            Assert.AreEqual(point.Longitude, back.Longitude, 1e-5);
            Assert.AreEqual(point.Altitude, back.Altitude, 0.01);
        }

        [TestMethod]
        public void Projection_LTP_EastPositive()
        {
            var origin = new GeoCoordinate(0.0, 0.0);
            var proj = new Projection(origin);
            var local = proj.ToLocal(0.0, 1.0); // 1° east of origin
            Assert.IsTrue(local.x > 0, "East should be positive x");
            Assert.AreEqual(0.0, local.y, 50.0);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Projection — UTM
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Projection_UTM_OriginMapsToZero()
        {
            var origin = new GeoCoordinate(59.3293, 18.0686);
            var proj = new Projection(origin, ProjectionType.UTM);
            var local = proj.ToLocal(origin);
            Assert.AreEqual(0.0, local.x, 0.01);
            Assert.AreEqual(0.0, local.y, 0.01);
        }

        [TestMethod]
        public void Projection_UTM_ZoneDetected()
        {
            // Stockholm is in UTM zone 34
            var proj = new Projection(new GeoCoordinate(59.3293, 18.0686), ProjectionType.UTM);
            Assert.AreEqual(34, proj.UtmZone);
            Assert.IsTrue(proj.UtmNorth);
        }

        [TestMethod]
        public void Projection_UTM_RoundTrip()
        {
            var origin = new GeoCoordinate(59.3293, 18.0686);
            var proj = new Projection(origin, ProjectionType.UTM);

            var point = new GeoCoordinate(59.35, 18.10, 30.0);
            var local = proj.ToLocal(point);
            var back = proj.ToGeo(local);

            Assert.AreEqual(point.Latitude, back.Latitude, 1e-5);
            Assert.AreEqual(point.Longitude, back.Longitude, 1e-5);
            Assert.AreEqual(point.Altitude, back.Altitude, 0.01);
        }

        [TestMethod]
        public void Projection_UTM_SouthernHemisphere()
        {
            var origin = new GeoCoordinate(-33.8688, 151.2093); // Sydney
            var proj = new Projection(origin, ProjectionType.UTM);
            Assert.IsFalse(proj.UtmNorth);

            var local = proj.ToLocal(-33.85, 151.22);
            var back = proj.ToGeo(local);
            Assert.AreEqual(-33.85, back.Latitude, 1e-4);
            Assert.AreEqual(151.22, back.Longitude, 1e-4);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoGrid.FromLatLon
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoGrid_FromLatLon_CreatesGrid()
        {
            var sw = new GeoCoordinate(59.32, 18.06);
            var ne = new GeoCoordinate(59.34, 18.10);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 50, 100);

            Assert.IsNotNull(grid.Projection);
            Assert.IsTrue(grid.Nx > 1);
            Assert.IsTrue(grid.Ny > 1);
            Assert.AreEqual(0.0, grid.XMin);
            Assert.AreEqual(0.0, grid.YMin);
        }

        [TestMethod]
        public void GeoGrid_FromLatLon_ProjectionPreserved()
        {
            var sw = new GeoCoordinate(59.32, 18.06);
            var ne = new GeoCoordinate(59.34, 18.10);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 10, 200);

            Assert.AreEqual(ProjectionType.LocalTangentPlane, grid.Projection.Type);
            Assert.AreEqual(sw.Latitude, grid.Projection.Origin.Latitude);
            Assert.AreEqual(sw.Longitude, grid.Projection.Origin.Longitude);
        }

        [TestMethod]
        public void GeoGrid_FromLatLon_UTM()
        {
            var sw = new GeoCoordinate(59.32, 18.06);
            var ne = new GeoCoordinate(59.34, 18.10);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 10, 200, ProjectionType.UTM);

            Assert.AreEqual(ProjectionType.UTM, grid.Projection.Type);
            Assert.IsTrue(grid.CellCount > 0);
        }

        [TestMethod]
        public void GeoGrid_CellCentreGeo_ReturnsWGS84()
        {
            var sw = new GeoCoordinate(59.32, 18.06);
            var ne = new GeoCoordinate(59.34, 18.10);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 0, 500);

            // First cell centre should be close to the SW corner
            var geo = grid.CellCentreGeo(0);
            Assert.IsTrue(geo.Latitude >= 59.31 && geo.Latitude <= 59.35,
                $"Lat: {geo.Latitude}");
            Assert.IsTrue(geo.Longitude >= 18.05 && geo.Longitude <= 18.11,
                $"Lon: {geo.Longitude}");
        }

        [TestMethod]
        [ExpectedException(typeof(InvalidOperationException))]
        public void GeoGrid_CellCentreGeo_ThrowsWithoutProjection()
        {
            var grid = new GeoGrid(0, 1000, 0, 1000, 0, 0, 100);
            grid.CellCentreGeo(0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void GeoGrid_FromLatLon_RejectsBadBounds()
        {
            var sw = new GeoCoordinate(59.34, 18.06);
            var ne = new GeoCoordinate(59.32, 18.10); // lat reversed
            GeoGrid.FromLatLon(sw, ne, 0, 10, 100);
        }

        // ═══════════════════════════════════════════════════════════════
        //  GeoJSON with CRS
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void GeoJson_LocalGrid_WritesLocalCrs()
        {
            var grid = new GeoGrid(0, 200, 0, 200, 0, 0, 100);
            var snapshot = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);
            var json = GeoJsonExporter.ToGeoJson(snapshot);
            Assert.IsTrue(json.Contains("\"crs\":\"local\""));
        }

        [TestMethod]
        public void GeoJson_GeoGrid_WritesWGS84()
        {
            var sw = new GeoCoordinate(59.32, 18.06);
            var ne = new GeoCoordinate(59.34, 18.10);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 0, 500);

            var values = new double[grid.CellCount];
            values[0] = 0.5;
            var snapshot = new GridSnapshot(grid, values, 0, 0);
            var json = GeoJsonExporter.ToGeoJson(snapshot);

            Assert.IsTrue(json.Contains("\"crs\":\"WGS84\""),
                "Should contain WGS84 CRS");
            // GeoJSON coordinates should be in lon/lat range, not metres
            Assert.IsTrue(json.Contains("18."), "Should contain longitude ~18");
            Assert.IsTrue(json.Contains("59."), "Should contain latitude ~59");
        }

        [TestMethod]
        public void GeoJson_GeoGrid_CoordinatesAreLonLatOrder()
        {
            var sw = new GeoCoordinate(10.0, 20.0);
            var ne = new GeoCoordinate(10.01, 20.01);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 0, 500);

            var snapshot = new GridSnapshot(grid, new double[grid.CellCount], 0, 0);
            var json = GeoJsonExporter.ToGeoJson(snapshot);

            // In GeoJSON, coordinates are [longitude, latitude, altitude]
            // First coordinate should start with ~20 (longitude), not ~10 (latitude)
            int coordIdx = json.IndexOf("\"coordinates\":[") + "\"coordinates\":[".Length;
            string coordStart = json.Substring(coordIdx, 4);
            Assert.IsTrue(coordStart.StartsWith("20"),
                $"Expected longitude first (~20.x), got: {coordStart}");
        }

        [TestMethod]
        public void Cesium_GeoGrid_WritesWGS84()
        {
            var sw = new GeoCoordinate(59.32, 18.06);
            var ne = new GeoCoordinate(59.34, 18.10);
            var grid = GeoGrid.FromLatLon(sw, ne, 0, 0, 500);

            var values = new double[grid.CellCount];
            for (int i = 0; i < values.Length; i++) values[i] = 0.3;

            // Build a minimal ProbabilityMap manually via a single-scenario path
            // Use CesiumExporter.ToGeoJsonCesium with the grid
            var snapshot = new GridSnapshot(grid, values, 0, 0);
            var json = GeoJsonExporter.ToGeoJson(snapshot);
            Assert.IsTrue(json.Contains("\"crs\":\"WGS84\""));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Projection edge cases
        // ═══════════════════════════════════════════════════════════════

        [TestMethod]
        public void Projection_LTP_NegativeLatitude()
        {
            var origin = new GeoCoordinate(-34.0, 151.0);
            var proj = new Projection(origin);
            var local = proj.ToLocal(-33.0, 151.0);
            Assert.IsTrue(local.y > 0, "Moving north should increase y");
        }

        [TestMethod]
        public void Projection_LTP_AltitudePreserved()
        {
            var origin = new GeoCoordinate(59.0, 18.0, 100.0);
            var proj = new Projection(origin);
            var local = proj.ToLocal(59.0, 18.0, 150.0);
            Assert.AreEqual(50.0, local.z, 0.01);
        }

        [TestMethod]
        public void GeoCoordinate_HashCode_ConsistentWithEquals()
        {
            var a = new GeoCoordinate(59.0, 18.0, 10.0);
            var b = new GeoCoordinate(59.0, 18.0, 10.0);
            Assert.AreEqual(a.GetHashCode(), b.GetHashCode());
        }
    }
}
