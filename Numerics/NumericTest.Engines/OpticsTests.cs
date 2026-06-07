using CSharpNumerics.Engines.Game;
using CSharpNumerics.Engines.Game.Objects;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Materials;
using CSharpNumerics.Physics.Materials.Optical;
using CSharpNumerics.Physics.Optics;
using CSharpNumerics.Physics.Optics.Enums;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;

namespace NumericTest;

[TestClass]
public class OpticsTests
{
    private const double Tol = 1e-6;

    // ══════════════════════════════════════════════════════
    //  Ray
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void Ray_Direction_Is_Normalised()
    {
        var ray = new Ray(new Vector(0, 0, 0), new Vector(3, 4, 0));
        Assert.AreEqual(1.0, ray.Direction.GetMagnitude(), Tol);
    }

    [TestMethod]
    public void Ray_PointAt_Returns_Correct_Position()
    {
        var ray = new Ray(new Vector(1, 0, 0), new Vector(1, 0, 0));
        var pt = ray.PointAt(5);
        Assert.AreEqual(6.0, pt.x, Tol);
        Assert.AreEqual(0.0, pt.y, Tol);
    }

    // ══════════════════════════════════════════════════════
    //  OpticalMedium & Dispersion
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void OpticalMedium_Dispersion_Returns_BaseIndex_When_No_AbbeNumber()
    {
        var m = new OpticalMedium("Test", 1.5);
        Assert.AreEqual(1.5, m.RefractiveIndexAt(400), Tol);
    }

    [TestMethod]
    public void OpticalMedium_Dispersion_Varies_With_Wavelength()
    {
        var bk7 = OpticalMaterialLibrary.CrownGlass; // Abbe = 64.17
        double nBlue = bk7.RefractiveIndexAt(450);
        double nRed = bk7.RefractiveIndexAt(650);
        Assert.IsTrue(nBlue > nRed, "Blue light should have higher n than red");
    }

    // ══════════════════════════════════════════════════════
    //  Snell's Law
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void Snell_Normal_Incidence_Gives_Zero_Refraction()
    {
        var angle = OpticsExtensions.RefractionAngle(1.0, 1.5, 0.0);
        Assert.IsTrue(angle.HasValue);
        Assert.AreEqual(0.0, angle.Value, Tol);
    }

    [TestMethod]
    public void Snell_Air_To_Glass_30Degrees()
    {
        double theta1 = Math.PI / 6; // 30°
        double n1 = 1.0, n2 = 1.5;
        var theta2 = OpticsExtensions.RefractionAngle(n1, n2, theta1);
        Assert.IsTrue(theta2.HasValue);
        // sin(30°)/1.5 = 0.333..., asin = ~19.47°
        double expected = Math.Asin(Math.Sin(theta1) * n1 / n2);
        Assert.AreEqual(expected, theta2.Value, Tol);
    }

    [TestMethod]
    public void Snell_TotalInternalReflection()
    {
        double n1 = 1.5, n2 = 1.0;
        double critical = OpticsExtensions.CriticalAngle(n1, n2);
        // Angle beyond critical
        var angle = OpticsExtensions.RefractionAngle(n1, n2, critical + 0.1);
        Assert.IsNull(angle);
    }

    [TestMethod]
    public void CriticalAngle_Glass_To_Air()
    {
        double n1 = 1.5, n2 = 1.0;
        double critical = OpticsExtensions.CriticalAngle(n1, n2);
        double expected = Math.Asin(n2 / n1); // ~41.8°
        Assert.AreEqual(expected, critical, Tol);
    }

    [TestMethod]
    public void IsTotalInternalReflection_Works()
    {
        Assert.IsTrue(OpticsExtensions.IsTotalInternalReflection(1.5, 1.0, Math.PI / 4));
        Assert.IsFalse(OpticsExtensions.IsTotalInternalReflection(1.0, 1.5, Math.PI / 4));
    }

    // ══════════════════════════════════════════════════════
    //  Fresnel Equations
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void Fresnel_NormalIncidence_MatchesFormula()
    {
        // At normal incidence R = ((n1-n2)/(n1+n2))²
        double n1 = 1.0, n2 = 1.5;
        double expected = Math.Pow((n1 - n2) / (n1 + n2), 2);
        double actual = OpticsExtensions.FresnelReflectance(n1, n2, 0.0);
        Assert.AreEqual(expected, actual, 1e-4);
    }

    [TestMethod]
    public void Fresnel_TIR_Returns_One()
    {
        double R = OpticsExtensions.FresnelReflectance(1.5, 1.0, Math.PI / 3);
        Assert.AreEqual(1.0, R, Tol);
    }

    [TestMethod]
    public void Schlick_NormalIncidence_MatchesFresnel()
    {
        double n1 = 1.0, n2 = 1.5;
        double cosTheta = 1.0; // normal incidence
        double schlick = OpticsExtensions.SchlickReflectance(n1, n2, cosTheta);
        double fresnel = OpticsExtensions.FresnelReflectance(n1, n2, 0.0);
        Assert.AreEqual(fresnel, schlick, 1e-4);
    }

    // ══════════════════════════════════════════════════════
    //  Beer–Lambert
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void BeerLambert_Zero_Absorption_Full_Transmission()
    {
        double T = OpticsExtensions.BeerLambertTransmittance(0.0, 100.0);
        Assert.AreEqual(1.0, T, Tol);
    }

    [TestMethod]
    public void BeerLambert_Exponential_Decay()
    {
        double alpha = 0.5; // 1/m
        double d = 2.0; // m
        double T = OpticsExtensions.BeerLambertTransmittance(alpha, d);
        Assert.AreEqual(Math.Exp(-1.0), T, Tol);
    }

    // ══════════════════════════════════════════════════════
    //  Reflect / Refract vectors
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void Reflect_At_45_Degrees()
    {
        var incident = new Vector(1, -1, 0).GetUnitVector();
        var normal = new Vector(0, 1, 0);
        var reflected = OpticsExtensions.Reflect(incident, normal);
        Assert.AreEqual(incident.x, reflected.x, Tol);
        Assert.AreEqual(-incident.y, reflected.y, Tol);
    }

    [TestMethod]
    public void Refract_Normal_Incidence_Same_Direction()
    {
        var incident = new Vector(0, -1, 0);
        var normal = new Vector(0, 1, 0);
        var refracted = OpticsExtensions.Refract(incident, normal, 1.0, 1.5);
        Assert.IsNotNull(refracted);
        // Expect straight through (same direction)
        Assert.AreEqual(0.0, refracted.Value.x, Tol);
        Assert.AreEqual(-1.0, refracted.Value.y, Tol);
    }

    [TestMethod]
    public void Refract_TIR_Returns_Null()
    {
        // Glass to air at steep angle
        var incident = new Vector(0.8, -0.6, 0).GetUnitVector();
        var normal = new Vector(0, 1, 0);
        var refracted = OpticsExtensions.Refract(incident, normal, 1.5, 1.0);
        Assert.IsNull(refracted);
    }

    // ══════════════════════════════════════════════════════
    //  Plane Mirror
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void PlaneMirror_Intersection()
    {
        var mirror = new PlaneMirror(new Vector(0, 0, 5), new Vector(0, 0, -1));
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var hit = mirror.Intersect(ray);
        Assert.IsNotNull(hit);
        Assert.AreEqual(5.0, hit.Value.Distance, Tol);
        Assert.AreEqual(5.0, hit.Value.Point.z, Tol);
    }

    [TestMethod]
    public void PlaneMirror_Parallel_Ray_No_Hit()
    {
        var mirror = new PlaneMirror(new Vector(0, 0, 5), new Vector(0, 0, 1));
        var ray = new Ray(new Vector(0, 0, 0), new Vector(1, 0, 0));
        var hit = mirror.Intersect(ray);
        Assert.IsNull(hit);
    }

    // ══════════════════════════════════════════════════════
    //  Spherical Mirror
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void SphericalMirror_ImageDistance_Concave()
    {
        var mirror = new SphericalMirror(new Vector(0, 0, 0), 20.0, MirrorType.Concave);
        // f = 10, object at 30 => 1/10 = 1/30 + 1/di => di = 15
        double di = mirror.ImageDistance(30);
        Assert.AreEqual(15.0, di, Tol);
    }

    [TestMethod]
    public void SphericalMirror_Magnification()
    {
        var mirror = new SphericalMirror(new Vector(0, 0, 0), 20.0, MirrorType.Concave);
        double m = mirror.Magnification(30);
        Assert.AreEqual(-0.5, m, Tol);
    }

    // ══════════════════════════════════════════════════════
    //  Thin Lens
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void ThinLens_ImageDistance_Converging()
    {
        var lens = new ThinLens(new Vector(0, 0, 0), new Vector(0, 0, 1), 10.0);
        // f=10, do=20 => di=20
        double di = lens.ImageDistance(20);
        Assert.AreEqual(20.0, di, Tol);
    }

    [TestMethod]
    public void ThinLens_Magnification()
    {
        var lens = new ThinLens(new Vector(0, 0, 0), new Vector(0, 0, 1), 10.0);
        double m = lens.Magnification(20);
        Assert.AreEqual(-1.0, m, Tol);
    }

    [TestMethod]
    public void ThinLens_Lensmaker_Equation()
    {
        double n = 1.5168;
        double r1 = 0.2, r2 = -0.2;
        double f = ThinLens.LensmakerFocalLength(n, r1, r2);
        double expected = 1.0 / ((n - 1.0) * (1.0 / r1 - 1.0 / r2));
        Assert.AreEqual(expected, f, Tol);
    }

    [TestMethod]
    public void ThinLens_Intersect_OnAxis()
    {
        var lens = new ThinLens(new Vector(0, 0, 5), new Vector(0, 0, 1), 10.0);
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var hit = lens.Intersect(ray);
        Assert.IsNotNull(hit);
        Assert.AreEqual(5.0, hit.Value.Distance, Tol);
    }

    [TestMethod]
    public void ThinLens_Aperture_Blocks_Ray()
    {
        var lens = new ThinLens(
            new Vector(0, 0, 5), new Vector(0, 0, 1), 10.0,
            apertureRadius: 0.5);
        // Ray aimed well outside the aperture
        var ray = new Ray(new Vector(10, 10, 0), new Vector(0, 0, 1));
        var hit = lens.Intersect(ray);
        Assert.IsNull(hit);
    }

    // ══════════════════════════════════════════════════════
    //  Prism
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void Prism_MinimumDeviation()
    {
        var prism = new Prism(Math.PI / 3, OpticalMaterialLibrary.CrownGlass);
        double dMin = prism.MinimumDeviation(589.3);
        // For BK7 n~1.5168, A=60°: δ_min = 2*asin(1.5168*sin(30°)) - 60°
        double expected = 2.0 * Math.Asin(1.5168 * Math.Sin(Math.PI / 6)) - Math.PI / 3;
        Assert.AreEqual(expected, dMin, 1e-4);
    }

    [TestMethod]
    public void Prism_Dispersion_Blue_Deviates_More_Than_Red()
    {
        var prism = new Prism(Math.PI / 3, OpticalMaterialLibrary.CrownGlass);
        double dBlue = prism.MinimumDeviation(450);
        double dRed = prism.MinimumDeviation(650);
        Assert.IsTrue(dBlue > dRed, "Blue should deviate more than red");
    }

    [TestMethod]
    public void Prism_DeviationAngle_Matches_MinDevAtSymmetry()
    {
        var prism = new Prism(Math.PI / 3, OpticalMaterialLibrary.CrownGlass);
        // At minimum deviation the incidence angle = (δ_min + A) / 2
        double dMin = prism.MinimumDeviation(589.3);
        double thetaSymmetric = (dMin + Math.PI / 3) / 2.0;
        double? dev = prism.DeviationAngle(thetaSymmetric, 589.3);
        Assert.IsNotNull(dev);
        Assert.AreEqual(dMin, dev.Value, 1e-3);
    }

    // ══════════════════════════════════════════════════════
    //  Circular Aperture
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void CircularAperture_Passes_Ray_Inside()
    {
        var ap = new CircularAperture(new Vector(0, 0, 5), new Vector(0, 0, 1), 2.0);
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var hit = ap.Intersect(ray);
        Assert.IsNotNull(hit);
    }

    [TestMethod]
    public void CircularAperture_Blocks_Ray_Outside()
    {
        var ap = new CircularAperture(new Vector(0, 0, 5), new Vector(0, 0, 1), 0.5);
        var ray = new Ray(new Vector(3, 3, 0), new Vector(0, 0, 1));
        var hit = ap.Intersect(ray);
        Assert.IsNull(hit);
    }

    // ══════════════════════════════════════════════════════
    //  Rectangular Aperture
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void RectangularAperture_Passes_CentreRay()
    {
        var ap = new RectangularAperture(
            new Vector(0, 0, 5), new Vector(0, 0, 1), new Vector(1, 0, 0), 4.0, 2.0);
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var hit = ap.Intersect(ray);
        Assert.IsNotNull(hit);
    }

    [TestMethod]
    public void RectangularAperture_Blocks_Outside()
    {
        var ap = new RectangularAperture(
            new Vector(0, 0, 5), new Vector(0, 0, 1), new Vector(1, 0, 0), 1.0, 1.0);
        var ray = new Ray(new Vector(5, 5, 0), new Vector(0, 0, 1));
        var hit = ap.Intersect(ray);
        Assert.IsNull(hit);
    }

    // ══════════════════════════════════════════════════════
    //  Image Sensor
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void ImageSensor_Records_Hits()
    {
        var sensor = new ImageSensor(
            new Vector(0, 0, 10), new Vector(0, 0, 1), new Vector(1, 0, 0),
            4.0, 4.0, 64, 64);
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var hit = sensor.Intersect(ray);
        Assert.IsNotNull(hit);
        Assert.AreEqual(1, sensor.HitCount);
    }

    [TestMethod]
    public void ImageSensor_Clear_Resets()
    {
        var sensor = new ImageSensor(
            new Vector(0, 0, 10), new Vector(0, 0, 1), new Vector(1, 0, 0),
            4.0, 4.0, 64, 64);
        sensor.Intersect(new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1)));
        sensor.Clear();
        Assert.AreEqual(0, sensor.HitCount);
    }

    // ══════════════════════════════════════════════════════
    //  Ray Tracer — mirror scene
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void RayTracer_Mirror_Reflects_Ray()
    {
        var scene = new OpticalScene();
        scene.Add(new PlaneMirror(new Vector(0, 0, 10), new Vector(0, 0, -1)));

        var tracer = new RayTracer(scene);
        var result = tracer.Trace(new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1)));

        Assert.IsTrue(result.Segments.Count >= 1);
        Assert.AreEqual(10.0, result.Segments[0].End.z, Tol);
    }

    [TestMethod]
    public void RayTracer_TwoMirrors_Bounces()
    {
        var scene = new OpticalScene();
        scene.Add(new PlaneMirror(new Vector(0, 0, 10), new Vector(0, 0, -1)));
        scene.Add(new PlaneMirror(new Vector(0, 0, 0), new Vector(0, 0, 1)));

        var tracer = new RayTracer(scene) { MaxBounces = 5 };
        var result = tracer.Trace(new Ray(new Vector(0, 0, 1), new Vector(0, 0, 1)));

        // Should bounce several times between the two mirrors
        Assert.IsTrue(result.Segments.Count >= 3);
    }

    [TestMethod]
    public void RayTracer_IntensityCutoff_TerminatesRay()
    {
        var scene = new OpticalScene();
        scene.Add(new PlaneMirror(new Vector(0, 0, 10), new Vector(0, 0, -1), reflectivity: 0.1));
        scene.Add(new PlaneMirror(new Vector(0, 0, 0), new Vector(0, 0, 1), reflectivity: 0.1));

        var tracer = new RayTracer(scene) { MaxBounces = 100, IntensityCutoff = 0.01 };
        var result = tracer.Trace(new Ray(new Vector(0, 0, 1), new Vector(0, 0, 1)));

        // Low reflectivity should cause early termination
        Assert.IsTrue(result.Terminated);
        Assert.IsTrue(result.Segments.Count < 10);
    }

    [TestMethod]
    public void RayTracer_Lens_And_Sensor()
    {
        var scene = new OpticalScene();
        var lens = new ThinLens(new Vector(0, 0, 5), new Vector(0, 0, 1), 10.0);
        var sensor = new ImageSensor(
            new Vector(0, 0, 15), new Vector(0, 0, 1), new Vector(1, 0, 0),
            4.0, 4.0, 32, 32);
        scene.Add(lens);
        scene.Add(sensor);

        var tracer = new RayTracer(scene);
        var result = tracer.Trace(new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1)));

        Assert.IsTrue(result.Segments.Count >= 2, "Ray should pass through lens then hit sensor");
        Assert.AreEqual(1, sensor.HitCount);
    }

    // ══════════════════════════════════════════════════════
    //  Material Library & Factory
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void OpticalMaterialLibrary_AllMedia_Have_Positive_N()
    {
        var media = new[]
        {
            OpticalMaterialLibrary.Vacuum, OpticalMaterialLibrary.Air,
            OpticalMaterialLibrary.Water, OpticalMaterialLibrary.CrownGlass,
            OpticalMaterialLibrary.FlintGlass, OpticalMaterialLibrary.Diamond,
            OpticalMaterialLibrary.Acrylic, OpticalMaterialLibrary.Polycarbonate,
        };
        foreach (var m in media)
            Assert.IsTrue(m.RefractiveIndex > 0, $"{m.Name} has invalid n");
    }

    [TestMethod]
    public void OpticalLibrary_Get_CaseInsensitive()
    {
        var m = OpticalLibrary.Get("crown glass");
        Assert.IsTrue(m.RefractiveIndex > 1.5);
    }

    [TestMethod]
    public void Materials_Optical_Factory()
    {
        var m = Materials.Optical("Diamond");
        Assert.AreEqual(2.417, m.RefractiveIndex, Tol);
    }

    [TestMethod]
    [ExpectedException(typeof(ArgumentException))]
    public void OpticalLibrary_Get_Throws_On_Unknown()
    {
        OpticalLibrary.Get("Unobtainium");
    }

    // ══════════════════════════════════════════════════════
    //  Game Engine — Raycast Extensions
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void Raycast_AABB_Hit()
    {
        var aabb = new AABB(new Vector(-1, -1, -1), new Vector(1, 1, 1));
        var ray = new Ray(new Vector(0, 0, -5), new Vector(0, 0, 1));
        var t = ray.IntersectAABB(aabb);
        Assert.IsNotNull(t);
        Assert.AreEqual(4.0, t.Value, Tol);
    }

    [TestMethod]
    public void Raycast_AABB_Miss()
    {
        var aabb = new AABB(new Vector(-1, -1, -1), new Vector(1, 1, 1));
        var ray = new Ray(new Vector(5, 5, -5), new Vector(0, 0, 1));
        var t = ray.IntersectAABB(aabb);
        Assert.IsNull(t);
    }

    [TestMethod]
    public void Raycast_Sphere_Hit()
    {
        var sphere = new BoundingSphere(new Vector(0, 0, 5), 1.0);
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var t = ray.IntersectSphere(sphere);
        Assert.IsNotNull(t);
        Assert.AreEqual(4.0, t.Value, Tol);
    }

    [TestMethod]
    public void Raycast_Sphere_Miss()
    {
        var sphere = new BoundingSphere(new Vector(0, 0, 5), 1.0);
        var ray = new Ray(new Vector(5, 0, 0), new Vector(0, 0, 1));
        var t = ray.IntersectSphere(sphere);
        Assert.IsNull(t);
    }

    [TestMethod]
    public void Raycast_Sphere_Detailed_Returns_Normal()
    {
        var sphere = new BoundingSphere(new Vector(0, 0, 5), 1.0);
        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var result = ray.IntersectSphereDetailed(sphere);
        Assert.IsNotNull(result);
        Assert.AreEqual(0.0, result.Value.normal.x, Tol);
        Assert.AreEqual(0.0, result.Value.normal.y, Tol);
        Assert.AreEqual(-1.0, result.Value.normal.z, Tol);
    }

    // ══════════════════════════════════════════════════════
    //  Optical Scene
    // ══════════════════════════════════════════════════════

    [TestMethod]
    public void OpticalScene_ClosestHit_Returns_Nearest()
    {
        var scene = new OpticalScene();
        scene.Add(new PlaneMirror(new Vector(0, 0, 10), new Vector(0, 0, -1)));
        scene.Add(new PlaneMirror(new Vector(0, 0, 5), new Vector(0, 0, -1)));

        var ray = new Ray(new Vector(0, 0, 0), new Vector(0, 0, 1));
        var hit = scene.ClosestHit(ray);
        Assert.IsNotNull(hit);
        Assert.AreEqual(5.0, hit.Value.Distance, Tol);
    }
}
