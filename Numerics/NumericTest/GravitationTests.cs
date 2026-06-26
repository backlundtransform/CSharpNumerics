using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Gravitation;

namespace NumericsTests;

[TestClass]
public class GravitationTests
{
    private const double G = PhysicsConstants.GravitationalConstant;
    private const double SolarMass = PhysicsConstants.SolarMass;
    private const double EarthMass = PhysicsConstants.EarthMass;
    private const double MoonMass = PhysicsConstants.MoonMass;
    private const double AU = PhysicsConstants.AstronomicalUnit;

    #region N-body

    [TestMethod]
    public void PairwiseAcceleration_MatchesNewtonInverseSquare()
    {
        // Body at origin, mass M at distance d along +x → a = GM/d² toward +x
        double d = 1.0e9, M = 1.0e24;
        var a = NBodyGravity.PairwiseAcceleration(new Vector(0, 0, 0), new Vector(d, 0, 0), M);

        Assert.AreEqual(G * M / (d * d), a.x, G * M / (d * d) * 1e-12);
        Assert.AreEqual(0.0, a.y, 1e-20);
    }

    [TestMethod]
    public void Accelerations_ConserveMomentum()
    {
        // Σ mᵢ·aᵢ = 0 for internal gravity (Newton's third law)
        var pos = new[] { new Vector(0, 0, 0), new Vector(1e9, 0, 0), new Vector(0, 2e9, 0) };
        var mass = new[] { 1e24, 2e24, 3e24 };

        var acc = NBodyGravity.Accelerations(pos, mass);

        double px = 0, py = 0, pz = 0;
        for (int i = 0; i < acc.Length; i++)
        {
            px += mass[i] * acc[i].x;
            py += mass[i] * acc[i].y;
            pz += mass[i] * acc[i].z;
        }
        // Reference scale for the tolerance
        double scale = mass[0] * NBodyGravity.Acceleration(pos, mass, 0).GetMagnitude();
        Assert.AreEqual(0.0, px, scale * 1e-10);
        Assert.AreEqual(0.0, py, scale * 1e-10);
        Assert.AreEqual(0.0, pz, scale * 1e-10);
    }

    [TestMethod]
    public void PotentialEnergy_TwoBody_MatchesAnalytic()
    {
        var pos = new[] { new Vector(0, 0, 0), new Vector(3e8, 0, 0) };
        var mass = new[] { EarthMass, MoonMass };
        double expected = -G * EarthMass * MoonMass / 3e8;
        Assert.AreEqual(expected, NBodyGravity.PotentialEnergy(pos, mass), Math.Abs(expected) * 1e-12);
    }

    [TestMethod]
    public void CenterOfMass_TwoEqualMasses_IsMidpoint()
    {
        var pos = new[] { new Vector(0, 0, 0), new Vector(10, 0, 0) };
        var mass = new[] { 5.0, 5.0 };
        var com = NBodyGravity.CenterOfMass(pos, mass);
        Assert.AreEqual(5.0, com.x, 1e-12);
    }

    #endregion

    #region Celestial mechanics

    [TestMethod]
    public void VisViva_CircularOrbit_EqualsCircularSpeed()
    {
        // r = a → v = √(μ/r)
        double mu = G * EarthMass;
        double r = PhysicsConstants.EarthRadius + 4.0e5;
        Assert.AreEqual(Math.Sqrt(mu / r), CelestialMechanics.VisVivaSpeed(mu, r, r), 1e-6);
    }

    [TestMethod]
    public void ReducedMass_IsSymmetricAndBelowSmaller()
    {
        double m1 = 5.0, m2 = 3.0;
        double mu = CelestialMechanics.ReducedMass(m1, m2);
        Assert.AreEqual(CelestialMechanics.ReducedMass(m2, m1), mu, 1e-12);
        Assert.IsTrue(mu < Math.Min(m1, m2));      // μ = 15/8 = 1.875
        Assert.AreEqual(1.875, mu, 1e-12);
    }

    [TestMethod]
    public void HillSphere_Earth_IsAbout1Point5MillionKm()
    {
        double rH = CelestialMechanics.HillSphereRadius(AU, 0.0, EarthMass, SolarMass);
        Assert.AreEqual(1.5e9, rH, 0.1e9); // ≈ 1.5×10⁶ km
    }

    [TestMethod]
    public void SphereOfInfluence_Earth_IsAbout924000Km()
    {
        double soi = CelestialMechanics.SphereOfInfluence(AU, EarthMass, SolarMass);
        Assert.AreEqual(9.2e8, soi, 0.4e8); // ≈ 924 000 km
    }

    [TestMethod]
    public void RocheLimit_Fluid_ExceedsRigid()
    {
        // Moon-like satellite around Earth
        double rEarth = PhysicsConstants.EarthRadius;
        double rhoEarth = 5514, rhoMoon = 3344;

        double fluid = CelestialMechanics.RocheLimitFluid(rEarth, rhoEarth, rhoMoon);
        double rigid = CelestialMechanics.RocheLimitRigid(rEarth, rhoEarth, rhoMoon);

        Assert.IsTrue(fluid > rigid);
        Assert.AreEqual(1.84e7, fluid, 0.1e7); // ≈ 18 400 km
    }

    [TestMethod]
    public void TidalAcceleration_FollowsInverseCube()
    {
        // Doubling the distance reduces the tidal acceleration by 2³ = 8
        double a1 = CelestialMechanics.TidalAcceleration(MoonMass, 3.844e8, PhysicsConstants.EarthRadius);
        double a2 = CelestialMechanics.TidalAcceleration(MoonMass, 2.0 * 3.844e8, PhysicsConstants.EarthRadius);
        Assert.AreEqual(a1 / 8.0, a2, a1 * 1e-12);
        Assert.IsTrue(a1 > 0);
    }

    #endregion

    #region Lagrange points

    [TestMethod]
    public void LagrangePoints_EarthSun_L1_IsAbout1Point5MillionKmFromEarth()
    {
        var (l1, l2, _, _, _) = LagrangePoints.All(SolarMass, EarthMass, AU);

        // L1 sits between the Sun (origin) and Earth (AU, 0, 0), ~1.5e9 m inside Earth
        double l1FromEarth = AU - l1.x;
        Assert.AreEqual(1.5e9, l1FromEarth, 0.1e9);
        Assert.IsTrue(l1.x < AU && l1.x > 0);

        // L2 lies just beyond Earth at a similar distance
        double l2FromEarth = l2.x - AU;
        Assert.AreEqual(1.5e9, l2FromEarth, 0.1e9);
    }

    [TestMethod]
    public void LagrangePoints_TriangularPointsFormEquilateralTriangles()
    {
        double sep = AU;
        var (_, _, _, l4, l5) = LagrangePoints.All(SolarMass, EarthMass, sep);

        // L4 equidistant (= separation) from both primaries
        double fromSun = l4.GetMagnitude();
        double dx = l4.x - sep, dy = l4.y;
        double fromEarth = Math.Sqrt(dx * dx + dy * dy);

        Assert.AreEqual(sep, fromSun, sep * 1e-12);
        Assert.AreEqual(sep, fromEarth, sep * 1e-12);

        // L5 is the mirror image across the axis
        Assert.AreEqual(l4.x, l5.x, 1e-6);
        Assert.AreEqual(-l4.y, l5.y, 1e-6);
    }

    #endregion
}
