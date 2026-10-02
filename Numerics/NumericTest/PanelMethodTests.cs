using CSharpNumerics.Physics.FluidDynamics.Aerodynamics;

namespace NumericTest;

/// <summary>
/// Characterisation tests for the 2-D panel method.
/// </summary>
/// <remarks>
/// Written to pin the solver's behaviour down before its hand-rolled Gaussian elimination is
/// replaced with the shared LU decomposition.
///
/// These assert invariants that hold whatever the solver computes — symmetry, invariance of Cp
/// under freestream scaling, net zero source strength, finiteness, sensitivity to incidence —
/// rather than physical correctness. That is deliberate: the one test here that did check against
/// a closed-form solution does not pass, and is left ignored with the question written down. So
/// read this file as a change detector, not as validation that the physics is right.
/// </remarks>
[TestClass]
public class PanelMethodTests
{
    /// <summary>
    /// Builds a closed circular contour, traversed clockwise so the panel normals point
    /// outward in the same sense the solver assumes for an airfoil.
    /// </summary>
    private static (double[] x, double[] y) Cylinder(int panels, double radius = 1.0)
    {
        var x = new double[panels + 1];
        var y = new double[panels + 1];

        for (var i = 0; i <= panels; i++)
        {
            var theta = -2.0 * Math.PI * i / panels;
            x[i] = radius * Math.Cos(theta);
            y[i] = radius * Math.Sin(theta);
        }

        return (x, y);
    }

    /// <summary>A symmetric diamond, the simplest closed body with sharp corners.</summary>
    private static (double[] x, double[] y) Diamond()
    {
        var x = new double[] { 1.0, 0.5, 0.0, 0.5, 1.0 };
        var y = new double[] { 0.0, -0.1, 0.0, 0.1, 0.0 };
        return (x, y);
    }

    /// <summary>
    /// Open question, not a passing assertion. Potential flow around a cylinder has
    /// Cp = 1 − 4 sin²θ, ranging over [−3, 1]. This solver instead returns Cp ≈ 1 on every
    /// panel (measured range [0.9999, 1.0000] at 200 panels, zero incidence), which means the
    /// tangential velocity comes out ≈ 0 everywhere — not a potential flow around a body.
    /// </summary>
    /// <remarks>
    /// Two readings are possible and this test does not settle which: the solver may be wrong,
    /// or a seamed cylinder may simply be outside its domain, since it imposes a Kutta condition
    /// at a trailing edge that a cylinder does not have. Only the clockwise traversal was
    /// measured. Left ignored rather than deleted so the question is not lost — investigating it
    /// is separate work from the v4.3 migration.
    /// </remarks>
    [TestMethod]
    [Ignore("Unresolved: the solver returns Cp ~= 1 everywhere for a cylinder. See the summary above.")]
    public void Cylinder_MatchesAnalyticPressureDistribution()
    {
        // Potential flow around a cylinder: Cp = 1 - 4 sin²θ, independent of radius.
        var (x, y) = Cylinder(200);

        var result = PanelMethod.Solve(x, y, alpha: 0.0);

        for (var i = 0; i < result.Cp.Length; i++)
        {
            var theta = Math.Atan2(result.Ym[i], result.Xm[i]);
            var expected = 1.0 - 4.0 * Math.Sin(theta) * Math.Sin(theta);

            Assert.AreEqual(expected, result.Cp[i], 0.05,
                $"Panel {i} at theta = {theta:F3}");
        }
    }

    [TestMethod]
    public void Cylinder_IsSymmetricAtZeroIncidence()
    {
        var (x, y) = Cylinder(120);

        var result = PanelMethod.Solve(x, y, alpha: 0.0);

        // A body symmetric about y = 0 at zero incidence must carry a symmetric
        // pressure distribution.
        for (var i = 0; i < result.Cp.Length; i++)
        {
            var mirror = NearestPanel(result, result.Xm[i], -result.Ym[i]);

            Assert.AreEqual(result.Cp[i], result.Cp[mirror], 1e-6,
                $"Panel {i} and its mirror {mirror}");
        }
    }

    [TestMethod]
    public void FreestreamScaling_LeavesPressureCoefficientUnchanged()
    {
        var (x, y) = Cylinder(80);

        var unit = PanelMethod.Solve(x, y, alpha: 0.0, freestream: 1.0);
        var scaled = PanelMethod.Solve(x, y, alpha: 0.0, freestream: 7.5);

        // Cp is normalised by dynamic pressure, so it cannot depend on the freestream speed.
        for (var i = 0; i < unit.Cp.Length; i++)
        {
            Assert.AreEqual(unit.Cp[i], scaled.Cp[i], 1e-9, $"Panel {i}");
        }
    }

    [TestMethod]
    public void SourceStrengths_SumToApproximatelyZeroForAClosedBody()
    {
        var (x, y) = Cylinder(160);

        var result = PanelMethod.Solve(x, y, alpha: 0.0);

        // A closed body in incompressible flow neither creates nor destroys mass.
        var net = 0.0;
        for (var i = 0; i < result.Sigma.Length; i++)
        {
            net += result.Sigma[i];
        }

        Assert.AreEqual(0.0, net, 1e-6);
    }

    [TestMethod]
    public void Diamond_ProducesFiniteResultsOnASharpBody()
    {
        var (x, y) = Diamond();

        var result = PanelMethod.Solve(x, y, alpha: 0.05);

        Assert.AreEqual(4, result.Cp.Length);

        foreach (var cp in result.Cp)
        {
            Assert.IsFalse(double.IsNaN(cp), "Cp must not be NaN");
            Assert.IsFalse(double.IsInfinity(cp), "Cp must not be infinite");
        }

        foreach (var vt in result.Vt)
        {
            Assert.IsFalse(double.IsNaN(vt), "Vt must not be NaN");
        }
    }

    [TestMethod]
    public void IncidenceChangesThePressureDistribution()
    {
        var (x, y) = Cylinder(120);

        var straight = PanelMethod.Solve(x, y, alpha: 0.0);
        var inclined = PanelMethod.Solve(x, y, alpha: 0.2);

        var changed = false;
        for (var i = 0; i < straight.Cp.Length; i++)
        {
            if (Math.Abs(straight.Cp[i] - inclined.Cp[i]) > 1e-6)
            {
                changed = true;
                break;
            }
        }

        Assert.IsTrue(changed, "Angle of attack must affect the solution.");
    }

    [TestMethod]
    public void TooFewPoints_Throws()
    {
        Assert.ThrowsException<ArgumentException>(
            () => PanelMethod.Solve(new double[] { 0, 1, 0 }, new double[] { 0, 0, 0 }, 0.0));
    }

    [TestMethod]
    public void NullCoordinates_Throws()
    {
        Assert.ThrowsException<ArgumentNullException>(
            () => PanelMethod.Solve(null!, null!, 0.0));
    }

    private static int NearestPanel(PanelMethodResult result, double x, double y)
    {
        var best = 0;
        var bestDistance = double.MaxValue;

        for (var i = 0; i < result.Xm.Length; i++)
        {
            var dx = result.Xm[i] - x;
            var dy = result.Ym[i] - y;
            var distance = dx * dx + dy * dy;

            if (distance < bestDistance)
            {
                bestDistance = distance;
                best = i;
            }
        }

        return best;
    }
}
