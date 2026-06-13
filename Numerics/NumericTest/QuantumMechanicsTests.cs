using System;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Constants;
using CSharpNumerics.Physics.Quantum;

namespace NumericsTests;

[TestClass]
public class QuantumMechanicsTests
{
    #region Time-independent Schrödinger equation

    [TestMethod]
    public void InfiniteSquareWell_NumericEnergies_MatchAnalytic()
    {
        // V = 0 inside [0, 1], infinite walls. Natural units (ħ = 1, m = 1).
        Func<double, double> freeParticle = _ => 0.0;
        var states = freeParticle.SolveStationaryStates(xMin: 0.0, xMax: 1.0, points: 150, states: 3);

        for (int n = 1; n <= 3; n++)
        {
            double analytic = n.InfiniteSquareWellEnergy(width: 1.0);   // n²π²/2
            double numeric = states.Energies[n - 1];
            Assert.AreEqual(analytic, numeric, analytic * 0.01,
                $"Level {n}: numeric {numeric:F4} should match analytic {analytic:F4} within 1%.");
        }

        // Ground-state probability is symmetric about the centre → ⟨x⟩ ≈ 0.5.
        var ground = states.WaveFunctions[0].ToComplexWaveFunction();
        Assert.AreEqual(0.5, ground.ExpectationPosition(states.Grid, states.GridSpacing), 0.02);
    }

    [TestMethod]
    public void HarmonicOscillator_NumericEnergies_MatchAnalytic()
    {
        // V = ½ x² (m = 1, ω = 1, ħ = 1) → Eₙ = n + ½.
        Func<double, double> harmonic = x => 0.5 * x * x;
        var states = harmonic.SolveStationaryStates(xMin: -8.0, xMax: 8.0, points: 200, states: 3);

        for (int n = 0; n < 3; n++)
        {
            double analytic = n.HarmonicOscillatorEnergy(angularFrequency: 1.0);   // n + 0.5
            Assert.AreEqual(analytic, states.Energies[n], 0.02,
                $"Level {n}: numeric {states.Energies[n]:F4} should match analytic {analytic:F4}.");
        }
    }

    #endregion

    #region Wavefunction observables

    [TestMethod]
    public void WaveFunction_Normalization_GivesUnitProbability()
    {
        int n = 400;
        double xMin = -10, xMax = 10, dx = (xMax - xMin) / (n - 1);
        var grid = new double[n];
        var psi = new ComplexNumber[n];
        for (int i = 0; i < n; i++)
        {
            double x = xMin + i * dx;
            grid[i] = x;
            double g = Math.Exp(-(x - 2.0) * (x - 2.0) / (2.0 * 1.0 * 1.0));  // Gaussian, centre 2, σ = 1
            psi[i] = new ComplexNumber(g, 0.0);
        }

        var normalized = psi.Normalize(dx);
        Assert.AreEqual(1.0, normalized.NormSquared(dx), 1e-9);

        var density = normalized.ProbabilityDensity();
        foreach (var d in density) Assert.IsTrue(d >= 0.0, "Probability density must be non-negative.");

        // ⟨x⟩ ≈ centre; Δx ≈ σ/√2 for a Gaussian |ψ|².
        Assert.AreEqual(2.0, normalized.ExpectationPosition(grid, dx), 0.01);
        Assert.AreEqual(1.0 / Math.Sqrt(2.0), normalized.PositionUncertainty(grid, dx), 0.02);
    }

    [TestMethod]
    public void ExpectationMomentum_MatchesPlaneWaveModulation()
    {
        int n = 600;
        double xMin = -12, xMax = 12, dx = (xMax - xMin) / (n - 1);
        const double k = 2.0;

        var realGaussian = new ComplexNumber[n];
        var modulated = new ComplexNumber[n];
        for (int i = 0; i < n; i++)
        {
            double x = xMin + i * dx;
            double g = Math.Exp(-(x * x) / (2.0 * 2.0 * 2.0));   // σ = 2
            realGaussian[i] = new ComplexNumber(g, 0.0);
            modulated[i] = new ComplexNumber(Math.Cos(k * x) * g, Math.Sin(k * x) * g);  // e^{ikx} g
        }

        // Real wavefunction carries no net momentum.
        Assert.AreEqual(0.0, realGaussian.ExpectationMomentum(dx), 1e-6);

        // A plane-wave-modulated packet has ⟨p⟩ ≈ ħk = k (ħ = 1).
        Assert.AreEqual(k, modulated.ExpectationMomentum(dx), 0.05);
    }

    #endregion

    #region Tunnelling

    [TestMethod]
    public void RectangularBarrier_Transmission_BehavesPhysically()
    {
        const double e = 1.0, v0 = 2.0;   // E < V₀ → tunnelling regime

        double tThin = e.RectangularBarrierTransmission(v0, barrierWidth: 0.5);
        double tMid = e.RectangularBarrierTransmission(v0, barrierWidth: 1.0);
        double tThick = e.RectangularBarrierTransmission(v0, barrierWidth: 2.0);

        foreach (var t in new[] { tThin, tMid, tThick })
            Assert.IsTrue(t > 0.0 && t < 1.0, $"Transmission {t} must be strictly between 0 and 1.");

        // Wider barrier → exponentially less tunnelling.
        Assert.IsTrue(tThin > tMid && tMid > tThick, "Transmission must decrease with barrier width.");

        // Vanishing barrier → full transmission.
        Assert.AreEqual(1.0, e.RectangularBarrierTransmission(v0, barrierWidth: 1e-6), 1e-3);

        // Reflection complements transmission.
        Assert.AreEqual(1.0, tMid + e.RectangularBarrierReflection(v0, 1.0), 1e-12);
    }

    #endregion

    #region de Broglie

    [TestMethod]
    public void DeBroglieWavelength_MatchesPlanckOverMomentum()
    {
        double momentum = 5.0e-25;   // kg·m/s
        Assert.AreEqual(PhysicsConstants.PlancksConstant / momentum, momentum.DeBroglieWavelength(), 1e-40);

        // λ from kinetic energy should equal λ from the corresponding momentum p = √(2mE).
        double mass = PhysicsConstants.ElectronMass;
        double energy = 1.0e-18;
        double p = Math.Sqrt(2.0 * mass * energy);
        Assert.AreEqual(p.DeBroglieWavelength(), energy.DeBroglieWavelengthFromEnergy(mass), 1e-18);
    }

    #endregion

    #region Time evolution

    [TestMethod]
    public void TimeEvolution_PreservesNorm_AndStationaryDensity()
    {
        Func<double, double> well = _ => 0.0;
        var states = well.SolveStationaryStates(xMin: 0.0, xMax: 1.0, points: 120, states: 4);
        double dx = states.GridSpacing;

        // Superposition of the first two states (Σ cₙ² = 1).
        double inv = 1.0 / Math.Sqrt(2.0);
        var coefficients = new[] { inv, inv, 0.0, 0.0 };

        var psi0 = states.Evolve(coefficients, time: 0.0);
        var psiT = states.Evolve(coefficients, time: 5.0);

        Assert.AreEqual(1.0, psi0.NormSquared(dx), 1e-9, "Norm at t=0 should be 1.");
        Assert.AreEqual(1.0, psiT.NormSquared(dx), 1e-9, "Norm must be conserved under time evolution.");

        // A single stationary state only accrues a global phase → its density is time-independent.
        var stationary0 = states.Evolve(new[] { 1.0, 0.0, 0.0, 0.0 }, time: 0.0).ProbabilityDensity();
        var stationaryT = states.Evolve(new[] { 1.0, 0.0, 0.0, 0.0 }, time: 3.7).ProbabilityDensity();
        for (int i = 0; i < stationary0.Length; i++)
            Assert.AreEqual(stationary0[i], stationaryT[i], 1e-9);
    }

    #endregion
}
