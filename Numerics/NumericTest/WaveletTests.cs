using System;
using System.Linq;
using CSharpNumerics.Numerics.SignalProcessing.Wavelets;

namespace NumericsTests;

[TestClass]
public class WaveletTests
{
    private static double NextGaussian(Random rng, double std)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return std * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double Energy(double[] x) => x.Sum(v => v * v);

    #region Round-trip

    [TestMethod]
    public void DwtIdwt_RoundTrip_ReconstructsWithinTolerance()
    {
        var families = new[]
        {
            WaveletFamily.Haar,
            WaveletFamily.Daubechies4,
            WaveletFamily.Daubechies8,
            WaveletFamily.Symlet4
        };

        var rng = new Random(31);
        int n = 256;
        var signal = new double[n];
        for (int i = 0; i < n; i++) signal[i] = rng.NextDouble() * 10.0 - 5.0;

        foreach (var family in families)
        {
            var dwt = new DiscreteWaveletTransform(family);
            var idwt = new InverseWaveletTransform(family);

            var decomposition = dwt.Decompose(signal, levels: 3);
            var reconstructed = idwt.Reconstruct(decomposition);

            double maxError = 0.0;
            for (int i = 0; i < n; i++)
                maxError = Math.Max(maxError, Math.Abs(reconstructed[i] - signal[i]));

            // Reconstruction is exact; the residual is the floating-point accumulation floor
            // (~1e-11 for 8-tap filters via the transpose inverse over a 256-point, 3-level transform).
            Assert.IsTrue(maxError < 1e-10,
                $"{family.Name}: reconstruction error {maxError:E3} exceeds the machine-precision floor.");
        }
    }

    #endregion

    #region Haar step function

    [TestMethod]
    public void Haar_OnStepFunction_LocalisesDetailAtDiscontinuity()
    {
        // Step within an even-indexed pair (samples 4,5): the discontinuity is captured by a
        // single Haar detail coefficient, with zeros across the flat regions.
        var signal = new double[] { 0, 0, 0, 0, 0, 1, 1, 1 };

        var dwt = new DiscreteWaveletTransform(WaveletFamily.Haar);
        var (_, detail) = dwt.SingleLevel(signal);

        Assert.AreEqual(0.0, detail[0], 1e-12);
        Assert.AreEqual(0.0, detail[1], 1e-12);
        Assert.AreEqual(-1.0 / Math.Sqrt(2.0), detail[2], 1e-12);   // straddles the 0→1 step
        Assert.AreEqual(0.0, detail[3], 1e-12);
    }

    #endregion

    #region Energy localisation

    [TestMethod]
    public void Daubechies4_OnLowFrequencySine_ConcentratesEnergyInApproximation()
    {
        int n = 256;
        var signal = new double[n];
        for (int i = 0; i < n; i++) signal[i] = Math.Sin(2.0 * Math.PI * i / 64.0);  // low frequency

        var dwt = new DiscreteWaveletTransform(WaveletFamily.Daubechies4);
        var decomposition = dwt.Decompose(signal, levels: 4);

        double totalEnergy = Energy(signal);
        double approxEnergy = Energy(decomposition.Approximation);
        double detailEnergy = decomposition.Details.Sum(Energy);

        // Parseval: orthonormal transform preserves energy.
        Assert.AreEqual(totalEnergy, approxEnergy + detailEnergy, 1e-9);

        // A smooth low-frequency sine puts almost all energy into the coarse approximation.
        Assert.IsTrue(approxEnergy / totalEnergy > 0.9,
            $"Approximation holds only {approxEnergy / totalEnergy:P1} of the energy.");

        // The finest detail (highest frequency) should be nearly empty.
        Assert.IsTrue(Energy(decomposition.Details[0]) / totalEnergy < 0.01,
            "Finest detail band should carry negligible energy for a low-frequency sine.");
    }

    #endregion

    #region Denoising

    [TestMethod]
    public void WaveletDenoising_ReducesRmse_ByMoreThanHalf()
    {
        int n = 2048;
        var rng = new Random(7);
        var clean = new double[n];
        var noisy = new double[n];
        for (int i = 0; i < n; i++)
        {
            clean[i] = Math.Sin(2.0 * Math.PI * i / 256.0);   // low-frequency → sparse in wavelet domain
            noisy[i] = clean[i] + NextGaussian(rng, 0.4);
        }

        var denoised = WaveletDenoising.Denoise(
            noisy, WaveletFamily.Symlet4, levels: 6,
            ThresholdType.Soft, ThresholdRule.VisuShrink);

        double rmseNoisy = Rmse(noisy, clean);
        double rmseDenoised = Rmse(denoised, clean);

        Assert.IsTrue(rmseDenoised < 0.5 * rmseNoisy,
            $"Denoised RMSE ({rmseDenoised:F4}) should be less than half of noisy RMSE ({rmseNoisy:F4}).");
    }

    private static double Rmse(double[] a, double[] b)
    {
        double sse = 0.0;
        for (int i = 0; i < a.Length; i++) sse += (a[i] - b[i]) * (a[i] - b[i]);
        return Math.Sqrt(sse / a.Length);
    }

    #endregion

    #region MODWT

    [TestMethod]
    public void Modwt_IsShiftInvariant()
    {
        int n = 64, shift = 5, levels = 3;
        var rng = new Random(99);
        var signal = new double[n];
        for (int i = 0; i < n; i++) signal[i] = rng.NextDouble() * 2.0 - 1.0;

        // Circularly shift the input by `shift`.
        var shifted = new double[n];
        for (int t = 0; t < n; t++) shifted[t] = signal[((t - shift) % n + n) % n];

        var modwt = new MaximalOverlapDWT(WaveletFamily.Daubechies4);
        var baseDecomp = modwt.Forward(signal, levels);
        var shiftedDecomp = modwt.Forward(shifted, levels);

        // Each shifted-input band must equal the shifted original band.
        for (int level = 0; level < levels; level++)
            for (int t = 0; t < n; t++)
            {
                double expected = baseDecomp.Details[level][((t - shift) % n + n) % n];
                Assert.AreEqual(expected, shiftedDecomp.Details[level][t], 1e-9,
                    $"Level {level}, t={t}: MODWT is not shift-invariant.");
            }

        for (int t = 0; t < n; t++)
        {
            double expected = baseDecomp.Smooth[((t - shift) % n + n) % n];
            Assert.AreEqual(expected, shiftedDecomp.Smooth[t], 1e-9);
        }
    }

    [TestMethod]
    public void Modwt_RoundTrip_ReconstructsWithinTolerance()
    {
        int n = 128;
        var rng = new Random(3);
        var signal = new double[n];
        for (int i = 0; i < n; i++) signal[i] = rng.NextDouble() * 4.0 - 2.0;

        var modwt = new MaximalOverlapDWT(WaveletFamily.Symlet4);
        var decomposition = modwt.Forward(signal, levels: 4);
        var reconstructed = modwt.Inverse(decomposition);

        double maxError = 0.0;
        for (int i = 0; i < n; i++)
            maxError = Math.Max(maxError, Math.Abs(reconstructed[i] - signal[i]));

        Assert.IsTrue(maxError < 1e-9, $"MODWT round-trip error {maxError:E3} exceeds 1e-9.");
    }

    #endregion
}
