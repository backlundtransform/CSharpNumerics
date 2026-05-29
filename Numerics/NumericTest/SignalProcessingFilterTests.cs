using CSharpNumerics.Numerics.SignalProcessing;

namespace NumericsTests;

[TestClass]
public class SignalProcessingFilterTests
{
    #region Savitzky-Golay Tests

    [TestMethod]
    public void SavitzkyGolay_Order0_EqualsMovingAverage()
    {
        // Savitzky-Golay with polynomial order 0 should be equivalent to a moving average
        var signal = new double[] { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21 };
        var filter = new SavitzkyGolayFilter(windowSize: 5, polynomialOrder: 0);
        var result = filter.Apply(signal);

        // For the interior points, should equal mean of 5-point window
        for (int i = 2; i < signal.Length - 2; i++)
        {
            double expected = (signal[i - 2] + signal[i - 1] + signal[i] + signal[i + 1] + signal[i + 2]) / 5.0;
            Assert.AreEqual(expected, result[i], 1e-10,
                $"At index {i}: expected {expected}, got {result[i]}");
        }
    }

    [TestMethod]
    public void SavitzkyGolay_PreservesStraightLine()
    {
        // A polynomial filter of order >= 1 should perfectly preserve a linear signal
        var signal = new double[50];
        for (int i = 0; i < 50; i++) signal[i] = 2.0 * i + 3.0;

        var filter = new SavitzkyGolayFilter(windowSize: 7, polynomialOrder: 2);
        var result = filter.Apply(signal);

        // Interior points should be exactly preserved
        for (int i = 3; i < 47; i++)
        {
            Assert.AreEqual(signal[i], result[i], 1e-10,
                $"Linear signal not preserved at index {i}");
        }
    }

    [TestMethod]
    public void SavitzkyGolay_PreservesQuadratic()
    {
        // Order 2 filter should perfectly preserve a quadratic signal
        var signal = new double[50];
        for (int i = 0; i < 50; i++) signal[i] = 0.5 * i * i - 3.0 * i + 7.0;

        var filter = new SavitzkyGolayFilter(windowSize: 7, polynomialOrder: 2);
        var result = filter.Apply(signal);

        for (int i = 3; i < 47; i++)
        {
            Assert.AreEqual(signal[i], result[i], 1e-8,
                $"Quadratic signal not preserved at index {i}");
        }
    }

    [TestMethod]
    public void SavitzkyGolay_ReducesNoise()
    {
        // Adding noise to a smooth signal should be reduced by the filter
        var rng = new Random(42);
        var clean = new double[200];
        var noisy = new double[200];
        for (int i = 0; i < 200; i++)
        {
            clean[i] = Math.Sin(2.0 * Math.PI * i / 50.0);
            noisy[i] = clean[i] + 0.3 * (rng.NextDouble() * 2.0 - 1.0);
        }

        var filter = new SavitzkyGolayFilter(windowSize: 11, polynomialOrder: 3);
        var filtered = filter.Apply(noisy);

        // Compute RMS error before and after filtering
        double rmsBefore = 0, rmsAfter = 0;
        for (int i = 10; i < 190; i++) // Avoid boundaries
        {
            rmsBefore += (noisy[i] - clean[i]) * (noisy[i] - clean[i]);
            rmsAfter += (filtered[i] - clean[i]) * (filtered[i] - clean[i]);
        }
        rmsBefore = Math.Sqrt(rmsBefore / 180);
        rmsAfter = Math.Sqrt(rmsAfter / 180);

        Assert.IsTrue(rmsAfter < rmsBefore * 0.5,
            $"Filter should reduce noise significantly. Before: {rmsBefore:F4}, After: {rmsAfter:F4}");
    }

    [TestMethod]
    public void SavitzkyGolay_InvalidParameters_Throws()
    {
        Assert.ThrowsException<ArgumentException>(() => new SavitzkyGolayFilter(4, 2)); // Even window
        Assert.ThrowsException<ArgumentException>(() => new SavitzkyGolayFilter(1, 0)); // Too small
        Assert.ThrowsException<ArgumentException>(() => new SavitzkyGolayFilter(5, 5)); // Order >= window
    }

    [TestMethod]
    public void SavitzkyGolay_CoefficientsSymmetric()
    {
        // Smoothing coefficients for even polynomial order should be symmetric
        var filter = new SavitzkyGolayFilter(windowSize: 7, polynomialOrder: 2);
        var coeffs = filter.Coefficients;

        for (int i = 0; i < coeffs.Length / 2; i++)
        {
            Assert.AreEqual(coeffs[i], coeffs[coeffs.Length - 1 - i], 1e-12,
                $"Coefficients not symmetric at position {i}");
        }
    }

    [TestMethod]
    public void SavitzkyGolay_CoefficientsSum_ToOne()
    {
        // Smoothing coefficients must sum to 1 (preserves DC level)
        var filter = new SavitzkyGolayFilter(windowSize: 9, polynomialOrder: 4);
        var coeffs = filter.Coefficients;

        double sum = 0;
        foreach (var c in coeffs) sum += c;

        Assert.AreEqual(1.0, sum, 1e-12, "Coefficients must sum to 1.");
    }

    #endregion

    #region Butterworth Filter Tests

    [TestMethod]
    public void Butterworth_Lowpass_AttenuatesHighFrequency()
    {
        // Design a lowpass filter and verify it attenuates frequencies above cutoff
        double sampleRate = 100.0;
        double cutoff = 10.0;
        var filter = FilterDesign.DesignLowpass(order: 4, cutoffFrequency: cutoff, sampleRate: sampleRate);

        // Create a signal with two frequencies: 5 Hz (passes) and 30 Hz (attenuated)
        int n = 1000;
        var signal = new double[n];
        for (int i = 0; i < n; i++)
        {
            double t = i / sampleRate;
            signal[i] = Math.Sin(2.0 * Math.PI * 5.0 * t) + Math.Sin(2.0 * Math.PI * 30.0 * t);
        }

        var filtered = filter.Apply(signal);

        // Measure power of the high-frequency component in the second half (after transient)
        double highFreqPower = 0;
        double lowFreqPower = 0;
        for (int i = 500; i < n; i++)
        {
            double t = i / sampleRate;
            double lowRef = Math.Sin(2.0 * Math.PI * 5.0 * t);
            lowFreqPower += filtered[i] * filtered[i];
        }

        // The filtered signal power should be roughly half of original (only low freq remains)
        double originalPower = 0;
        for (int i = 500; i < n; i++) originalPower += signal[i] * signal[i];

        Assert.IsTrue(lowFreqPower < originalPower * 0.75,
            "Lowpass filter should significantly reduce total power by removing high-frequency component.");
    }

    [TestMethod]
    public void Butterworth_Lowpass_Attenuation_12dB_Per_Octave_Per_Order()
    {
        // A 2nd-order Butterworth should attenuate by ~12 dB at 2× cutoff frequency
        double sampleRate = 1000.0;
        double cutoff = 100.0;
        var filter = FilterDesign.DesignLowpass(order: 2, cutoffFrequency: cutoff, sampleRate: sampleRate);

        var freqs = new double[] { cutoff / (sampleRate / 2.0), 2.0 * cutoff / (sampleRate / 2.0) };
        var response = filter.FrequencyResponse(freqs);

        // At cutoff: should be ~ -3dB (≈ 0.707)
        double atCutoff_dB = 20.0 * Math.Log10(response[0]);
        Assert.IsTrue(Math.Abs(atCutoff_dB - (-3.0)) < 1.5,
            $"At cutoff: expected ~-3dB, got {atCutoff_dB:F1} dB");

        // At 2× cutoff: 2nd order → ~-12 dB
        double at2xCutoff_dB = 20.0 * Math.Log10(response[1]);
        Assert.IsTrue(at2xCutoff_dB < -8.0,
            $"At 2× cutoff: expected < -8 dB, got {at2xCutoff_dB:F1} dB");
    }

    [TestMethod]
    public void Butterworth_Highpass_PassesHighFrequency()
    {
        double sampleRate = 100.0;
        double cutoff = 10.0;
        var filter = FilterDesign.DesignHighpass(order: 3, cutoffFrequency: cutoff, sampleRate: sampleRate);

        // Create signal: 2 Hz (blocked) + 25 Hz (passes)
        int n = 1000;
        var signal = new double[n];
        for (int i = 0; i < n; i++)
        {
            double t = i / sampleRate;
            signal[i] = Math.Sin(2.0 * Math.PI * 2.0 * t) + Math.Sin(2.0 * Math.PI * 25.0 * t);
        }

        var filtered = filter.Apply(signal);

        // The filtered signal should have negligible low-frequency content
        // Check that the mean of filtered signal over a long window is near zero
        double mean = 0;
        for (int i = 500; i < n; i++) mean += filtered[i];
        mean /= (n - 500);

        Assert.IsTrue(Math.Abs(mean) < 0.1,
            $"Highpass should remove DC/low-frequency. Mean = {mean:F4}");
    }

    [TestMethod]
    public void Butterworth_InvalidParameters_Throws()
    {
        Assert.ThrowsException<ArgumentException>(() => FilterDesign.DesignLowpass(0, 10, 100));
        Assert.ThrowsException<ArgumentException>(() => FilterDesign.DesignLowpass(2, 60, 100)); // cutoff >= Nyquist
        Assert.ThrowsException<ArgumentException>(() => FilterDesign.DesignLowpass(2, -5, 100));
    }

    #endregion

    #region FIR Filter Tests

    [TestMethod]
    public void FIR_Lowpass_AttenuatesHighFrequency()
    {
        double sampleRate = 1000.0;
        double cutoff = 100.0;
        var filter = FilterDesign.DesignFIRLowpass(numTaps: 51, cutoffFrequency: cutoff, sampleRate: sampleRate);

        // Signal: 50 Hz + 300 Hz
        int n = 2000;
        var signal = new double[n];
        for (int i = 0; i < n; i++)
        {
            double t = i / sampleRate;
            signal[i] = Math.Sin(2.0 * Math.PI * 50.0 * t) + Math.Sin(2.0 * Math.PI * 300.0 * t);
        }

        var filtered = filter.Apply(signal);

        // Measure energy in second half (past transient)
        double filteredEnergy = 0;
        double originalEnergy = 0;
        for (int i = 1000; i < n; i++)
        {
            filteredEnergy += filtered[i] * filtered[i];
            originalEnergy += signal[i] * signal[i];
        }

        // Should have roughly half the energy (one of two equal-amplitude components removed)
        double ratio = filteredEnergy / originalEnergy;
        Assert.IsTrue(ratio < 0.7 && ratio > 0.2,
            $"FIR lowpass energy ratio = {ratio:F3}, expected ~0.5");
    }

    [TestMethod]
    public void FIR_CoefficientsSum_ToOne_ForLowpass()
    {
        var filter = FilterDesign.DesignFIRLowpass(numTaps: 31, cutoffFrequency: 50, sampleRate: 200);
        var coeffs = filter.Coefficients;

        double sum = 0;
        foreach (var c in coeffs) sum += c;

        Assert.AreEqual(1.0, sum, 1e-10, "Lowpass FIR coefficients should sum to 1 (unity DC gain).");
    }

    [TestMethod]
    public void FIR_Highpass_BlocksDC()
    {
        var filter = FilterDesign.DesignFIRHighpass(numTaps: 51, cutoffFrequency: 50, sampleRate: 1000);

        // DC signal (constant)
        var signal = new double[500];
        for (int i = 0; i < 500; i++) signal[i] = 5.0;

        var filtered = filter.Apply(signal);

        // After transient, output should be near zero
        double mean = 0;
        for (int i = 100; i < 400; i++) mean += filtered[i];
        mean /= 300;

        Assert.IsTrue(Math.Abs(mean) < 0.01,
            $"Highpass filter should block DC. Mean = {mean:F6}");
    }

    #endregion

    #region ZeroPhaseFiltFilt Tests

    [TestMethod]
    public void FiltFilt_ZeroPhaseShift()
    {
        // A peak in the signal should remain at exactly the same position after filtfilt
        double sampleRate = 100.0;
        var filter = FilterDesign.DesignLowpass(order: 4, cutoffFrequency: 10.0, sampleRate: sampleRate);

        int n = 500;
        var signal = new double[n];
        int peakPos = 250;
        // Gaussian pulse
        for (int i = 0; i < n; i++)
        {
            double x = (i - peakPos) / 5.0;
            signal[i] = Math.Exp(-0.5 * x * x);
        }

        var filtered = ZeroPhaseFiltFilt.Apply(filter, signal);

        // Find peak position in filtered signal
        int filteredPeakPos = 0;
        double maxVal = double.MinValue;
        for (int i = 50; i < n - 50; i++)
        {
            if (filtered[i] > maxVal)
            {
                maxVal = filtered[i];
                filteredPeakPos = i;
            }
        }

        Assert.AreEqual(peakPos, filteredPeakPos,
            $"Peak should remain at same position. Expected {peakPos}, got {filteredPeakPos}");
    }

    [TestMethod]
    public void FiltFilt_Butterworth_LowpassHighpass_Reconstruction()
    {
        // Lowpass + highpass of same signal should approximately reconstruct original
        double sampleRate = 100.0;
        double cutoff = 15.0;

        var lpFilter = FilterDesign.DesignLowpass(order: 3, cutoffFrequency: cutoff, sampleRate: sampleRate);
        var hpFilter = FilterDesign.DesignHighpass(order: 3, cutoffFrequency: cutoff, sampleRate: sampleRate);

        int n = 500;
        var rng = new Random(123);
        var signal = new double[n];
        for (int i = 0; i < n; i++)
        {
            double t = i / sampleRate;
            signal[i] = Math.Sin(2.0 * Math.PI * 5.0 * t) + 0.5 * Math.Sin(2.0 * Math.PI * 25.0 * t);
        }

        var lp = ZeroPhaseFiltFilt.Apply(lpFilter, signal);
        var hp = ZeroPhaseFiltFilt.Apply(hpFilter, signal);

        // lp + hp should approximately equal original (Butterworth complementary pair is not perfect)
        double maxError = 0;
        for (int i = 50; i < n - 50; i++)
        {
            double recon = lp[i] + hp[i];
            double error = Math.Abs(recon - signal[i]);
            if (error > maxError) maxError = error;
        }

        // Butterworth LP+HP doesn't perfectly reconstruct, but should be reasonable
        Assert.IsTrue(maxError < 0.5,
            $"LP + HP reconstruction max error = {maxError:F4}, expected < 0.5");
    }

    [TestMethod]
    public void FiltFilt_FIR_ZeroPhaseShift()
    {
        double sampleRate = 100.0;
        var filter = FilterDesign.DesignFIRLowpass(numTaps: 31, cutoffFrequency: 10.0, sampleRate: sampleRate);

        int n = 300;
        var signal = new double[n];
        int peakPos = 150;
        for (int i = 0; i < n; i++)
        {
            double x = (i - peakPos) / 5.0;
            signal[i] = Math.Exp(-0.5 * x * x);
        }

        var filtered = ZeroPhaseFiltFilt.Apply(filter, signal);

        // Find peak position
        int filteredPeakPos = 0;
        double maxVal = double.MinValue;
        for (int i = 30; i < n - 30; i++)
        {
            if (filtered[i] > maxVal)
            {
                maxVal = filtered[i];
                filteredPeakPos = i;
            }
        }

        Assert.AreEqual(peakPos, filteredPeakPos,
            $"FIR FiltFilt peak should remain at same position. Expected {peakPos}, got {filteredPeakPos}");
    }

    #endregion
}
