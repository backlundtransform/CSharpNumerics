using CSharpNumerics.Engines.Audio;
using CSharpNumerics.Engines.Common;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.IO;

namespace NumericsTests
{
    [TestClass]
    public class AudioEngineTests
    {
        #region 3.0 Common Infrastructure

        [TestMethod]
        public void SimulationClock_Accumulate_ReturnsCorrectSteps()
        {
            var clock = new SimulationClock { FixedDt = 1.0 / 60.0 };
            int steps = clock.Accumulate(1.0 / 30.0); // 2 fixed steps
            Assert.AreEqual(2, steps);
            Assert.AreEqual(2, clock.TickCount);
        }

        [TestMethod]
        public void SimulationClock_Alpha_PartialStep()
        {
            var clock = new SimulationClock { FixedDt = 0.01 };
            clock.Accumulate(0.015); // 1 step + 0.005 remainder
            Assert.AreEqual(1, clock.TickCount);
            Assert.AreEqual(0.5, clock.Alpha, 0.01);
        }

        [TestMethod]
        public void SimulationClock_Reset()
        {
            var clock = new SimulationClock();
            clock.Accumulate(1.0);
            clock.Reset();
            Assert.AreEqual(0, clock.TickCount);
            Assert.AreEqual(0.0, clock.ElapsedTime, 1e-12);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void SimulationClock_NegativeDt_Throws()
        {
            var clock = new SimulationClock();
            clock.Accumulate(-0.1);
        }

        [TestMethod]
        public void EventBus_PublishSubscribe()
        {
            var bus = new EventBus();
            string received = "";
            bus.Subscribe<string>(s => received = s);
            bus.Publish("hello");
            Assert.AreEqual("hello", received);
        }

        [TestMethod]
        public void EventBus_Unsubscribe()
        {
            var bus = new EventBus();
            int count = 0;
            Action<int> handler = x => count += x;
            bus.Subscribe(handler);
            bus.Publish(1);
            bus.Unsubscribe(handler);
            bus.Publish(1);
            Assert.AreEqual(1, count);
        }

        [TestMethod]
        public void EventBus_Clear()
        {
            var bus = new EventBus();
            int count = 0;
            bus.Subscribe<int>(x => count++);
            bus.Clear();
            bus.Publish(42);
            Assert.AreEqual(0, count);
        }

        [TestMethod]
        public void FieldSerializer_BinaryRoundTrip()
        {
            var data = new double[] { 1.1, 2.2, 3.3, 4.4 };
            string path = Path.Combine(Path.GetTempPath(), "test_field.bin");
            try
            {
                FieldSerializer.SaveBinary(data, path);
                var loaded = FieldSerializer.LoadBinary(path);
                Assert.AreEqual(data.Length, loaded.Length);
                for (int i = 0; i < data.Length; i++)
                    Assert.AreEqual(data[i], loaded[i], 1e-15);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void FieldSerializer_CsvRoundTrip()
        {
            var data = new double[] { 1.5, -2.7, 0.0, 100.01 };
            string path = Path.Combine(Path.GetTempPath(), "test_field.csv");
            try
            {
                FieldSerializer.SaveCsv(data, path);
                var loaded = FieldSerializer.LoadCsv(path);
                Assert.AreEqual(data.Length, loaded.Length);
                for (int i = 0; i < data.Length; i++)
                    Assert.AreEqual(data[i], loaded[i], 1e-10);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [TestMethod]
        public void FieldSerializer_JsonRoundTrip()
        {
            var data = new double[] { 3.14, -1.0, 0.0, 42.5 };
            string json = FieldSerializer.ToJson(data);
            var loaded = FieldSerializer.FromJson(json);
            Assert.AreEqual(data.Length, loaded.Length);
            for (int i = 0; i < data.Length; i++)
                Assert.AreEqual(data[i], loaded[i], 1e-15);
        }

        [TestMethod]
        public void PlatformAdapter_DefaultIsStandalone()
        {
            var adapter = PlatformAdapter.Default;
            Assert.AreEqual(RuntimePlatform.Standalone, adapter.Platform);
        }

        #endregion

        #region 3.1 AudioBuffer

        [TestMethod]
        public void AudioBuffer_Constructor_SetsProperties()
        {
            var buf = new AudioBuffer(44100, 1, 1.0);
            Assert.AreEqual(44100, buf.SampleRate);
            Assert.AreEqual(1, buf.Channels);
            Assert.AreEqual(44100, buf.FrameCount);
            Assert.AreEqual(1.0, buf.Duration, 0.001);
        }

        [TestMethod]
        public void AudioBuffer_Stereo()
        {
            var buf = new AudioBuffer(44100, 2, 0.5);
            Assert.AreEqual(2, buf.Channels);
            Assert.AreEqual(22050, buf.FrameCount);
            buf[0, 0] = 0.5;
            buf[0, 1] = -0.5;
            Assert.AreEqual(0.5, buf[0, 0]);
            Assert.AreEqual(-0.5, buf[0, 1]);
        }

        [TestMethod]
        public void AudioBuffer_ToMono_AveragesChannels()
        {
            var stereo = new AudioBuffer(new double[] { 0.8, 0.2, -0.4, 0.6 }, 44100, 2);
            var mono = stereo.ToMono();
            Assert.AreEqual(1, mono.Channels);
            Assert.AreEqual(2, mono.FrameCount);
            Assert.AreEqual(0.5, mono.Samples[0], 1e-10);
            Assert.AreEqual(0.1, mono.Samples[1], 1e-10);
        }

        [TestMethod]
        public void AudioBuffer_MixIn_Adds()
        {
            var a = new AudioBuffer(new double[] { 0.5, 0.3 }, 44100, 1);
            var b = new AudioBuffer(new double[] { 0.1, 0.2 }, 44100, 1);
            a.MixIn(b);
            Assert.AreEqual(0.6, a.Samples[0], 1e-10);
            Assert.AreEqual(0.5, a.Samples[1], 1e-10);
        }

        [TestMethod]
        public void AudioBuffer_Normalize_PeakAt1()
        {
            var buf = new AudioBuffer(new double[] { 0.5, -2.0, 1.0 }, 44100, 1);
            buf.Normalize();
            Assert.AreEqual(0.25, buf.Samples[0], 1e-10);
            Assert.AreEqual(-1.0, buf.Samples[1], 1e-10);
            Assert.AreEqual(0.5, buf.Samples[2], 1e-10);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void AudioBuffer_InvalidSampleRate_Throws()
        {
            new AudioBuffer(0, 1, 1.0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void AudioBuffer_InvalidChannels_Throws()
        {
            new AudioBuffer(44100, 3, 1.0);
        }

        #endregion

        #region 3.1 SignalGenerator

        [TestMethod]
        public void SignalGenerator_Sine_CorrectFrequency()
        {
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 1.0, 0.1, 44100);
            Assert.AreEqual(44100, buf.SampleRate);
            Assert.AreEqual(4410, buf.FrameCount);
            // At t=0 phase=0, sin(0)=0
            Assert.AreEqual(0.0, buf.Samples[0], 1e-10);
        }

        [TestMethod]
        public void SignalGenerator_Square_Values()
        {
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Square, 1, 1.0, 1.0, 100);
            // First half = +1, second half = -1
            Assert.AreEqual(1.0, buf.Samples[0], 1e-10);
            Assert.AreEqual(-1.0, buf.Samples[50], 1e-10);
        }

        [TestMethod]
        public void SignalGenerator_Sawtooth_RampsUp()
        {
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sawtooth, 1, 1.0, 1.0, 100);
            // At t=0: 2*0-1 = -1
            Assert.AreEqual(-1.0, buf.Samples[0], 1e-10);
            // At t~0.5: 2*0.5-1 = 0
            Assert.AreEqual(0.0, buf.Samples[50], 0.05);
        }

        [TestMethod]
        public void SignalGenerator_WhiteNoise_HasVariance()
        {
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.WhiteNoise, 0, 1.0, 0.1, 44100);
            double mean = 0;
            foreach (var s in buf.Samples) mean += s;
            mean /= buf.FrameCount;
            // Mean should be close to 0 for uniform noise
            Assert.AreEqual(0.0, mean, 0.1);
        }

        [TestMethod]
        public void SignalGenerator_Triangle_Symmetric()
        {
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Triangle, 1, 1.0, 1.0, 1000);
            // At t=0: frac=0 → 4*0-1 = -1
            Assert.AreEqual(-1.0, buf.Samples[0], 1e-10);
            // At t=0.25: frac=0.25 → 4*0.25-1 = 0
            Assert.AreEqual(0.0, buf.Samples[250], 0.01);
            // At t≈0.5: frac≈0.5 → peak = 1
            Assert.AreEqual(1.0, buf.Samples[499], 0.02);
        }

        #endregion

        #region 3.1 AudioOscillator

        [TestMethod]
        public void AudioOscillator_GeneratesBuffer()
        {
            var osc = new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 0.8);
            var buf = osc.GenerateBuffer(0.05, 44100);
            Assert.AreEqual(2205, buf.FrameCount);
            Assert.IsTrue(Math.Abs(buf.Samples[0]) < 0.01); // sin starts near 0
        }

        [TestMethod]
        public void AudioOscillator_ContinuousPhase()
        {
            var osc = new AudioOscillator(SignalGenerator.Waveform.Sine, 100, 1.0);
            double s1 = osc.NextSample(1000);
            double s2 = osc.NextSample(1000);
            // Phase advances smoothly
            Assert.IsTrue(s2 > s1); // sine increasing from 0
        }

        [TestMethod]
        public void AudioOscillator_Reset_ClearsPhase()
        {
            var osc = new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 1.0);
            osc.NextSample(44100);
            osc.NextSample(44100);
            osc.Reset();
            Assert.AreEqual(0, osc.Phase, 1e-12);
        }

        #endregion

        #region 3.1 Envelope (ADSR)

        [TestMethod]
        public void Envelope_AttackPhase()
        {
            var env = new Envelope(0.1, 0.1, 0.5, 0.2);
            // Halfway through attack: ~0.5
            Assert.AreEqual(0.5, env.Evaluate(0.05), 0.01);
            // End of attack: 1.0
            Assert.AreEqual(1.0, env.Evaluate(0.1), 0.01);
        }

        [TestMethod]
        public void Envelope_DecayPhase()
        {
            var env = new Envelope(0.1, 0.1, 0.5, 0.2);
            // End of decay: sustain level
            Assert.AreEqual(0.5, env.Evaluate(0.2), 0.01);
        }

        [TestMethod]
        public void Envelope_SustainPhase()
        {
            var env = new Envelope(0.1, 0.1, 0.7, 0.2);
            // During sustain
            Assert.AreEqual(0.7, env.Evaluate(0.5), 0.01);
            Assert.AreEqual(0.7, env.Evaluate(1.0), 0.01);
        }

        [TestMethod]
        public void Envelope_ReleasePhase()
        {
            var env = new Envelope(0.01, 0.01, 0.8, 0.1);
            double noteOff = 0.5;
            // Just after note off: starts from sustain
            Assert.AreEqual(0.8, env.Evaluate(noteOff, noteOff), 0.02);
            // End of release: 0
            Assert.AreEqual(0.0, env.Evaluate(noteOff + 0.1, noteOff), 0.01);
        }

        [TestMethod]
        public void Envelope_Apply_ModulatesBuffer()
        {
            var env = new Envelope(0.01, 0.01, 1.0, 0.01);
            // Fill with constant 1.0
            var buf = new AudioBuffer(new double[] { 1.0, 1.0, 1.0, 1.0 }, 100, 1);
            env.Apply(buf, double.MaxValue);
            // First sample at t=0: attack start = 0
            Assert.AreEqual(0.0, buf.Samples[0], 0.01);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Envelope_NegativeAttack_Throws()
        {
            new Envelope(-0.1, 0.1, 0.5, 0.1);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Envelope_InvalidSustain_Throws()
        {
            new Envelope(0.1, 0.1, 1.5, 0.1);
        }

        #endregion

        #region 3.1 Synthesizer

        [TestMethod]
        public void Synthesizer_SingleVoice_MatchesOscillator()
        {
            var synth = new Synthesizer { SampleRate = 8000 };
            synth.AddVoice(new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 1.0));
            var buf = synth.Render(0.01);
            Assert.AreEqual(80, buf.FrameCount);
            Assert.IsTrue(Math.Abs(buf.Samples[0]) < 0.01); // sin starts near 0
        }

        [TestMethod]
        public void Synthesizer_TwoVoices_AreAdditive()
        {
            var synth = new Synthesizer { SampleRate = 8000 };
            synth.AddVoice(new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 0.5));
            synth.AddVoice(new AudioOscillator(SignalGenerator.Waveform.Sine, 880, 0.3));
            Assert.AreEqual(2, synth.VoiceCount);
            var buf = synth.Render(0.01);
            Assert.AreEqual(80, buf.FrameCount);
        }

        [TestMethod]
        public void Synthesizer_WithEnvelope_ShapesOutput()
        {
            var synth = new Synthesizer { SampleRate = 1000 };
            synth.AddVoice(new AudioOscillator(SignalGenerator.Waveform.Sine, 100, 1.0));
            synth.Envelope = new Envelope(0.1, 0.1, 0.5, 0.2);
            var buf = synth.Render(0.5, 0.4);
            // At start (attack phase), amplitude should be small
            Assert.IsTrue(Math.Abs(buf.Samples[0]) < 0.01);
        }

        [TestMethod]
        public void Synthesizer_Reset()
        {
            var synth = new Synthesizer { SampleRate = 8000 };
            var osc = new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 1.0);
            synth.AddVoice(osc);
            synth.Render(0.01);
            synth.Reset();
            Assert.AreEqual(0, osc.Phase, 1e-12);
        }

        #endregion

        #region 3.2 AudioFilter

        [TestMethod]
        public void AudioFilter_LowPass_AttenuatesHighFreqs()
        {
            // Generate a signal with two components: 100 Hz + 5000 Hz
            int sr = 44100;
            double dur = 0.1;
            var low = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 100, 0.5, dur, sr);
            var high = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 5000, 0.5, dur, sr);
            low.MixIn(high);

            var filtered = AudioFilter.Apply(low, AudioFilter.FilterType.LowPass, 500);

            // After low-pass at 500 Hz, the 5000 Hz component should be gone
            // Check energy: filtered should have ~half the energy of original
            double energyFiltered = 0, energyOriginal = 0;
            for (int i = 0; i < filtered.FrameCount; i++)
            {
                energyFiltered += filtered.Samples[i] * filtered.Samples[i];
            }
            // Regenerate original for comparison
            var orig = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 100, 0.5, dur, sr);
            var origHigh = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 5000, 0.5, dur, sr);
            orig.MixIn(origHigh);
            for (int i = 0; i < orig.FrameCount; i++)
                energyOriginal += orig.Samples[i] * orig.Samples[i];

            Assert.IsTrue(energyFiltered < energyOriginal * 0.75);
        }

        [TestMethod]
        public void AudioFilter_HighPass_AttenuatesLowFreqs()
        {
            int sr = 44100;
            double dur = 0.1;
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 100, 1.0, dur, sr);
            var filtered = AudioFilter.Apply(buf, AudioFilter.FilterType.HighPass, 0, 500);

            double energy = 0;
            for (int i = 0; i < filtered.FrameCount; i++)
                energy += filtered.Samples[i] * filtered.Samples[i];

            // 100 Hz should be almost completely removed by a 500 Hz high-pass
            Assert.IsTrue(energy < 0.1);
        }

        [TestMethod]
        public void AudioFilter_BandPass_PassesBand()
        {
            int sr = 44100;
            double dur = 0.1;
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 1000, 1.0, dur, sr);

            var filtered = AudioFilter.Apply(buf, AudioFilter.FilterType.BandPass, 500, 2000);

            double energyFiltered = 0, energyOriginal = 0;
            for (int i = 0; i < filtered.FrameCount; i++)
                energyFiltered += filtered.Samples[i] * filtered.Samples[i];
            for (int i = 0; i < buf.FrameCount; i++)
                energyOriginal += buf.Samples[i] * buf.Samples[i];

            // 1000 Hz is within [500, 2000], so energy should be preserved
            Assert.IsTrue(energyFiltered > energyOriginal * 0.8);
        }

        #endregion

        #region 3.2 Reverb

        [TestMethod]
        public void Reverb_Process_ReturnsBuffer()
        {
            var reverb = new Reverb(0.5, 0.5, 0.3);
            var input = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 0.5, 0.1, 44100);
            var output = reverb.Process(input);
            Assert.AreEqual(input.FrameCount, output.FrameCount);
            Assert.AreEqual(input.SampleRate, output.SampleRate);
        }

        [TestMethod]
        public void Reverb_DrySignal_WhenWetMixZero()
        {
            var reverb = new Reverb(0.5, 0.5, 0.0);
            var input = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 0.5, 0.05, 44100);
            var output = reverb.Process(input);
            // With wetMix=0, output should equal input
            for (int i = 0; i < input.FrameCount; i++)
                Assert.AreEqual(input.Samples[i], output.Samples[i], 1e-10);
        }

        #endregion

        #region 3.2 Delay

        [TestMethod]
        public void Delay_Process_ReturnsBuffer()
        {
            var delay = new Delay(0.1, 0.5, 0.5);
            var input = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 0.5, 0.5, 44100);
            var output = delay.Process(input);
            Assert.AreEqual(input.FrameCount, output.FrameCount);
        }

        [TestMethod]
        public void Delay_DrySignal_WhenWetMixZero()
        {
            var delay = new Delay(0.1, 0.5, 0.0);
            var input = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 0.5, 0.1, 44100);
            var output = delay.Process(input);
            for (int i = 0; i < input.FrameCount; i++)
                Assert.AreEqual(input.Samples[i], output.Samples[i], 1e-10);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Delay_InvalidDelayTime_Throws()
        {
            new Delay(0);
        }

        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Delay_InvalidFeedback_Throws()
        {
            new Delay(0.1, 1.0);
        }

        #endregion

        #region 3.2 Compressor

        [TestMethod]
        public void Compressor_ReducesLoudSignals()
        {
            var comp = new Compressor(0.3, 4.0, 0.001, 0.01);
            var input = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 0.9, 0.1, 44100);
            var output = comp.Process(input);

            // Peak of output should be less than peak of input
            double peakIn = 0, peakOut = 0;
            for (int i = 0; i < input.FrameCount; i++)
            {
                peakIn = Math.Max(peakIn, Math.Abs(input.Samples[i]));
                peakOut = Math.Max(peakOut, Math.Abs(output.Samples[i]));
            }

            Assert.IsTrue(peakOut < peakIn);
        }

        [TestMethod]
        public void Compressor_PassesQuietSignals()
        {
            var comp = new Compressor(0.8, 4.0, 0.001, 0.01);
            var input = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 0.1, 0.1, 44100);
            var output = comp.Process(input);

            // Signal below threshold should pass through mostly unchanged
            double diffEnergy = 0;
            for (int i = 0; i < input.FrameCount; i++)
            {
                double diff = output.Samples[i] - input.Samples[i];
                diffEnergy += diff * diff;
            }

            Assert.IsTrue(diffEnergy < 0.01);
        }

        #endregion

        #region 3.2 SpatialAudio

        [TestMethod]
        public void SpatialAudio_Pan_Center_EqualChannels()
        {
            var mono = new AudioBuffer(new double[] { 1.0, 1.0, 1.0 }, 44100, 1);
            var stereo = SpatialAudio.Pan(mono, 0.0);
            Assert.AreEqual(2, stereo.Channels);
            // Center pan: both channels should be equal (constant-power)
            Assert.AreEqual(stereo[0, 0], stereo[0, 1], 0.01);
        }

        [TestMethod]
        public void SpatialAudio_Pan_Left_LeftLouder()
        {
            var mono = new AudioBuffer(new double[] { 1.0, 1.0 }, 44100, 1);
            var stereo = SpatialAudio.Pan(mono, -1.0);
            Assert.IsTrue(Math.Abs(stereo[0, 0]) > Math.Abs(stereo[0, 1]));
        }

        [TestMethod]
        public void SpatialAudio_Pan_Right_RightLouder()
        {
            var mono = new AudioBuffer(new double[] { 1.0, 1.0 }, 44100, 1);
            var stereo = SpatialAudio.Pan(mono, 1.0);
            Assert.IsTrue(Math.Abs(stereo[0, 1]) > Math.Abs(stereo[0, 0]));
        }

        [TestMethod]
        public void SpatialAudio_Distance_AttenuatesWithDistance()
        {
            var buf = new AudioBuffer(new double[] { 1.0, 1.0, 1.0 }, 44100, 1);
            var near = SpatialAudio.AttenuateByDistance(buf, 1.0);
            var far = SpatialAudio.AttenuateByDistance(buf, 10.0);
            Assert.IsTrue(Math.Abs(far.Samples[0]) < Math.Abs(near.Samples[0]));
        }

        [TestMethod]
        public void SpatialAudio_Spatialize_CombinesPanAndDistance()
        {
            var mono = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 1.0, 0.01, 8000);
            var result = SpatialAudio.Spatialize(mono, 5, 3, 0, 0);
            Assert.AreEqual(2, result.Channels);
            // Source is to the right → right channel louder
            Assert.IsTrue(Math.Abs(result[0, 1]) >= Math.Abs(result[0, 0]) - 0.01);
        }

        #endregion

        #region 3.3 SpectrumAnalyzer

        [TestMethod]
        public void SpectrumAnalyzer_DetectsPeakFrequency()
        {
            int sr = 8192;
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 1.0, 0.5, sr);
            var analyzer = new SpectrumAnalyzer(2048, SpectrumAnalyzer.WindowType.Hann);
            var spectrum = analyzer.Analyze(buf.Samples, sr);

            // Find peak
            double maxMag = 0;
            double peakFreq = 0;
            foreach (var s in spectrum)
            {
                if (s.Value > maxMag)
                {
                    maxMag = s.Value;
                    peakFreq = s.Index;
                }
            }

            // Peak should be near 440 Hz (within FFT resolution = sr/fftSize = 4 Hz)
            Assert.AreEqual(440, peakFreq, 10);
        }

        [TestMethod]
        public void SpectrumAnalyzer_AnalyzeBuffer_ReturnsAveragedSpectrum()
        {
            int sr = 8192;
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 1000, 1.0, 0.5, sr);
            var analyzer = new SpectrumAnalyzer(1024, SpectrumAnalyzer.WindowType.Hamming);
            var spectrum = analyzer.AnalyzeBuffer(buf);
            Assert.IsTrue(spectrum.Count > 0);
        }

        [TestMethod]
        public void SpectrumAnalyzer_WindowFunctions()
        {
            // Hann window: endpoints = 0, middle = 1
            Assert.AreEqual(0.0, SpectrumAnalyzer.ApplyWindow(SpectrumAnalyzer.WindowType.Hann, 0, 256), 1e-10);
            Assert.AreEqual(1.0, SpectrumAnalyzer.ApplyWindow(SpectrumAnalyzer.WindowType.Hann, 127, 255), 0.01);

            // Hamming: endpoints ≈ 0.08
            Assert.AreEqual(0.08, SpectrumAnalyzer.ApplyWindow(SpectrumAnalyzer.WindowType.Hamming, 0, 256), 0.01);

            // Rectangular: always 1
            Assert.AreEqual(1.0, SpectrumAnalyzer.ApplyWindow(SpectrumAnalyzer.WindowType.Rectangular, 0, 256), 1e-10);
        }

        #endregion

        #region 3.3 PitchDetector

        [TestMethod]
        public void PitchDetector_Autocorrelation_DetectsA440()
        {
            int sr = 44100;
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 1.0, 0.1, sr);
            var detector = new PitchDetector(4096);
            double freq = detector.Detect(buf, PitchDetector.Method.Autocorrelation);
            Assert.AreEqual(440, freq, 5);
        }

        [TestMethod]
        public void PitchDetector_HPS_DetectsA440()
        {
            int sr = 44100;
            // HPS needs harmonics — use a square wave (has odd harmonics)
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Square, 440, 1.0, 0.2, sr);
            var detector = new PitchDetector(8192);
            double freq = detector.Detect(buf, PitchDetector.Method.HarmonicProductSpectrum);
            Assert.AreEqual(440, freq, 20);
        }

        [TestMethod]
        public void PitchDetector_DetectsLowPitch()
        {
            int sr = 44100;
            var buf = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 100, 1.0, 0.2, sr);
            var detector = new PitchDetector(8192) { MinFrequency = 50, MaxFrequency = 500 };
            double freq = detector.Detect(buf);
            Assert.AreEqual(100, freq, 10);
        }

        #endregion

        #region 3.3 BeatDetector

        [TestMethod]
        public void BeatDetector_DetectsOnsetsInPulsedSignal()
        {
            // Create a signal with 4 "beats" = short bursts at regular intervals
            int sr = 8000;
            double totalDuration = 2.0;
            var buf = new AudioBuffer(sr, 1, totalDuration);
            double interval = 0.5; // 120 BPM
            int burstSamples = 200;

            for (int beat = 0; beat < 4; beat++)
            {
                int start = (int)(beat * interval * sr);
                for (int i = 0; i < burstSamples && start + i < buf.FrameCount; i++)
                {
                    double t = (double)i / sr;
                    buf.Samples[start + i] = 0.9 * Math.Sin(2 * Math.PI * 440 * t);
                }
            }

            var detector = new BeatDetector(512, 256) { Threshold = 1.0 };
            var onsets = detector.Detect(buf);

            // Should detect at least 2 onsets
            Assert.IsTrue(onsets.Count >= 2, $"Expected at least 2 onsets, got {onsets.Count}");
        }

        [TestMethod]
        public void BeatDetector_EstimateTempo()
        {
            // 4 equidistant onsets at 120 BPM
            var onsets = new List<double> { 0.0, 0.5, 1.0, 1.5 };
            double bpm = BeatDetector.EstimateTempoFromOnsets(onsets);
            Assert.AreEqual(120, bpm, 1);
        }

        [TestMethod]
        public void BeatDetector_SilentSignal_NoOnsets()
        {
            var buf = new AudioBuffer(8000, 1, 1.0);
            var detector = new BeatDetector();
            var onsets = detector.Detect(buf);
            Assert.AreEqual(0, onsets.Count);
        }

        #endregion
    }
}
