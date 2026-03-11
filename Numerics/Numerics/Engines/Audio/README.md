## 🔊 Audio Engine

Digital audio synthesis, effects processing, and analysis — built on the Numerics FFT and Physics Waves foundations.

**Namespace:** `CSharpNumerics.Engines.Audio`

---

### AudioBuffer

Core audio container with sample data, sample rate, and channel info.

```csharp
// Create a 1-second mono buffer at 44.1 kHz
var buf = new AudioBuffer(44100, 1, 1.0);

// Wrap existing data
var stereo = new AudioBuffer(new double[] { 0.5, -0.3, 0.8, 0.2 }, 44100, 2);

// Access / properties
double sample = stereo[0, 1];        // frame 0, right channel
int frames = stereo.FrameCount;      // 2
double dur = stereo.Duration;         // ~0.0000453 s

// Operations
var mono = stereo.ToMono();           // average channels
buf.MixIn(otherBuffer, gain: 0.5);   // additive mix
buf.Normalize();                      // peak → 1.0
```

### SignalGenerator

Generate standard waveforms: Sine, Square, Sawtooth, Triangle, WhiteNoise.

```csharp
var sine = SignalGenerator.Generate(SignalGenerator.Waveform.Sine, 440, 1.0, 2.0);
var square = SignalGenerator.Generate(SignalGenerator.Waveform.Square, 220, 0.8, 1.0, sampleRate: 48000);
var noise = SignalGenerator.Generate(SignalGenerator.Waveform.WhiteNoise, 0, 0.5, 0.5);
```

### AudioOscillator

Stateful oscillator with continuous phase — click-free, suitable for real-time synthesis.

```csharp
var osc = new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 0.8);

// Sample-by-sample
double s = osc.NextSample(44100);

// Batch generate
var buf = osc.GenerateBuffer(duration: 1.0, sampleRate: 44100);

// Modulation: change frequency on the fly
osc.Frequency = 880;
```

### Envelope (ADSR)

Attack-Decay-Sustain-Release amplitude shaping.

```csharp
var env = new Envelope(attack: 0.05, decay: 0.1, sustain: 0.7, release: 0.3);

double amp = env.Evaluate(t: 0.02);                  // during attack
double ampR = env.Evaluate(t: 0.8, noteOffTime: 0.5); // during release

// Apply to a buffer (note released at 0.5 s)
env.Apply(buffer, noteOffTime: 0.5);
```

### Synthesizer

Additive synthesis: combine multiple oscillators + ADSR envelope → AudioBuffer.

```csharp
var synth = new Synthesizer { SampleRate = 44100 };
synth.AddVoice(new AudioOscillator(SignalGenerator.Waveform.Sine, 440, 1.0));
synth.AddVoice(new AudioOscillator(SignalGenerator.Waveform.Sine, 880, 0.5), gain: 0.3);
synth.Envelope = new Envelope(0.01, 0.05, 0.8, 0.2);

var buf = synth.Render(duration: 1.0, noteOffTime: 0.8);
```

---

### AudioFilter

Frequency-domain filtering via FFT → mask → IFFT.

```csharp
var lp = AudioFilter.Apply(buffer, AudioFilter.FilterType.LowPass, cutoffLow: 1000);
var hp = AudioFilter.Apply(buffer, AudioFilter.FilterType.HighPass, 0, cutoffHigh: 500);
var bp = AudioFilter.Apply(buffer, AudioFilter.FilterType.BandPass, cutoffLow: 200, cutoffHigh: 4000);
```

### Reverb

Schroeder-model reverb: parallel comb filters + series all-pass filters.

```csharp
var reverb = new Reverb(roomSize: 0.7, damping: 0.4, wetMix: 0.3);
var wet = reverb.Process(buffer);
```

### Delay

Circular-buffer delay with feedback.

```csharp
var delay = new Delay(delayTime: 0.25, feedback: 0.6, wetMix: 0.4);
var echoed = delay.Process(buffer);
```

### Compressor

Dynamic range compression with attack/release envelope.

```csharp
var comp = new Compressor(threshold: 0.5, ratio: 4.0, attack: 0.01, release: 0.1);
var compressed = comp.Process(buffer);
```

### SpatialAudio

Stereo panning and distance attenuation.

```csharp
// Constant-power panning (-1 left, 0 center, +1 right)
var stereo = SpatialAudio.Pan(monoBuffer, pan: -0.3);

// Inverse-distance attenuation
var attenuated = SpatialAudio.AttenuateByDistance(buffer, distance: 5.0);

// Combined: position-based spatialization
var spatial = SpatialAudio.Spatialize(monoBuffer,
    sourceX: 10, sourceY: 5, listenerX: 0, listenerY: 0);
```

---

### SpectrumAnalyzer

Windowed FFT analysis with Hann, Hamming, Blackman, or rectangular windows.

```csharp
var analyzer = new SpectrumAnalyzer(fftSize: 2048, window: SpectrumAnalyzer.WindowType.Hann);

// Single frame → (frequency, magnitude) pairs
var spectrum = analyzer.Analyze(samples, sampleRate: 44100);

// Averaged over whole buffer (50% overlap)
var avgSpectrum = analyzer.AnalyzeBuffer(buffer);
```

### PitchDetector

Fundamental frequency detection via autocorrelation or Harmonic Product Spectrum.

```csharp
var detector = new PitchDetector(fftSize: 4096);
double f0 = detector.Detect(buffer, PitchDetector.Method.Autocorrelation);

// HPS — best with harmonic-rich signals
double f0hps = detector.Detect(buffer, PitchDetector.Method.HarmonicProductSpectrum);
```

### BeatDetector

Onset detection via spectral flux, with tempo estimation.

```csharp
var beat = new BeatDetector(frameSize: 1024) { Threshold = 1.5 };
List<double> onsets = beat.Detect(buffer);       // onset times in seconds
double bpm = beat.EstimateTempo(buffer);          // estimated BPM
```
