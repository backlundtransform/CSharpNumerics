using System;
using System.Collections.Generic;
using CSharpNumerics.Engines.Common;
using CSharpNumerics.Engines.Exoplanet.Data;
using CSharpNumerics.Engines.Exoplanet.Pipeline;

namespace CSharpNumerics.Engines.Exoplanet;

/// <summary>
/// Simulation engine for exoplanet transit detection and classification.
/// Implements <see cref="ISimulationEngine"/>: call <see cref="Init"/> to prepare the engine,
/// then <see cref="Enqueue(LightCurve)"/> light curves and <see cref="Step(double)"/> to process them.
/// Detected transits are published as <see cref="TransitDetectedEvent"/> via the <see cref="Bus"/>.
/// </summary>
public class ExoplanetEngine : ISimulationEngine
{
    private readonly Queue<LightCurve> _queue = new Queue<LightCurve>();
    private TransitInferencePipeline _inference;

    /// <summary>Engine configuration.</summary>
    public ExoplanetEngineConfig Config { get; }

    /// <summary>Event bus for publishing <see cref="TransitDetectedEvent"/>.</summary>
    public EventBus Bus { get; } = new EventBus();

    /// <inheritdoc />
    public double Time { get; private set; }

    /// <inheritdoc />
    public bool IsInitialized { get; private set; }

    /// <summary>Total number of light curves processed since last <see cref="Reset"/>.</summary>
    public int ProcessedCount { get; private set; }

    /// <summary>All transit candidates detected since last <see cref="Reset"/>.</summary>
    public List<TransitCandidate> Detections { get; } = new List<TransitCandidate>();

    public ExoplanetEngine(ExoplanetEngineConfig config)
    {
        Config = config ?? throw new ArgumentNullException(nameof(config));
    }

    /// <summary>
    /// Initializes the engine for classical (non-ML) detection.
    /// To use ML-assisted classification, call <see cref="Init(TrainedTransitModel)"/> instead,
    /// or use <see cref="ModelSerializer.LoadFromFile"/> to load a pre-trained model from disk.
    /// </summary>
    public void Init()
    {
        IsInitialized = true;
    }

    /// <summary>
    /// Initializes the engine with a pre-trained model for ML-assisted classification.
    /// </summary>
    public void Init(TrainedTransitModel trainedModel)
    {
        if (trainedModel != null)
            _inference = TransitInferencePipeline.FromModel(trainedModel);

        IsInitialized = true;
    }

    /// <summary>
    /// Enqueue a light curve for processing on the next <see cref="Step"/> call.
    /// </summary>
    public void Enqueue(LightCurve lc)
    {
        if (lc == null) throw new ArgumentNullException(nameof(lc));
        _queue.Enqueue(lc);
    }

    /// <summary>
    /// Advances simulation time by <paramref name="dt"/> seconds and processes all queued light curves.
    /// For each light curve, runs transit detection (and ML classification if a model is loaded).
    /// Publishes a <see cref="TransitDetectedEvent"/> for each detected candidate.
    /// </summary>
    public void Step(double dt)
    {
        if (!IsInitialized)
            throw new InvalidOperationException("Engine must be initialized before stepping.");

        Time += dt;

        while (_queue.Count > 0)
        {
            var lc = _queue.Dequeue();
            ProcessLightCurve(lc);
            ProcessedCount++;
        }
    }

    /// <inheritdoc />
    public void Reset()
    {
        Time = 0;
        ProcessedCount = 0;
        _queue.Clear();
        Detections.Clear();
        Bus.Clear();
        IsInitialized = false;
    }

    private void ProcessLightCurve(LightCurve lc)
    {
        if (_inference != null)
        {
            // ML-assisted path: detection + classification
            var predictions = _inference.Predict(lc);
            foreach (var pred in predictions)
            {
                Detections.Add(pred.Candidate);
                Bus.Publish(new TransitDetectedEvent(pred.Candidate, Time));
            }
        }
        else
        {
            // Classical detection path
            var candidates = TransitDetectionPipeline.Detect(lc, Config.Detection);
            foreach (var candidate in candidates)
            {
                Detections.Add(candidate);
                Bus.Publish(new TransitDetectedEvent(candidate, Time));
            }
        }
    }
}
