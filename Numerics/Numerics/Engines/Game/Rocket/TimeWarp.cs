using System;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Configurable simulation time warp (1x–1000x) for interactive rocket simulations.
/// Maintains consistent physics by adjusting the number of sub-steps per render frame
/// rather than increasing the physics timestep.
///
/// At high warp factors, multiple physics steps are computed per visual frame.
/// Interpolation between sub-steps provides smooth rendering.
/// </summary>
public class TimeWarp
{
    /// <summary>Current warp factor (1.0 = real-time). Range [1, MaxWarpFactor].</summary>
    public double WarpFactor
    {
        get => _warpFactor;
        set => _warpFactor = Math.Max(1.0, Math.Min(value, MaxWarpFactor));
    }

    /// <summary>Maximum allowed warp factor. Default 1000.</summary>
    public double MaxWarpFactor { get; set; } = 1000.0;

    /// <summary>Fixed physics timestep (seconds). Default 0.01 (100 Hz).</summary>
    public double PhysicsTimestep { get; set; } = 0.01;

    /// <summary>Render frame interval (seconds). Default 1/60.</summary>
    public double RenderInterval { get; set; } = 1.0 / 60.0;

    /// <summary>Number of physics steps computed in the last frame update.</summary>
    public int StepsLastFrame { get; private set; }

    /// <summary>Interpolation factor (0–1) between last two physics states for smooth rendering.</summary>
    public double InterpolationAlpha { get; private set; }

    /// <summary>Total simulation time advanced (seconds).</summary>
    public double TotalSimTime { get; private set; }

    /// <summary>Total number of physics steps computed.</summary>
    public long TotalSteps { get; private set; }

    private double _warpFactor = 1.0;
    private double _accumulator;

    /// <summary>
    /// Advances the simulation by one render frame. Returns the number of physics steps
    /// that should be executed for this frame.
    /// </summary>
    /// <param name="renderDt">Actual render frame time (seconds). If 0, uses RenderInterval.</param>
    /// <returns>Number of physics sub-steps to execute at PhysicsTimestep each.</returns>
    public int AdvanceFrame(double renderDt = 0)
    {
        if (renderDt <= 0) renderDt = RenderInterval;

        double simDt = renderDt * _warpFactor;
        _accumulator += simDt;

        int steps = (int)(_accumulator / PhysicsTimestep);
        _accumulator -= steps * PhysicsTimestep;

        // Clamp steps to prevent spiral of death at extreme warp + lag
        int maxSteps = (int)(MaxWarpFactor * RenderInterval / PhysicsTimestep) + 1;
        if (steps > maxSteps) steps = maxSteps;

        InterpolationAlpha = _accumulator / PhysicsTimestep;
        StepsLastFrame = steps;
        TotalSteps += steps;
        TotalSimTime += steps * PhysicsTimestep;

        return steps;
    }

    /// <summary>
    /// Steps the simulation engine the appropriate number of sub-steps for this render frame.
    /// </summary>
    /// <param name="engine">The simulation engine to step.</param>
    /// <param name="renderDt">Actual render frame time (seconds).</param>
    public void StepEngine(RocketSimulationEngine engine, double renderDt = 0)
    {
        int steps = AdvanceFrame(renderDt);
        for (int i = 0; i < steps; i++)
        {
            engine.Step(PhysicsTimestep);
        }
    }

    /// <summary>
    /// Increases warp factor to the next standard level (1, 2, 5, 10, 50, 100, 1000).
    /// </summary>
    public void IncreaseWarp()
    {
        double[] levels = { 1, 2, 5, 10, 50, 100, 500, 1000 };
        for (int i = 0; i < levels.Length; i++)
        {
            if (levels[i] > _warpFactor + 0.01)
            {
                WarpFactor = levels[i];
                return;
            }
        }
        WarpFactor = MaxWarpFactor;
    }

    /// <summary>
    /// Decreases warp factor to the previous standard level.
    /// </summary>
    public void DecreaseWarp()
    {
        double[] levels = { 1, 2, 5, 10, 50, 100, 500, 1000 };
        for (int i = levels.Length - 1; i >= 0; i--)
        {
            if (levels[i] < _warpFactor - 0.01)
            {
                WarpFactor = levels[i];
                return;
            }
        }
        WarpFactor = 1.0;
    }

    /// <summary>Resets the time warp to 1x and clears accumulators.</summary>
    public void Reset()
    {
        _warpFactor = 1.0;
        _accumulator = 0;
        StepsLastFrame = 0;
        InterpolationAlpha = 0;
        TotalSimTime = 0;
        TotalSteps = 0;
    }
}
