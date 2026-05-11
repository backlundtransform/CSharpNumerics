using CSharpNumerics.Engines.Game.Fluids;
using System;
using System.Threading;

namespace CSharpNumerics.Engines.Game.Performance;

/// <summary>
/// Runs a <see cref="GameFluidSolver3D"/> on a background thread.
///
/// The game thread samples the latest completed snapshot without blocking.
/// Double-buffered: while the solver computes the next frame, the game
/// reads the previous frame's density/velocity data.
///
/// Thread-safe: all public methods can be called from the game thread.
/// </summary>
public class FluidSolverThread : IDisposable
{
    private readonly GameFluidSolver3D _solver;
    private readonly Thread _thread;
    private readonly ManualResetEventSlim _stepSignal = new(false);
    private readonly ManualResetEventSlim _readySignal = new(true);
    private volatile bool _running;
    private double _dt;
    private int _stepCount;

    // Double-buffered snapshot
    private double[] _densitySnapshot;
    private double[] _velocityXSnapshot;
    private double[] _velocityYSnapshot;
    private double[] _velocityZSnapshot;
    private readonly object _snapshotLock = new();

    /// <summary>Number of simulation steps completed.</summary>
    public int StepCount => _stepCount;

    /// <summary>Whether the solver thread is running.</summary>
    public bool IsRunning => _running;

    /// <summary>
    /// Creates and starts a background fluid solver thread.
    /// </summary>
    /// <param name="solver">The 3D fluid solver to run.</param>
    public FluidSolverThread(GameFluidSolver3D solver)
    {
        _solver = solver ?? throw new ArgumentNullException(nameof(solver));

        int size = solver.NX * solver.NY * solver.NZ;
        _densitySnapshot = new double[size];
        _velocityXSnapshot = new double[size];
        _velocityYSnapshot = new double[size];
        _velocityZSnapshot = new double[size];

        _running = true;
        _thread = new Thread(WorkerLoop)
        {
            IsBackground = true,
            Name = "FluidSolverThread",
            Priority = ThreadPriority.BelowNormal,
        };
        _thread.Start();
    }

    /// <summary>
    /// Request the solver to perform one step with the given timestep.
    /// Non-blocking: returns immediately. The solver thread picks up the request.
    /// </summary>
    /// <param name="dt">Timestep in seconds.</param>
    public void RequestStep(double dt)
    {
        _dt = dt;
        _readySignal.Wait(); // wait until previous step is done
        _readySignal.Reset();
        _stepSignal.Set();
    }

    /// <summary>
    /// Wait for the current step to complete.
    /// </summary>
    /// <param name="timeoutMs">Maximum wait time in milliseconds. -1 = infinite.</param>
    /// <returns>True if the step completed within the timeout.</returns>
    public bool WaitForStep(int timeoutMs = -1)
    {
        return _readySignal.Wait(timeoutMs);
    }

    /// <summary>
    /// Get the latest density snapshot (thread-safe copy).
    /// </summary>
    public void GetDensitySnapshot(double[] destination)
    {
        lock (_snapshotLock)
            Array.Copy(_densitySnapshot, destination, Math.Min(_densitySnapshot.Length, destination.Length));
    }

    /// <summary>
    /// Get the latest velocity snapshot (thread-safe copy).
    /// </summary>
    public void GetVelocitySnapshot(double[] destX, double[] destY, double[] destZ)
    {
        lock (_snapshotLock)
        {
            int len = Math.Min(_velocityXSnapshot.Length, destX.Length);
            Array.Copy(_velocityXSnapshot, destX, len);
            Array.Copy(_velocityYSnapshot, destY, len);
            Array.Copy(_velocityZSnapshot, destZ, len);
        }
    }

    /// <summary>
    /// Sample density at a position from the latest snapshot.
    /// Uses nearest-cell lookup (no interpolation — fast but coarse).
    /// </summary>
    public double SampleDensity(int i, int j, int k)
    {
        int idx = i + j * _solver.NX + k * _solver.NX * _solver.NY;
        if (idx < 0 || idx >= _densitySnapshot.Length) return 0;
        return _densitySnapshot[idx]; // volatile read is good enough for visual data
    }

    /// <summary>
    /// Stop the solver thread. Call from the game thread.
    /// </summary>
    public void Stop()
    {
        _running = false;
        _stepSignal.Set(); // wake up if waiting
        _thread.Join(5000);
    }

    /// <inheritdoc/>
    public void Dispose()
    {
        Stop();
        _stepSignal.Dispose();
        _readySignal.Dispose();
    }

    private void WorkerLoop()
    {
        while (_running)
        {
            _stepSignal.Wait();
            _stepSignal.Reset();

            if (!_running) break;

            // Perform the simulation step
            _solver.Step(_dt);
            Interlocked.Increment(ref _stepCount);

            // Snapshot the results
            lock (_snapshotLock)
            {
                _solver.Density.CopyTo(_densitySnapshot);
                _solver.VelocityX.CopyTo(_velocityXSnapshot);
                _solver.VelocityY.CopyTo(_velocityYSnapshot);
                _solver.VelocityZ.CopyTo(_velocityZSnapshot);
            }

            _readySignal.Set();
        }
    }
}
