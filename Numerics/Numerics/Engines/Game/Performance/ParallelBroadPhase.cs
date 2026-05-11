using CSharpNumerics.Physics.Mechanics.Objects;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace CSharpNumerics.Engines.Game.Performance;

/// <summary>
/// Multi-threaded broad phase that partitions the body array across threads.
///
/// Splits bodies into spatial chunks along the X-axis, runs broad phase
/// on each chunk in parallel, then cross-tests boundary pairs sequentially.
///
/// Uses <see cref="Task.WhenAll"/> for parallelism.
/// Falls back to single-threaded when body count is below a threshold.
/// </summary>
public class ParallelBroadPhase : BroadPhase.IBroadPhase
{
    private readonly int _threadCount;
    private readonly int _singleThreadThreshold;
    private int[] _sorted;

    /// <summary>
    /// Creates a parallel broad phase.
    /// </summary>
    /// <param name="threadCount">Number of parallel workers. Default: processor count.</param>
    /// <param name="singleThreadThreshold">Below this body count, run single-threaded. Default 128.</param>
    public ParallelBroadPhase(int threadCount = 0, int singleThreadThreshold = 128)
    {
        _threadCount = threadCount > 0 ? threadCount : Environment.ProcessorCount;
        _singleThreadThreshold = singleThreadThreshold;
        _sorted = Array.Empty<int>();
    }

    /// <inheritdoc/>
    public void FindPairs(RigidBody[] bodies, double[] radii, int count, List<(int a, int b)> results)
    {
        results.Clear();
        if (count < 2) return;

        // Below threshold: single-threaded sweep-and-prune
        if (count < _singleThreadThreshold)
        {
            SingleThreadSweep(bodies, radii, count, results);
            return;
        }

        // Sort all bodies by X
        if (_sorted.Length < count)
            _sorted = new int[count];
        for (int i = 0; i < count; i++)
            _sorted[i] = i;

        Array.Sort(_sorted, 0, count, Comparer<int>.Create((a, b) =>
            (bodies[a].Position.x - radii[a]).CompareTo(bodies[b].Position.x - radii[b])));

        // Partition into chunks
        int chunkSize = (count + _threadCount - 1) / _threadCount;
        var chunkResults = new List<(int a, int b)>[_threadCount];
        var tasks = new Task[Math.Min(_threadCount, (count + chunkSize - 1) / chunkSize)];

        for (int t = 0; t < tasks.Length; t++)
        {
            int start = t * chunkSize;
            int end = Math.Min(start + chunkSize, count);
            chunkResults[t] = new List<(int a, int b)>();
            var localResults = chunkResults[t];
            int localStart = start;
            int localEnd = end;

            tasks[t] = Task.Run(() =>
                SweepChunk(bodies, radii, count, _sorted, localStart, localEnd, localResults));
        }

        Task.WaitAll(tasks);

        // Merge results
        for (int t = 0; t < tasks.Length; t++)
            results.AddRange(chunkResults[t]);
    }

    private void SweepChunk(RigidBody[] bodies, double[] radii, int totalCount,
        int[] sorted, int start, int end, List<(int a, int b)> results)
    {
        for (int i = start; i < end; i++)
        {
            int a = sorted[i];
            double maxXa = bodies[a].Position.x + radii[a];

            // Search forward from i+1 through the entire sorted array
            // (not just our chunk — we need cross-chunk pairs too)
            for (int j = i + 1; j < totalCount; j++)
            {
                int b = sorted[j];
                double minXb = bodies[b].Position.x - radii[b];

                if (minXb > maxXa) break;

                if (bodies[a].IsStatic && bodies[b].IsStatic) continue;

                double rSum = radii[a] + radii[b];
                if (Math.Abs(bodies[a].Position.y - bodies[b].Position.y) > rSum) continue;
                if (Math.Abs(bodies[a].Position.z - bodies[b].Position.z) > rSum) continue;

                var d = bodies[b].Position - bodies[a].Position;
                if (d.Dot(d) <= rSum * rSum)
                {
                    lock (results)
                        results.Add((a, b));
                }
            }
        }
    }

    private void SingleThreadSweep(RigidBody[] bodies, double[] radii, int count,
        List<(int a, int b)> results)
    {
        if (_sorted.Length < count)
            _sorted = new int[count];
        for (int i = 0; i < count; i++)
            _sorted[i] = i;

        Array.Sort(_sorted, 0, count, Comparer<int>.Create((a, b) =>
            (bodies[a].Position.x - radii[a]).CompareTo(bodies[b].Position.x - radii[b])));

        for (int i = 0; i < count; i++)
        {
            int a = _sorted[i];
            double maxXa = bodies[a].Position.x + radii[a];

            for (int j = i + 1; j < count; j++)
            {
                int b = _sorted[j];
                double minXb = bodies[b].Position.x - radii[b];
                if (minXb > maxXa) break;

                if (bodies[a].IsStatic && bodies[b].IsStatic) continue;

                double rSum = radii[a] + radii[b];
                if (Math.Abs(bodies[a].Position.y - bodies[b].Position.y) > rSum) continue;
                if (Math.Abs(bodies[a].Position.z - bodies[b].Position.z) > rSum) continue;

                var d = bodies[b].Position - bodies[a].Position;
                if (d.Dot(d) <= rSum * rSum)
                    results.Add((a, b));
            }
        }
    }
}
