using CSharpNumerics.Physics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Applied.BroadPhase
{
    /// <summary>
    /// Sweep-and-prune broad phase along the X axis.
    /// O(n log n) due to the sort; much faster than brute force for large scenes
    /// where objects are spread out spatially.
    /// </summary>
    public class SweepAndPruneBroadPhase : IBroadPhase
    {
        private int[] _sorted = [];

        public void FindPairs(RigidBody[] bodies, double[] radii, int count, List<(int a, int b)> results)
        {
            results.Clear();
            if (count < 2) return;

            // Ensure sorted index array is large enough
            if (_sorted.Length < count)
                _sorted = new int[count];
            for (int i = 0; i < count; i++)
                _sorted[i] = i;

            // Sort by AABB min X (position.x - radius)
            Array.Sort(_sorted, 0, count, Comparer<int>.Create((a, b) =>
                (bodies[a].Position.x - radii[a]).CompareTo(bodies[b].Position.x - radii[b])));

            // Sweep along X, then filter by Y/Z overlap + full sphere check
            for (int i = 0; i < count; i++)
            {
                int a = _sorted[i];
                double maxXa = bodies[a].Position.x + radii[a];

                for (int j = i + 1; j < count; j++)
                {
                    int b = _sorted[j];
                    double minXb = bodies[b].Position.x - radii[b];

                    // Past the X overlap window — no more candidates for this body
                    if (minXb > maxXa) break;

                    // Skip static-static pairs
                    if (bodies[a].IsStatic && bodies[b].IsStatic) continue;

                    // Quick Y/Z AABB rejection
                    double rSum = radii[a] + radii[b];
                    if (Math.Abs(bodies[a].Position.y - bodies[b].Position.y) > rSum) continue;
                    if (Math.Abs(bodies[a].Position.z - bodies[b].Position.z) > rSum) continue;

                    // Full sphere overlap check
                    var d = bodies[b].Position - bodies[a].Position;
                    if (d.Dot(d) <= rSum * rSum)
                        results.Add((a, b));
                }
            }
        }
    }
}
