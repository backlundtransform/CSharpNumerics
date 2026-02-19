using CSharpNumerics.Physics.Objects;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Applied.BroadPhase
{
    /// <summary>
    /// O(n²) brute-force broad phase. Tests every pair.
    /// Simple and allocation-free — suitable for small scenes (&lt; 100 bodies).
    /// </summary>
    public class BruteForceBroadPhase : IBroadPhase
    {
        public void FindPairs(RigidBody[] bodies, double[] radii, int count, List<(int a, int b)> results)
        {
            results.Clear();
            for (int i = 0; i < count; i++)
            {
                for (int j = i + 1; j < count; j++)
                {
                    // Skip static-static pairs (neither can move)
                    if (bodies[i].IsStatic && bodies[j].IsStatic) continue;

                    var d = bodies[j].Position - bodies[i].Position;
                    double distSq = d.Dot(d);
                    double rSum = radii[i] + radii[j];

                    if (distSq <= rSum * rSum)
                        results.Add((i, j));
                }
            }
        }
    }
}
