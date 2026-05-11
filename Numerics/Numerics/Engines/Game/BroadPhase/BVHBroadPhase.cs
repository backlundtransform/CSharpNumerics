using CSharpNumerics.Physics.Mechanics.Objects;
using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Engines.Game.BroadPhase;

/// <summary>
/// Bounding Volume Hierarchy broad phase using a binary tree of AABBs.
///
/// Top-down construction: recursively splits bodies along the longest axis
/// of the enclosing AABB using the median position.
///
/// Query: traverse the tree and collect overlapping leaf pairs.
/// Rebuilt each frame (suitable for dynamic scenes at moderate body counts).
/// </summary>
public class BVHBroadPhase : IBroadPhase
{
    private int[] _nodeBodyIndex;  // -1 for internal nodes, body index for leaves
    private Vector[] _nodeMin;
    private Vector[] _nodeMax;
    private int[] _nodeLeft;       // left child index, -1 for leaf
    private int[] _nodeRight;      // right child index, -1 for leaf
    private int _nodeCount;
    private int _capacity;

    public BVHBroadPhase(int initialCapacity = 256)
    {
        _capacity = initialCapacity;
        _nodeBodyIndex = new int[_capacity];
        _nodeMin = new Vector[_capacity];
        _nodeMax = new Vector[_capacity];
        _nodeLeft = new int[_capacity];
        _nodeRight = new int[_capacity];
    }

    public void FindPairs(RigidBody[] bodies, double[] radii, int count, List<(int a, int b)> results)
    {
        results.Clear();
        if (count < 2) return;

        // Build BVH
        _nodeCount = 0;
        var indices = new int[count];
        for (int i = 0; i < count; i++) indices[i] = i;

        EnsureCapacity(2 * count);
        BuildNode(bodies, radii, indices, 0, count);

        // Query all leaf-leaf pairs
        CollectPairs(0, bodies, radii, results);
    }

    private int BuildNode(RigidBody[] bodies, double[] radii, int[] indices, int start, int end)
    {
        int nodeIdx = _nodeCount++;
        EnsureCapacity(_nodeCount);

        // Compute AABB for all bodies in [start, end)
        double minX = double.MaxValue, minY = double.MaxValue, minZ = double.MaxValue;
        double maxX = double.MinValue, maxY = double.MinValue, maxZ = double.MinValue;

        for (int i = start; i < end; i++)
        {
            int b = indices[i];
            double r = radii[b];
            var p = bodies[b].Position;
            minX = Math.Min(minX, p.x - r);
            minY = Math.Min(minY, p.y - r);
            minZ = Math.Min(minZ, p.z - r);
            maxX = Math.Max(maxX, p.x + r);
            maxY = Math.Max(maxY, p.y + r);
            maxZ = Math.Max(maxZ, p.z + r);
        }

        _nodeMin[nodeIdx] = new Vector(minX, minY, minZ);
        _nodeMax[nodeIdx] = new Vector(maxX, maxY, maxZ);

        int count = end - start;
        if (count == 1)
        {
            // Leaf
            _nodeBodyIndex[nodeIdx] = indices[start];
            _nodeLeft[nodeIdx] = -1;
            _nodeRight[nodeIdx] = -1;
            return nodeIdx;
        }

        _nodeBodyIndex[nodeIdx] = -1;

        // Split along longest axis
        double extX = maxX - minX;
        double extY = maxY - minY;
        double extZ = maxZ - minZ;
        int axis = extX >= extY && extX >= extZ ? 0 : extY >= extZ ? 1 : 2;

        // Sort by axis (simple insertion sort for small arrays, works well for game body counts)
        Array.Sort(indices, start, count, Comparer<int>.Create((a, b) =>
        {
            double va = axis == 0 ? bodies[a].Position.x : axis == 1 ? bodies[a].Position.y : bodies[a].Position.z;
            double vb = axis == 0 ? bodies[b].Position.x : axis == 1 ? bodies[b].Position.y : bodies[b].Position.z;
            return va.CompareTo(vb);
        }));

        int mid = start + count / 2;

        _nodeLeft[nodeIdx] = BuildNode(bodies, radii, indices, start, mid);
        _nodeRight[nodeIdx] = BuildNode(bodies, radii, indices, mid, end);

        return nodeIdx;
    }

    private void CollectPairs(int nodeIdx, RigidBody[] bodies, double[] radii, List<(int a, int b)> results)
    {
        if (_nodeLeft[nodeIdx] == -1) return; // leaf has nothing to self-test

        int left = _nodeLeft[nodeIdx];
        int right = _nodeRight[nodeIdx];

        // Recurse into children
        CollectPairs(left, bodies, radii, results);
        CollectPairs(right, bodies, radii, results);

        // Cross-test left vs right subtrees
        CrossTest(left, right, bodies, radii, results);
    }

    private void CrossTest(int nodeA, int nodeB, RigidBody[] bodies, double[] radii, List<(int a, int b)> results)
    {
        // AABB overlap test
        if (_nodeMin[nodeA].x > _nodeMax[nodeB].x || _nodeMax[nodeA].x < _nodeMin[nodeB].x) return;
        if (_nodeMin[nodeA].y > _nodeMax[nodeB].y || _nodeMax[nodeA].y < _nodeMin[nodeB].y) return;
        if (_nodeMin[nodeA].z > _nodeMax[nodeB].z || _nodeMax[nodeA].z < _nodeMin[nodeB].z) return;

        bool leafA = _nodeLeft[nodeA] == -1;
        bool leafB = _nodeLeft[nodeB] == -1;

        if (leafA && leafB)
        {
            int a = _nodeBodyIndex[nodeA];
            int b = _nodeBodyIndex[nodeB];

            if (bodies[a].IsStatic && bodies[b].IsStatic) return;

            // Sphere overlap check
            var d = bodies[b].Position - bodies[a].Position;
            double rSum = radii[a] + radii[b];
            if (d.Dot(d) <= rSum * rSum)
                results.Add((a, b));
            return;
        }

        if (leafA)
        {
            CrossTest(nodeA, _nodeLeft[nodeB], bodies, radii, results);
            CrossTest(nodeA, _nodeRight[nodeB], bodies, radii, results);
        }
        else if (leafB)
        {
            CrossTest(_nodeLeft[nodeA], nodeB, bodies, radii, results);
            CrossTest(_nodeRight[nodeA], nodeB, bodies, radii, results);
        }
        else
        {
            CrossTest(_nodeLeft[nodeA], _nodeLeft[nodeB], bodies, radii, results);
            CrossTest(_nodeLeft[nodeA], _nodeRight[nodeB], bodies, radii, results);
            CrossTest(_nodeRight[nodeA], _nodeLeft[nodeB], bodies, radii, results);
            CrossTest(_nodeRight[nodeA], _nodeRight[nodeB], bodies, radii, results);
        }
    }

    private void EnsureCapacity(int needed)
    {
        if (needed <= _capacity) return;
        _capacity = Math.Max(_capacity * 2, needed);
        Array.Resize(ref _nodeBodyIndex, _capacity);
        Array.Resize(ref _nodeMin, _capacity);
        Array.Resize(ref _nodeMax, _capacity);
        Array.Resize(ref _nodeLeft, _capacity);
        Array.Resize(ref _nodeRight, _capacity);
    }
}
