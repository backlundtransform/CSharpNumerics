using CSharpNumerics.Numerics.Objects;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Physics.Mechanics.SoftBody;

/// <summary>
/// Mass-spring network on a triangulated mesh for soft-body deformation.
///
/// Each vertex is a point mass connected to neighbours by springs.
/// Springs provide structural (edge), shear (diagonal), and bend (skip-one) forces.
/// Integration uses Verlet for stability without explicit velocity storage.
///
/// The mesh can collide with spherical obstacles via simple projection.
/// </summary>
public class DeformableMesh
{
    /// <summary>Vertex positions.</summary>
    public Vector[] Positions;

    /// <summary>Previous positions (for Verlet integration).</summary>
    public Vector[] PrevPositions;

    /// <summary>Per-vertex inverse mass (0 = pinned).</summary>
    public double[] InverseMass;

    /// <summary>External acceleration (e.g. gravity).</summary>
    public Vector Gravity { get; set; } = new Vector(0, 0, -9.81);

    /// <summary>Velocity damping factor [0,1]. 1 = no damping.</summary>
    public double Damping { get; set; } = 0.99;

    /// <summary>Number of constraint iterations per step.</summary>
    public int Iterations { get; set; } = 5;

    private readonly List<Spring> _springs = new();
    private readonly int _vertexCount;

    /// <summary>Number of vertices.</summary>
    public int VertexCount => _vertexCount;

    /// <summary>Number of springs.</summary>
    public int SpringCount => _springs.Count;

    /// <summary>
    /// Creates a deformable mesh from vertex positions and a mass per vertex.
    /// </summary>
    public DeformableMesh(Vector[] positions, double mass)
    {
        _vertexCount = positions.Length;
        Positions = new Vector[_vertexCount];
        PrevPositions = new Vector[_vertexCount];
        InverseMass = new double[_vertexCount];

        double invM = mass > 0 ? 1.0 / mass : 0;
        for (int i = 0; i < _vertexCount; i++)
        {
            Positions[i] = positions[i];
            PrevPositions[i] = positions[i];
            InverseMass[i] = invM;
        }
    }

    /// <summary>Pin a vertex so it cannot move.</summary>
    public void Pin(int index) => InverseMass[index] = 0;

    /// <summary>Unpin a vertex.</summary>
    public void Unpin(int index, double mass) => InverseMass[index] = mass > 0 ? 1.0 / mass : 0;

    /// <summary>
    /// Add a structural spring between two vertices.
    /// Rest length is computed from current positions.
    /// </summary>
    public void AddSpring(int a, int b, double stiffness)
    {
        double restLength = (Positions[a] - Positions[b]).GetMagnitude();
        _springs.Add(new Spring(a, b, restLength, stiffness));
    }

    /// <summary>
    /// Add a spring with explicit rest length.
    /// </summary>
    public void AddSpring(int a, int b, double restLength, double stiffness)
    {
        _springs.Add(new Spring(a, b, restLength, stiffness));
    }

    /// <summary>
    /// Automatically generate springs from triangle indices.
    /// Each edge of each triangle gets a structural spring. Duplicate edges are skipped.
    /// </summary>
    /// <param name="triangles">Flat array of triangle vertex indices (length must be multiple of 3).</param>
    /// <param name="stiffness">Spring stiffness.</param>
    public void GenerateSpringsFromTriangles(int[] triangles, double stiffness)
    {
        var edges = new HashSet<(int, int)>();
        for (int i = 0; i < triangles.Length; i += 3)
        {
            AddEdge(triangles[i], triangles[i + 1], stiffness, edges);
            AddEdge(triangles[i + 1], triangles[i + 2], stiffness, edges);
            AddEdge(triangles[i + 2], triangles[i], stiffness, edges);
        }
    }

    /// <summary>
    /// Advance the simulation by dt using Verlet integration + constraint projection.
    /// </summary>
    public void Step(double dt)
    {
        double dt2 = dt * dt;

        // Verlet integration
        for (int i = 0; i < _vertexCount; i++)
        {
            if (InverseMass[i] == 0) continue;

            var vel = Positions[i] - PrevPositions[i];
            var newPos = Positions[i] + Damping * vel + dt2 * Gravity;
            PrevPositions[i] = Positions[i];
            Positions[i] = newPos;
        }

        // Constraint projection (springs)
        for (int iter = 0; iter < Iterations; iter++)
        {
            for (int s = 0; s < _springs.Count; s++)
            {
                var sp = _springs[s];
                var delta = Positions[sp.B] - Positions[sp.A];
                double dist = delta.GetMagnitude();
                if (dist < 1e-15) continue;

                double diff = (dist - sp.RestLength) / dist;
                double wA = InverseMass[sp.A];
                double wB = InverseMass[sp.B];
                double wSum = wA + wB;
                if (wSum < 1e-15) continue;

                double correction = diff * sp.Stiffness;
                var offset = correction / wSum * delta;

                if (wA > 0)
                    Positions[sp.A] = Positions[sp.A] + wA * offset;
                if (wB > 0)
                    Positions[sp.B] = Positions[sp.B] - wB * offset;
            }
        }
    }

    /// <summary>
    /// Resolve collision with a sphere (project vertices out of the sphere).
    /// </summary>
    public void CollideWithSphere(Vector center, double radius)
    {
        for (int i = 0; i < _vertexCount; i++)
        {
            if (InverseMass[i] == 0) continue;

            var d = Positions[i] - center;
            double dist = d.GetMagnitude();
            if (dist < radius && dist > 1e-15)
            {
                Positions[i] = center + (radius / dist) * d;
            }
        }
    }

    /// <summary>
    /// Resolve collision with a ground plane at z = height.
    /// </summary>
    public void CollideWithGround(double height = 0)
    {
        for (int i = 0; i < _vertexCount; i++)
        {
            if (InverseMass[i] == 0) continue;
            if (Positions[i].z < height)
                Positions[i] = new Vector(Positions[i].x, Positions[i].y, height);
        }
    }

    /// <summary>
    /// Create a rectangular cloth mesh (grid of vertices with triangle springs).
    /// </summary>
    /// <param name="width">Width in world units.</param>
    /// <param name="height">Height in world units.</param>
    /// <param name="resX">Number of vertices along X.</param>
    /// <param name="resY">Number of vertices along Y.</param>
    /// <param name="mass">Mass per vertex.</param>
    /// <param name="stiffness">Spring stiffness.</param>
    /// <param name="origin">Bottom-left corner position.</param>
    public static DeformableMesh CreateGrid(double width, double height, int resX, int resY,
        double mass, double stiffness, Vector origin = default)
    {
        var verts = new Vector[resX * resY];
        double dx = width / (resX - 1);
        double dy = height / (resY - 1);

        for (int j = 0; j < resY; j++)
            for (int i = 0; i < resX; i++)
                verts[i + j * resX] = origin + new Vector(i * dx, j * dy, 0);

        var mesh = new DeformableMesh(verts, mass);

        // Structural springs (horizontal + vertical)
        for (int j = 0; j < resY; j++)
        {
            for (int i = 0; i < resX; i++)
            {
                int idx = i + j * resX;
                if (i < resX - 1) mesh.AddSpring(idx, idx + 1, stiffness);
                if (j < resY - 1) mesh.AddSpring(idx, idx + resX, stiffness);
                // Shear springs (diagonals)
                if (i < resX - 1 && j < resY - 1)
                {
                    mesh.AddSpring(idx, idx + resX + 1, stiffness * 0.5);
                    mesh.AddSpring(idx + 1, idx + resX, stiffness * 0.5);
                }
                // Bend springs (skip one)
                if (i < resX - 2) mesh.AddSpring(idx, idx + 2, stiffness * 0.3);
                if (j < resY - 2) mesh.AddSpring(idx, idx + 2 * resX, stiffness * 0.3);
            }
        }

        return mesh;
    }

    private void AddEdge(int a, int b, double stiffness, HashSet<(int, int)> edges)
    {
        var key = a < b ? (a, b) : (b, a);
        if (edges.Add(key))
            AddSpring(a, b, stiffness);
    }

    private struct Spring
    {
        public int A, B;
        public double RestLength;
        public double Stiffness;

        public Spring(int a, int b, double restLength, double stiffness)
        {
            A = a; B = b; RestLength = restLength; Stiffness = stiffness;
        }
    }
}
