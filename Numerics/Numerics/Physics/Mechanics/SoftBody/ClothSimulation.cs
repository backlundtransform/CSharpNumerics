using CSharpNumerics.Numerics.Objects;
using System;

namespace CSharpNumerics.Physics.Mechanics.SoftBody;

/// <summary>
/// Cloth simulation using a constrained particle system.
/// Built on top of <see cref="DeformableMesh"/> with additional features:
/// self-collision avoidance and wind force application.
///
/// The cloth is a rectangular grid of particles connected by structural,
/// shear, and bend springs. Verlet integration with constraint projection.
/// </summary>
public class ClothSimulation
{
    private readonly DeformableMesh _mesh;
    private readonly int _resX, _resY;

    /// <summary>The underlying deformable mesh.</summary>
    public DeformableMesh Mesh => _mesh;

    /// <summary>Resolution along X (columns).</summary>
    public int ResX => _resX;

    /// <summary>Resolution along Y (rows).</summary>
    public int ResY => _resY;

    /// <summary>Wind force applied to each particle.</summary>
    public Vector Wind { get; set; } = new Vector(0, 0, 0);

    /// <summary>Self-collision distance threshold.</summary>
    public double SelfCollisionRadius { get; set; } = 0.05;

    /// <summary>Whether self-collision is enabled.</summary>
    public bool EnableSelfCollision { get; set; } = false;

    /// <summary>
    /// Creates a cloth simulation.
    /// </summary>
    /// <param name="width">Cloth width in world units.</param>
    /// <param name="height">Cloth height in world units.</param>
    /// <param name="resX">Number of particles along X.</param>
    /// <param name="resY">Number of particles along Y.</param>
    /// <param name="mass">Mass per particle.</param>
    /// <param name="stiffness">Spring stiffness.</param>
    /// <param name="origin">Top-left corner (cloth hangs in -Z).</param>
    public ClothSimulation(double width, double height, int resX, int resY,
        double mass = 0.1, double stiffness = 0.9, Vector origin = default)
    {
        _resX = resX;
        _resY = resY;
        _mesh = DeformableMesh.CreateGrid(width, height, resX, resY, mass, stiffness, origin);
    }

    /// <summary>Pin the top row of the cloth (fixed edge).</summary>
    public void PinTopEdge()
    {
        int topRow = _resY - 1;
        for (int i = 0; i < _resX; i++)
            _mesh.Pin(i + topRow * _resX);
    }

    /// <summary>Pin the top-left and top-right corners only.</summary>
    public void PinTopCorners()
    {
        int topRow = _resY - 1;
        _mesh.Pin(topRow * _resX);              // top-left
        _mesh.Pin(topRow * _resX + _resX - 1);  // top-right
    }

    /// <summary>Pin a specific vertex.</summary>
    public void Pin(int x, int y) => _mesh.Pin(x + y * _resX);

    /// <summary>Get the position of a vertex.</summary>
    public Vector GetPosition(int x, int y) => _mesh.Positions[x + y * _resX];

    /// <summary>
    /// Advance the cloth simulation by dt.
    /// </summary>
    public void Step(double dt)
    {
        // Apply wind force as a velocity impulse via position adjustment
        if (Wind.GetMagnitude() > 1e-10)
        {
            double dt2 = dt * dt;
            for (int i = 0; i < _mesh.VertexCount; i++)
            {
                if (_mesh.InverseMass[i] == 0) continue;
                // Wind force proportional to triangle normal (pressure model)
                _mesh.Positions[i] = _mesh.Positions[i] + dt2 * _mesh.InverseMass[i] * Wind;
            }
        }

        _mesh.Step(dt);

        // Self-collision
        if (EnableSelfCollision)
            ResolveSelfCollision();
    }

    /// <summary>
    /// Collide the cloth with a sphere.
    /// </summary>
    public void CollideWithSphere(Vector center, double radius)
    {
        _mesh.CollideWithSphere(center, radius);
    }

    /// <summary>
    /// Collide the cloth with a ground plane.
    /// </summary>
    public void CollideWithGround(double height = 0)
    {
        _mesh.CollideWithGround(height);
    }

    /// <summary>
    /// Compute the total kinetic energy of the cloth.
    /// </summary>
    public double KineticEnergy()
    {
        double ke = 0;
        for (int i = 0; i < _mesh.VertexCount; i++)
        {
            if (_mesh.InverseMass[i] == 0) continue;
            var vel = _mesh.Positions[i] - _mesh.PrevPositions[i];
            double mass = 1.0 / _mesh.InverseMass[i];
            ke += 0.5 * mass * vel.Dot(vel);
        }
        return ke;
    }

    private void ResolveSelfCollision()
    {
        double r2 = SelfCollisionRadius * SelfCollisionRadius;
        for (int i = 0; i < _mesh.VertexCount; i++)
        {
            if (_mesh.InverseMass[i] == 0) continue;
            for (int j = i + 1; j < _mesh.VertexCount; j++)
            {
                if (_mesh.InverseMass[j] == 0) continue;

                // Skip adjacent vertices (within 2 grid cells)
                int ix = i % _resX, iy = i / _resX;
                int jx = j % _resX, jy = j / _resX;
                if (Math.Abs(ix - jx) <= 1 && Math.Abs(iy - jy) <= 1) continue;

                var delta = _mesh.Positions[j] - _mesh.Positions[i];
                double dist2 = delta.Dot(delta);
                if (dist2 < r2 && dist2 > 1e-20)
                {
                    double dist = Math.Sqrt(dist2);
                    double overlap = SelfCollisionRadius - dist;
                    var correction = (overlap / (2 * dist)) * delta;
                    _mesh.Positions[i] = _mesh.Positions[i] - correction;
                    _mesh.Positions[j] = _mesh.Positions[j] + correction;
                }
            }
        }
    }
}
