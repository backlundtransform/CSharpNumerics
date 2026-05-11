using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics.Objects;
using System;

namespace CSharpNumerics.Engines.Game.Unity;

/// <summary>
/// Synchronizes a <see cref="PhysicsWorld"/> state with Unity Transforms.
///
/// Each frame, the sync layer reads body positions/orientations from the physics
/// world and writes them to arrays that a Unity MonoBehaviour can consume.
/// Supports interpolation between physics steps for smooth rendering.
///
/// This class does not depend on Unity — it produces float arrays that the
/// Unity shim reads.
/// </summary>
public class PhysicsSync
{
    private readonly PhysicsWorld _world;

    // Snapshot buffers for interpolation
    private UnityAdapter.UnityVector3[] _prevPositions;
    private UnityAdapter.UnityVector3[] _currPositions;
    private UnityAdapter.UnityQuaternion[] _prevRotations;
    private UnityAdapter.UnityQuaternion[] _currRotations;
    private int _capacity;

    /// <summary>Number of synced bodies.</summary>
    public int BodyCount => _world.BodyCount;

    /// <summary>The physics world being synced.</summary>
    public PhysicsWorld World => _world;

    /// <summary>
    /// Creates a physics sync for the given world.
    /// </summary>
    public PhysicsSync(PhysicsWorld world)
    {
        _world = world ?? throw new ArgumentNullException(nameof(world));
        _capacity = Math.Max(world.BodyCount, 16);
        _prevPositions = new UnityAdapter.UnityVector3[_capacity];
        _currPositions = new UnityAdapter.UnityVector3[_capacity];
        _prevRotations = new UnityAdapter.UnityQuaternion[_capacity];
        _currRotations = new UnityAdapter.UnityQuaternion[_capacity];
        SnapshotCurrent();
    }

    /// <summary>
    /// Step the physics world and capture a new snapshot.
    /// Call this from FixedUpdate-equivalent timing.
    /// </summary>
    /// <param name="dt">Physics timestep.</param>
    public void StepAndSync(double dt)
    {
        EnsureCapacity();

        // Shift current → previous
        Array.Copy(_currPositions, _prevPositions, _world.BodyCount);
        Array.Copy(_currRotations, _prevRotations, _world.BodyCount);

        _world.Step(dt);

        SnapshotCurrent();
    }

    /// <summary>
    /// Use the fixed-timestep accumulator and capture snapshots.
    /// </summary>
    /// <param name="elapsed">Wall-clock time since last call.</param>
    /// <returns>Number of physics steps taken.</returns>
    public int UpdateAndSync(double elapsed)
    {
        EnsureCapacity();

        // Shift current → previous
        Array.Copy(_currPositions, _prevPositions, _world.BodyCount);
        Array.Copy(_currRotations, _prevRotations, _world.BodyCount);

        int steps = _world.Update(elapsed);
        SnapshotCurrent();
        return steps;
    }

    /// <summary>
    /// Get the interpolated position for body at given index.
    /// Use <see cref="PhysicsWorld.Alpha"/> for the interpolation factor.
    /// </summary>
    public UnityAdapter.UnityVector3 GetInterpolatedPosition(int bodyIndex, double alpha)
    {
        var prev = _prevPositions[bodyIndex];
        var curr = _currPositions[bodyIndex];
        float a = (float)alpha;
        return new UnityAdapter.UnityVector3(
            prev.x + a * (curr.x - prev.x),
            prev.y + a * (curr.y - prev.y),
            prev.z + a * (curr.z - prev.z));
    }

    /// <summary>
    /// Get the interpolated rotation (SLERP) for body at given index.
    /// </summary>
    public UnityAdapter.UnityQuaternion GetInterpolatedRotation(int bodyIndex, double alpha)
    {
        var q0 = _prevRotations[bodyIndex];
        var q1 = _currRotations[bodyIndex];
        return Slerp(q0, q1, (float)alpha);
    }

    /// <summary>
    /// Get the current (non-interpolated) position in Unity coordinates.
    /// </summary>
    public UnityAdapter.UnityVector3 GetCurrentPosition(int bodyIndex)
    {
        return _currPositions[bodyIndex];
    }

    /// <summary>
    /// Get the current (non-interpolated) rotation in Unity coordinates.
    /// </summary>
    public UnityAdapter.UnityQuaternion GetCurrentRotation(int bodyIndex)
    {
        return _currRotations[bodyIndex];
    }

    /// <summary>
    /// Verify frame coherence: check that body positions haven't jumped unreasonably
    /// between previous and current snapshot.
    /// </summary>
    /// <param name="maxDelta">Maximum allowed position change per step.</param>
    /// <returns>True if all bodies are coherent.</returns>
    public bool CheckCoherence(double maxDelta = 100.0)
    {
        float maxD = (float)maxDelta;
        for (int i = 0; i < _world.BodyCount; i++)
        {
            float dx = _currPositions[i].x - _prevPositions[i].x;
            float dy = _currPositions[i].y - _prevPositions[i].y;
            float dz = _currPositions[i].z - _prevPositions[i].z;
            float dist = (float)Math.Sqrt(dx * dx + dy * dy + dz * dz);
            if (dist > maxD) return false;
        }
        return true;
    }

    private void SnapshotCurrent()
    {
        for (int i = 0; i < _world.BodyCount; i++)
        {
            var body = _world.Body(i);
            _currPositions[i] = UnityAdapter.ToUnityVector3(body.Position);
            _currRotations[i] = UnityAdapter.ToUnityQuaternion(body.Orientation);
        }
    }

    private void EnsureCapacity()
    {
        if (_world.BodyCount > _capacity)
        {
            _capacity = _world.BodyCount * 2;
            Array.Resize(ref _prevPositions, _capacity);
            Array.Resize(ref _currPositions, _capacity);
            Array.Resize(ref _prevRotations, _capacity);
            Array.Resize(ref _currRotations, _capacity);
        }
    }

    private static UnityAdapter.UnityQuaternion Slerp(
        UnityAdapter.UnityQuaternion a, UnityAdapter.UnityQuaternion b, float t)
    {
        float dot = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;

        // Ensure shortest path
        if (dot < 0)
        {
            b = new UnityAdapter.UnityQuaternion(-b.x, -b.y, -b.z, -b.w);
            dot = -dot;
        }

        if (dot > 0.9995f)
        {
            // Linear interpolation for very close quaternions
            return NormalizeQuat(new UnityAdapter.UnityQuaternion(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                a.z + t * (b.z - a.z),
                a.w + t * (b.w - a.w)));
        }

        float theta = (float)Math.Acos(Math.Min(1.0, dot));
        float sinTheta = (float)Math.Sin(theta);
        float wA = (float)Math.Sin((1 - t) * theta) / sinTheta;
        float wB = (float)Math.Sin(t * theta) / sinTheta;

        return new UnityAdapter.UnityQuaternion(
            wA * a.x + wB * b.x,
            wA * a.y + wB * b.y,
            wA * a.z + wB * b.z,
            wA * a.w + wB * b.w);
    }

    private static UnityAdapter.UnityQuaternion NormalizeQuat(UnityAdapter.UnityQuaternion q)
    {
        float mag = (float)Math.Sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
        if (mag < 1e-10f) return UnityAdapter.UnityQuaternion.Identity;
        return new UnityAdapter.UnityQuaternion(q.x / mag, q.y / mag, q.z / mag, q.w / mag);
    }
}
