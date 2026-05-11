using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Physics.Mechanics.Objects;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace CSharpNumerics.Engines.Game.Performance;

/// <summary>
/// Records and replays deterministic simulation state for network sync and replays.
///
/// Captures full state snapshots (body positions, velocities, orientations)
/// at each frame. Snapshots can be serialized to a binary stream and
/// replayed to reproduce the exact simulation.
///
/// Supports:
///   - Record mode: capture state each frame
///   - Playback mode: read state from recorded data
///   - Binary serialization for compact storage
///   - Snapshot-based (no delta compression — simple and reliable)
/// </summary>
public class SimulationRecorder
{
    /// <summary>
    /// Per-body state captured in a snapshot.
    /// </summary>
    public struct BodySnapshot
    {
        public Vector Position;
        public Vector Velocity;
        public Vector AngularVelocity;
        public Matrix Orientation;
    }

    /// <summary>
    /// A complete simulation state at one point in time.
    /// </summary>
    public class FrameSnapshot
    {
        /// <summary>Frame index.</summary>
        public int Frame;

        /// <summary>Simulation time.</summary>
        public double Time;

        /// <summary>Per-body state.</summary>
        public BodySnapshot[] Bodies;
    }

    private readonly List<FrameSnapshot> _frames = new();
    private int _currentFrame;

    /// <summary>Number of recorded frames.</summary>
    public int FrameCount => _frames.Count;

    /// <summary>Current playback frame index.</summary>
    public int CurrentFrame => _currentFrame;

    /// <summary>Total recorded time in seconds.</summary>
    public double Duration => _frames.Count > 0 ? _frames[_frames.Count - 1].Time : 0;

    /// <summary>
    /// Capture the current state of a PhysicsWorld.
    /// Call this once per simulation step during recording.
    /// </summary>
    public void RecordFrame(PhysicsWorld world, double time)
    {
        var snapshot = new FrameSnapshot
        {
            Frame = _frames.Count,
            Time = time,
            Bodies = new BodySnapshot[world.BodyCount],
        };

        for (int i = 0; i < world.BodyCount; i++)
        {
            var body = world.Body(i);
            snapshot.Bodies[i] = new BodySnapshot
            {
                Position = body.Position,
                Velocity = body.Velocity,
                AngularVelocity = body.AngularVelocity,
                Orientation = body.Orientation,
            };
        }

        _frames.Add(snapshot);
    }

    /// <summary>
    /// Apply a recorded frame to a PhysicsWorld (playback).
    /// </summary>
    /// <param name="world">The world to restore state into.</param>
    /// <param name="frameIndex">Frame index to apply.</param>
    public void ApplyFrame(PhysicsWorld world, int frameIndex)
    {
        if (frameIndex < 0 || frameIndex >= _frames.Count)
            throw new ArgumentOutOfRangeException(nameof(frameIndex));

        var snapshot = _frames[frameIndex];
        int count = Math.Min(snapshot.Bodies.Length, world.BodyCount);

        for (int i = 0; i < count; i++)
        {
            ref var body = ref world.Body(i);
            body.Position = snapshot.Bodies[i].Position;
            body.Velocity = snapshot.Bodies[i].Velocity;
            body.AngularVelocity = snapshot.Bodies[i].AngularVelocity;
            body.Orientation = snapshot.Bodies[i].Orientation;
        }

        _currentFrame = frameIndex;
    }

    /// <summary>
    /// Advance to the next frame in playback.
    /// </summary>
    /// <returns>True if there are more frames; false if at the end.</returns>
    public bool NextFrame(PhysicsWorld world)
    {
        int next = _currentFrame + 1;
        if (next >= _frames.Count) return false;

        ApplyFrame(world, next);
        return true;
    }

    /// <summary>
    /// Reset playback to the first frame.
    /// </summary>
    public void Rewind()
    {
        _currentFrame = 0;
    }

    /// <summary>
    /// Get a snapshot at a specific frame.
    /// </summary>
    public FrameSnapshot GetFrame(int index) => _frames[index];

    /// <summary>
    /// Clear all recorded data.
    /// </summary>
    public void Clear()
    {
        _frames.Clear();
        _currentFrame = 0;
    }

    /// <summary>
    /// Serialize all recorded frames to a binary stream.
    /// Format: [int frameCount] [frames...]
    /// Each frame: [int frame] [double time] [int bodyCount] [bodies...]
    /// Each body: [6 doubles: pos xyz, vel xyz] [3 doubles: angVel xyz]
    /// </summary>
    public void WriteTo(Stream stream)
    {
        using var writer = new BinaryWriter(stream, Encoding.UTF8, leaveOpen: true);

        writer.Write(_frames.Count);

        foreach (var frame in _frames)
        {
            writer.Write(frame.Frame);
            writer.Write(frame.Time);
            writer.Write(frame.Bodies.Length);

            foreach (var body in frame.Bodies)
            {
                writer.Write(body.Position.x);
                writer.Write(body.Position.y);
                writer.Write(body.Position.z);
                writer.Write(body.Velocity.x);
                writer.Write(body.Velocity.y);
                writer.Write(body.Velocity.z);
                writer.Write(body.AngularVelocity.x);
                writer.Write(body.AngularVelocity.y);
                writer.Write(body.AngularVelocity.z);
            }
        }
    }

    /// <summary>
    /// Deserialize recorded frames from a binary stream.
    /// </summary>
    public void ReadFrom(Stream stream)
    {
        _frames.Clear();
        _currentFrame = 0;

        using var reader = new BinaryReader(stream, Encoding.UTF8, leaveOpen: true);

        int frameCount = reader.ReadInt32();

        for (int f = 0; f < frameCount; f++)
        {
            var snapshot = new FrameSnapshot
            {
                Frame = reader.ReadInt32(),
                Time = reader.ReadDouble(),
            };

            int bodyCount = reader.ReadInt32();
            snapshot.Bodies = new BodySnapshot[bodyCount];

            for (int b = 0; b < bodyCount; b++)
            {
                snapshot.Bodies[b] = new BodySnapshot
                {
                    Position = new Vector(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble()),
                    Velocity = new Vector(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble()),
                    AngularVelocity = new Vector(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble()),
                };
            }

            _frames.Add(snapshot);
        }
    }
}
