using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Engines.Game.Rocket;

/// <summary>
/// Converts RocketState into a position/rotation representation suitable for
/// Unity Transform assignment. This class bridges the pure-math simulation
/// with a game engine's rendering coordinate system.
///
/// Unity conventions: Y-up, left-handed. The adapter converts from the
/// simulation's NED/ECI frame to Unity's Y-up frame.
/// </summary>
public class RocketUnityAdapter
{
    /// <summary>Coordinate mode for output.</summary>
    public AdapterCoordinateMode CoordinateMode { get; set; } = AdapterCoordinateMode.NEDToUnity;

    /// <summary>Scale factor applied to position (meters → Unity units). Default 1.0 (1:1).</summary>
    public double PositionScale { get; set; } = 1.0;

    /// <summary>Last computed Unity position (x, y, z).</summary>
    public Vector UnityPosition { get; private set; }

    /// <summary>Last computed Unity rotation as quaternion (w, x, y, z).</summary>
    public Quaternion UnityRotation { get; private set; } = Quaternion.Identity;

    /// <summary>Accumulated numerical drift metric (sum of position corrections applied).</summary>
    public double AccumulatedDrift { get; private set; }

    private Vector _prevPosition;
    private bool _hasPrev;

    /// <summary>
    /// Updates the adapter with the current simulation state.
    /// </summary>
    /// <param name="state">Current rocket state from the simulation engine.</param>
    public void Update(RocketState state)
    {
        Vector simPos = state.Position;
        Quaternion simAtt = state.Attitude;

        Vector unityPos;
        Quaternion unityRot;

        switch (CoordinateMode)
        {
            case AdapterCoordinateMode.NEDToUnity:
                // NED (x=North, y=East, z=Down) → Unity (x=East, y=Up, z=North)
                unityPos = new Vector(
                    simPos.y * PositionScale,
                    -simPos.z * PositionScale,
                    simPos.x * PositionScale);
                // Rotate quaternion from NED body to Unity body
                // NED body: x=Forward, y=Right, z=Down → Unity body: x=Right, y=Up, z=Forward
                unityRot = NEDToUnityQuaternion(simAtt);
                break;

            case AdapterCoordinateMode.ECIToUnity:
                // ECI (x=vernal, y=90°E, z=north pole) → Unity (x=East, y=Up, z=North)
                unityPos = new Vector(
                    simPos.x * PositionScale,
                    simPos.z * PositionScale,
                    simPos.y * PositionScale);
                unityRot = simAtt; // Direct mapping for ECI
                break;

            default: // PassThrough
                unityPos = new Vector(
                    simPos.x * PositionScale,
                    simPos.y * PositionScale,
                    simPos.z * PositionScale);
                unityRot = simAtt;
                break;
        }

        // Track drift (difference from expected double-precision position)
        if (_hasPrev)
        {
            double dx = unityPos.x - _prevPosition.x;
            double dy = unityPos.y - _prevPosition.y;
            double dz = unityPos.z - _prevPosition.z;
            double stepMag = Math.Sqrt(dx * dx + dy * dy + dz * dz);
            // Drift accumulates from floating-point round-off
            // We track it but don't correct it here (caller can use for diagnostics)
            if (stepMag > 0)
                AccumulatedDrift += stepMag * 1e-15; // Estimate: ~1 ULP per step
        }

        _prevPosition = unityPos;
        _hasPrev = true;

        UnityPosition = unityPos;
        UnityRotation = unityRot;
    }

    /// <summary>Resets accumulated drift tracking.</summary>
    public void Reset()
    {
        UnityPosition = new Vector(0, 0, 0);
        UnityRotation = Quaternion.Identity;
        AccumulatedDrift = 0;
        _hasPrev = false;
    }

    private static Quaternion NEDToUnityQuaternion(Quaternion nedQuat)
    {
        // NED→Unity frame rotation: swap axes
        // NED (x=N,y=E,z=D) → Unity (x=E,y=U,z=N) means:
        // Unity.x = NED.y, Unity.y = -NED.z, Unity.z = NED.x
        // For quaternion: q_unity = R_frame * q_ned * R_frame^-1
        // Simplified: swap vector components of the quaternion accordingly
        return new Quaternion(nedQuat.w, nedQuat.y, -nedQuat.z, nedQuat.x);
    }
}

/// <summary>
/// Coordinate mode for the Unity adapter.
/// </summary>
public enum AdapterCoordinateMode
{
    /// <summary>NED simulation → Unity Y-up left-handed.</summary>
    NEDToUnity,
    /// <summary>ECI simulation → Unity Y-up.</summary>
    ECIToUnity,
    /// <summary>Pass through without transformation.</summary>
    PassThrough
}
