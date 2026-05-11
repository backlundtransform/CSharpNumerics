using System;
using System.Buffers;

namespace CSharpNumerics.Engines.Game.Performance;

/// <summary>
/// Memory pooling for game engine temporary arrays.
///
/// Wraps <see cref="ArrayPool{T}"/> with a convenient Rent/Return pattern
/// to eliminate per-frame allocations in hot paths (fluid solver, collision,
/// particle system).
///
/// Usage: rent an array, use it, return it. Returned arrays are reused
/// on the next rent, avoiding GC pressure.
/// </summary>
public static class GameArrayPool
{
    private static readonly ArrayPool<double> _doublePool = ArrayPool<double>.Shared;
    private static readonly ArrayPool<int> _intPool = ArrayPool<int>.Shared;
    private static readonly ArrayPool<float> _floatPool = ArrayPool<float>.Shared;

    /// <summary>
    /// Rent a double array of at least the requested length.
    /// The returned array may be larger than requested.
    /// </summary>
    /// <param name="minimumLength">Minimum number of elements.</param>
    /// <returns>A pooled array. Must be returned via <see cref="Return(double[])"/>.</returns>
    public static double[] RentDouble(int minimumLength)
    {
        return _doublePool.Rent(minimumLength);
    }

    /// <summary>
    /// Return a rented double array to the pool.
    /// </summary>
    /// <param name="array">The array to return.</param>
    /// <param name="clearArray">Whether to zero the array before returning (security). Default false for performance.</param>
    public static void Return(double[] array, bool clearArray = false)
    {
        if (array != null)
            _doublePool.Return(array, clearArray);
    }

    /// <summary>
    /// Rent an int array.
    /// </summary>
    public static int[] RentInt(int minimumLength)
    {
        return _intPool.Rent(minimumLength);
    }

    /// <summary>
    /// Return a rented int array.
    /// </summary>
    public static void Return(int[] array, bool clearArray = false)
    {
        if (array != null)
            _intPool.Return(array, clearArray);
    }

    /// <summary>
    /// Rent a float array.
    /// </summary>
    public static float[] RentFloat(int minimumLength)
    {
        return _floatPool.Rent(minimumLength);
    }

    /// <summary>
    /// Return a rented float array.
    /// </summary>
    public static void Return(float[] array, bool clearArray = false)
    {
        if (array != null)
            _floatPool.Return(array, clearArray);
    }

    /// <summary>
    /// Rent a zeroed double array (common pattern for accumulators).
    /// </summary>
    public static double[] RentDoubleCleared(int minimumLength)
    {
        var arr = _doublePool.Rent(minimumLength);
        Array.Clear(arr, 0, minimumLength);
        return arr;
    }

    /// <summary>
    /// RAII-style handle that returns the array when disposed.
    /// Usage: using var buf = GameArrayPool.RentScope(size);
    /// </summary>
    public static PooledArray<double> RentScope(int minimumLength)
    {
        return new PooledArray<double>(_doublePool, minimumLength);
    }

    /// <summary>
    /// RAII-style handle for int arrays.
    /// </summary>
    public static PooledArray<int> RentIntScope(int minimumLength)
    {
        return new PooledArray<int>(_intPool, minimumLength);
    }
}

/// <summary>
/// Disposable handle that returns a pooled array when disposed.
/// Use with `using` statement for automatic cleanup.
/// </summary>
public struct PooledArray<T> : IDisposable
{
    private readonly ArrayPool<T> _pool;
    private T[] _array;

    /// <summary>The rented array.</summary>
    public T[] Array => _array;

    /// <summary>The usable length (requested, not actual array length).</summary>
    public int Length { get; }

    /// <summary>Access array as a Span of the requested length.</summary>
    public Span<T> Span => _array.AsSpan(0, Length);

    internal PooledArray(ArrayPool<T> pool, int length)
    {
        _pool = pool;
        Length = length;
        _array = pool.Rent(length);
    }

    /// <summary>Return the array to the pool.</summary>
    public void Dispose()
    {
        if (_array != null)
        {
            _pool.Return(_array);
            _array = null;
        }
    }
}
