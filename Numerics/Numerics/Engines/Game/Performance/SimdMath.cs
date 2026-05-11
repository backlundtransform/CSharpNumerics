using System;
#if NET8_0_OR_GREATER
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
#endif

namespace CSharpNumerics.Engines.Game.Performance;

/// <summary>
/// SIMD-accelerated bulk math operations for double arrays.
///
/// On .NET 8+ with AVX2/SSE2, uses hardware SIMD (Vector128/256).
/// Falls back to scalar loops on older runtimes (netstandard2.1).
///
/// Targets the hot inner loops of the fluid solver, particle system,
/// and collision detection where arrays of doubles are processed in bulk.
/// </summary>
public static class SimdMath
{
    /// <summary>
    /// Add two double arrays element-wise: result[i] = a[i] + b[i].
    /// </summary>
    public static void Add(ReadOnlySpan<double> a, ReadOnlySpan<double> b, Span<double> result)
    {
        int len = Math.Min(a.Length, Math.Min(b.Length, result.Length));
        int i = 0;

#if NET8_0_OR_GREATER
        if (Avx2.IsSupported)
        {
            for (; i + 4 <= len; i += 4)
            {
                var va = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vb = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(b.Slice(i)));
                var vr = Avx.Add(va, vb);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
        else if (Sse2.IsSupported)
        {
            for (; i + 2 <= len; i += 2)
            {
                var va = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vb = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(b.Slice(i)));
                var vr = Sse2.Add(va, vb);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
#endif
        for (; i < len; i++)
            result[i] = a[i] + b[i];
    }

    /// <summary>
    /// Multiply all elements by a scalar: result[i] = scalar * a[i].
    /// </summary>
    public static void Scale(ReadOnlySpan<double> a, double scalar, Span<double> result)
    {
        int len = Math.Min(a.Length, result.Length);
        int i = 0;

#if NET8_0_OR_GREATER
        if (Avx2.IsSupported)
        {
            var vs = Vector256.Create(scalar);
            for (; i + 4 <= len; i += 4)
            {
                var va = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vr = Avx.Multiply(va, vs);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
        else if (Sse2.IsSupported)
        {
            var vs = Vector128.Create(scalar);
            for (; i + 2 <= len; i += 2)
            {
                var va = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vr = Sse2.Multiply(va, vs);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
#endif
        for (; i < len; i++)
            result[i] = scalar * a[i];
    }

    /// <summary>
    /// Multiply and add: result[i] += scalar * a[i] (fused multiply-add pattern).
    /// </summary>
    public static void ScaleAdd(ReadOnlySpan<double> a, double scalar, Span<double> result)
    {
        int len = Math.Min(a.Length, result.Length);
        int i = 0;

#if NET8_0_OR_GREATER
        if (Avx2.IsSupported)
        {
            var vs = Vector256.Create(scalar);
            for (; i + 4 <= len; i += 4)
            {
                var va = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vr = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
                var prod = Avx.Multiply(va, vs);
                vr = Avx.Add(vr, prod);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
        else if (Sse2.IsSupported)
        {
            var vs = Vector128.Create(scalar);
            for (; i + 2 <= len; i += 2)
            {
                var va = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vr = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
                var prod = Sse2.Multiply(va, vs);
                vr = Sse2.Add(vr, prod);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
#endif
        for (; i < len; i++)
            result[i] += scalar * a[i];
    }

    /// <summary>
    /// Compute dot product of two double arrays.
    /// </summary>
    public static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        int len = Math.Min(a.Length, b.Length);
        double sum = 0;
        int i = 0;

#if NET8_0_OR_GREATER
        if (Avx2.IsSupported && len >= 4)
        {
            var acc = Vector256<double>.Zero;
            for (; i + 4 <= len; i += 4)
            {
                var va = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vb = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(b.Slice(i)));
                acc = Avx.Add(acc, Avx.Multiply(va, vb));
            }
            sum = acc[0] + acc[1] + acc[2] + acc[3];
        }
        else if (Sse2.IsSupported && len >= 2)
        {
            var acc = Vector128<double>.Zero;
            for (; i + 2 <= len; i += 2)
            {
                var va = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vb = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(b.Slice(i)));
                acc = Sse2.Add(acc, Sse2.Multiply(va, vb));
            }
            sum = acc[0] + acc[1];
        }
#endif
        for (; i < len; i++)
            sum += a[i] * b[i];
        return sum;
    }

    /// <summary>
    /// Decay all elements: a[i] *= factor.
    /// Used for density decay in fluid solvers.
    /// </summary>
    public static void Decay(Span<double> a, double factor)
    {
        int len = a.Length;
        int i = 0;

#if NET8_0_OR_GREATER
        if (Avx2.IsSupported)
        {
            var vf = Vector256.Create(factor);
            for (; i + 4 <= len; i += 4)
            {
                var va = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vr = Avx.Multiply(va, vf);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
            }
        }
        else if (Sse2.IsSupported)
        {
            var vf = Vector128.Create(factor);
            for (; i + 2 <= len; i += 2)
            {
                var va = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vr = Sse2.Multiply(va, vf);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
            }
        }
#endif
        for (; i < len; i++)
            a[i] *= factor;
    }

    /// <summary>
    /// Subtract two arrays: result[i] = a[i] - b[i].
    /// </summary>
    public static void Subtract(ReadOnlySpan<double> a, ReadOnlySpan<double> b, Span<double> result)
    {
        int len = Math.Min(a.Length, Math.Min(b.Length, result.Length));
        int i = 0;

#if NET8_0_OR_GREATER
        if (Avx2.IsSupported)
        {
            for (; i + 4 <= len; i += 4)
            {
                var va = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vb = Vector256.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(b.Slice(i)));
                var vr = Avx.Subtract(va, vb);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
        else if (Sse2.IsSupported)
        {
            for (; i + 2 <= len; i += 2)
            {
                var va = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(a.Slice(i)));
                var vb = Vector128.LoadUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(b.Slice(i)));
                var vr = Sse2.Subtract(va, vb);
                vr.StoreUnsafe(ref System.Runtime.InteropServices.MemoryMarshal.GetReference(result.Slice(i)));
            }
        }
#endif
        for (; i < len; i++)
            result[i] = a[i] - b[i];
    }

    /// <summary>
    /// Copy array: dst[i] = src[i]. Uses optimized memory copy.
    /// </summary>
    public static void Copy(ReadOnlySpan<double> src, Span<double> dst)
    {
        src.Slice(0, Math.Min(src.Length, dst.Length)).CopyTo(dst);
    }

    /// <summary>
    /// Fill array with a constant value.
    /// </summary>
    public static void Fill(Span<double> a, double value)
    {
        a.Fill(value);
    }

    /// <summary>
    /// Zero out an array.
    /// </summary>
    public static void Clear(Span<double> a)
    {
        a.Clear();
    }
}
