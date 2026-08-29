using System;
using CSharpNumerics.Numerics.Objects;

namespace CSharpNumerics.Numerics.LinearAlgebra.Decompositions;

/// <summary>
/// Eigenvalue decomposition A·V = V·D.
/// For symmetric matrices: Householder tridiagonalization followed by the implicit QL algorithm —
/// all eigenvalues are real, returned in ascending order with orthonormal eigenvectors.
/// For non-symmetric matrices: Hessenberg reduction followed by the shifted QR algorithm —
/// eigenvalues may be complex conjugate pairs (real parts in <see cref="RealEigenvalues"/>,
/// imaginary parts in <see cref="ImaginaryEigenvalues"/>) and D is block diagonal.
/// Based on the EISPACK/JAMA reference implementation.
/// </summary>
public sealed class EigenDecomposition
{
    private const double SymmetryTolerance = 1e-10;

    private readonly int n;
    private readonly double[] d;
    private readonly double[] e;
    private readonly double[,] v;
    private double[,] hess;
    private double[] ort;
    private double cdivr, cdivi;

    public EigenDecomposition(Matrix matrix)
    {
        if (matrix.rowLength != matrix.columnLength)
        {
            throw new Exception("Is not a NxN matrix");
        }

        n = matrix.rowLength;
        d = new double[n];
        e = new double[n];
        v = new double[n, n];
        IsSymmetric = CheckSymmetric(matrix);

        if (IsSymmetric)
        {
            for (var i = 0; i < n; i++)
            {
                for (var j = 0; j < n; j++)
                {
                    v[i, j] = matrix.values[i, j];
                }
            }

            Tridiagonalize();
            TridiagonalQl();
        }
        else
        {
            hess = (double[,])matrix.values.Clone();
            ort = new double[n];

            ReduceToHessenberg();
            HessenbergQr();
        }
    }

    /// <summary>
    /// True if the matrix was detected as symmetric (all eigenvalues real, eigenvectors orthonormal).
    /// </summary>
    public bool IsSymmetric { get; }

    /// <summary>
    /// Real parts of the eigenvalues. For symmetric matrices these are all eigenvalues, in ascending order.
    /// </summary>
    public double[] RealEigenvalues => (double[])d.Clone();

    /// <summary>
    /// Imaginary parts of the eigenvalues (all zero for symmetric matrices).
    /// </summary>
    public double[] ImaginaryEigenvalues => (double[])e.Clone();

    /// <summary>
    /// The eigenvector matrix V (eigenvector k in column k, matching eigenvalue k).
    /// </summary>
    public Matrix EigenVectors => new Matrix((double[,])v.Clone());

    /// <summary>
    /// The block diagonal eigenvalue matrix D such that A·V = V·D.
    /// Complex conjugate pairs appear as 2×2 blocks [[re, im], [-im, re]].
    /// </summary>
    public Matrix DiagonalMatrix
    {
        get
        {
            var diagonal = new double[n, n];
            for (var i = 0; i < n; i++)
            {
                diagonal[i, i] = d[i];
                if (e[i] > 0)
                {
                    diagonal[i, i + 1] = e[i];
                }
                else if (e[i] < 0)
                {
                    diagonal[i, i - 1] = e[i];
                }
            }
            return new Matrix(diagonal);
        }
    }

    private bool CheckSymmetric(Matrix matrix)
    {
        for (var i = 0; i < n; i++)
        {
            for (var j = i + 1; j < n; j++)
            {
                var a = matrix.values[i, j];
                var b = matrix.values[j, i];
                var scale = Math.Max(1.0, Math.Max(Math.Abs(a), Math.Abs(b)));

                if (Math.Abs(a - b) > SymmetryTolerance * scale)
                {
                    return false;
                }
            }
        }
        return true;
    }

    // Householder reduction of a symmetric matrix to tridiagonal form (EISPACK tred2).
    private void Tridiagonalize()
    {
        for (var j = 0; j < n; j++)
        {
            d[j] = v[n - 1, j];
        }

        for (var i = n - 1; i > 0; i--)
        {
            var scale = 0.0;
            var h = 0.0;
            for (var k = 0; k < i; k++)
            {
                scale += Math.Abs(d[k]);
            }

            if (scale == 0.0)
            {
                e[i] = d[i - 1];
                for (var j = 0; j < i; j++)
                {
                    d[j] = v[i - 1, j];
                    v[i, j] = 0.0;
                    v[j, i] = 0.0;
                }
            }
            else
            {
                for (var k = 0; k < i; k++)
                {
                    d[k] /= scale;
                    h += d[k] * d[k];
                }

                var f = d[i - 1];
                var g = Math.Sqrt(h);
                if (f > 0)
                {
                    g = -g;
                }

                e[i] = scale * g;
                h -= f * g;
                d[i - 1] = f - g;
                for (var j = 0; j < i; j++)
                {
                    e[j] = 0.0;
                }

                for (var j = 0; j < i; j++)
                {
                    f = d[j];
                    v[j, i] = f;
                    g = e[j] + v[j, j] * f;
                    for (var k = j + 1; k <= i - 1; k++)
                    {
                        g += v[k, j] * d[k];
                        e[k] += v[k, j] * f;
                    }
                    e[j] = g;
                }

                f = 0.0;
                for (var j = 0; j < i; j++)
                {
                    e[j] /= h;
                    f += e[j] * d[j];
                }

                var hh = f / (h + h);
                for (var j = 0; j < i; j++)
                {
                    e[j] -= hh * d[j];
                }

                for (var j = 0; j < i; j++)
                {
                    f = d[j];
                    g = e[j];
                    for (var k = j; k <= i - 1; k++)
                    {
                        v[k, j] -= f * e[k] + g * d[k];
                    }
                    d[j] = v[i - 1, j];
                    v[i, j] = 0.0;
                }
            }

            d[i] = h;
        }

        for (var i = 0; i < n - 1; i++)
        {
            v[n - 1, i] = v[i, i];
            v[i, i] = 1.0;
            var h = d[i + 1];

            if (h != 0.0)
            {
                for (var k = 0; k <= i; k++)
                {
                    d[k] = v[k, i + 1] / h;
                }
                for (var j = 0; j <= i; j++)
                {
                    var g = 0.0;
                    for (var k = 0; k <= i; k++)
                    {
                        g += v[k, i + 1] * v[k, j];
                    }
                    for (var k = 0; k <= i; k++)
                    {
                        v[k, j] -= g * d[k];
                    }
                }
            }

            for (var k = 0; k <= i; k++)
            {
                v[k, i + 1] = 0.0;
            }
        }

        for (var j = 0; j < n; j++)
        {
            d[j] = v[n - 1, j];
            v[n - 1, j] = 0.0;
        }
        v[n - 1, n - 1] = 1.0;
        e[0] = 0.0;
    }

    // Implicit QL algorithm with shifts for a symmetric tridiagonal matrix (EISPACK tql2).
    // Eigenvalues end up sorted ascending with matching eigenvector columns.
    private void TridiagonalQl()
    {
        for (var i = 1; i < n; i++)
        {
            e[i - 1] = e[i];
        }
        e[n - 1] = 0.0;

        var f = 0.0;
        var tst1 = 0.0;
        var eps = Math.Pow(2.0, -52.0);

        for (var l = 0; l < n; l++)
        {
            tst1 = Math.Max(tst1, Math.Abs(d[l]) + Math.Abs(e[l]));
            var m = l;
            while (m < n)
            {
                if (Math.Abs(e[m]) <= eps * tst1)
                {
                    break;
                }
                m++;
            }

            if (m > l)
            {
                do
                {
                    var g = d[l];
                    var p = (d[l + 1] - g) / (2.0 * e[l]);
                    var r = Hypot(p, 1.0);
                    if (p < 0)
                    {
                        r = -r;
                    }

                    d[l] = e[l] / (p + r);
                    d[l + 1] = e[l] * (p + r);
                    var dl1 = d[l + 1];
                    var h = g - d[l];
                    for (var i = l + 2; i < n; i++)
                    {
                        d[i] -= h;
                    }
                    f += h;

                    p = d[m];
                    var c = 1.0;
                    var c2 = c;
                    var c3 = c;
                    var el1 = e[l + 1];
                    var s = 0.0;
                    var s2 = 0.0;

                    for (var i = m - 1; i >= l; i--)
                    {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = Hypot(p, e[i]);
                        e[i + 1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i + 1] = h + s * (c * g + s * d[i]);

                        for (var k = 0; k < n; k++)
                        {
                            h = v[k, i + 1];
                            v[k, i + 1] = s * v[k, i] + c * h;
                            v[k, i] = c * v[k, i] - s * h;
                        }
                    }

                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;
                } while (Math.Abs(e[l]) > eps * tst1);
            }

            d[l] += f;
            e[l] = 0.0;
        }

        for (var i = 0; i < n - 1; i++)
        {
            var k = i;
            var p = d[i];
            for (var j = i + 1; j < n; j++)
            {
                if (d[j] < p)
                {
                    k = j;
                    p = d[j];
                }
            }

            if (k != i)
            {
                d[k] = d[i];
                d[i] = p;
                for (var j = 0; j < n; j++)
                {
                    (v[j, i], v[j, k]) = (v[j, k], v[j, i]);
                }
            }
        }
    }

    // Householder reduction of a general matrix to Hessenberg form (EISPACK orthes).
    private void ReduceToHessenberg()
    {
        var low = 0;
        var high = n - 1;

        for (var m = low + 1; m <= high - 1; m++)
        {
            var scale = 0.0;
            for (var i = m; i <= high; i++)
            {
                scale += Math.Abs(hess[i, m - 1]);
            }

            if (scale != 0.0)
            {
                var h = 0.0;
                for (var i = high; i >= m; i--)
                {
                    ort[i] = hess[i, m - 1] / scale;
                    h += ort[i] * ort[i];
                }

                var g = Math.Sqrt(h);
                if (ort[m] > 0)
                {
                    g = -g;
                }
                h -= ort[m] * g;
                ort[m] -= g;

                for (var j = m; j < n; j++)
                {
                    var f = 0.0;
                    for (var i = high; i >= m; i--)
                    {
                        f += ort[i] * hess[i, j];
                    }
                    f /= h;
                    for (var i = m; i <= high; i++)
                    {
                        hess[i, j] -= f * ort[i];
                    }
                }

                for (var i = 0; i <= high; i++)
                {
                    var f = 0.0;
                    for (var j = high; j >= m; j--)
                    {
                        f += ort[j] * hess[i, j];
                    }
                    f /= h;
                    for (var j = m; j <= high; j++)
                    {
                        hess[i, j] -= f * ort[j];
                    }
                }

                ort[m] *= scale;
                hess[m, m - 1] = scale * g;
            }
        }

        for (var i = 0; i < n; i++)
        {
            for (var j = 0; j < n; j++)
            {
                v[i, j] = i == j ? 1.0 : 0.0;
            }
        }

        for (var m = high - 1; m >= low + 1; m--)
        {
            if (hess[m, m - 1] != 0.0)
            {
                for (var i = m + 1; i <= high; i++)
                {
                    ort[i] = hess[i, m - 1];
                }

                for (var j = m; j <= high; j++)
                {
                    var g = 0.0;
                    for (var i = m; i <= high; i++)
                    {
                        g += ort[i] * v[i, j];
                    }
                    g = (g / ort[m]) / hess[m, m - 1];
                    for (var i = m; i <= high; i++)
                    {
                        v[i, j] += g * ort[i];
                    }
                }
            }
        }
    }

    private void ComplexDivide(double xr, double xi, double yr, double yi)
    {
        if (Math.Abs(yr) > Math.Abs(yi))
        {
            var r = yi / yr;
            var dd = yr + r * yi;
            cdivr = (xr + r * xi) / dd;
            cdivi = (xi - r * xr) / dd;
        }
        else
        {
            var r = yr / yi;
            var dd = yi + r * yr;
            cdivr = (r * xr + xi) / dd;
            cdivi = (r * xi - xr) / dd;
        }
    }

    // Shifted QR iteration on a Hessenberg matrix with eigenvector back-substitution (EISPACK hqr2).
    private void HessenbergQr()
    {
        var nn = n;
        var current = nn - 1;
        var low = 0;
        var high = nn - 1;
        var eps = Math.Pow(2.0, -52.0);
        var exshift = 0.0;
        double p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;

        var norm = 0.0;
        for (var i = 0; i < nn; i++)
        {
            if (i < low || i > high)
            {
                d[i] = hess[i, i];
                e[i] = 0.0;
            }
            for (var j = Math.Max(i - 1, 0); j < nn; j++)
            {
                norm += Math.Abs(hess[i, j]);
            }
        }

        var iter = 0;
        while (current >= low)
        {
            var l = current;
            while (l > low)
            {
                s = Math.Abs(hess[l - 1, l - 1]) + Math.Abs(hess[l, l]);
                if (s == 0.0)
                {
                    s = norm;
                }
                if (Math.Abs(hess[l, l - 1]) < eps * s)
                {
                    break;
                }
                l--;
            }

            if (l == current)
            {
                hess[current, current] += exshift;
                d[current] = hess[current, current];
                e[current] = 0.0;
                current--;
                iter = 0;
            }
            else if (l == current - 1)
            {
                w = hess[current, current - 1] * hess[current - 1, current];
                p = (hess[current - 1, current - 1] - hess[current, current]) / 2.0;
                q = p * p + w;
                z = Math.Sqrt(Math.Abs(q));
                hess[current, current] += exshift;
                hess[current - 1, current - 1] += exshift;
                x = hess[current, current];

                if (q >= 0)
                {
                    z = p >= 0 ? p + z : p - z;
                    d[current - 1] = x + z;
                    d[current] = d[current - 1];
                    if (z != 0.0)
                    {
                        d[current] = x - w / z;
                    }
                    e[current - 1] = 0.0;
                    e[current] = 0.0;
                    x = hess[current, current - 1];
                    s = Math.Abs(x) + Math.Abs(z);
                    p = x / s;
                    q = z / s;
                    r = Math.Sqrt(p * p + q * q);
                    p /= r;
                    q /= r;

                    for (var j = current - 1; j < nn; j++)
                    {
                        z = hess[current - 1, j];
                        hess[current - 1, j] = q * z + p * hess[current, j];
                        hess[current, j] = q * hess[current, j] - p * z;
                    }

                    for (var i = 0; i <= current; i++)
                    {
                        z = hess[i, current - 1];
                        hess[i, current - 1] = q * z + p * hess[i, current];
                        hess[i, current] = q * hess[i, current] - p * z;
                    }

                    for (var i = low; i <= high; i++)
                    {
                        z = v[i, current - 1];
                        v[i, current - 1] = q * z + p * v[i, current];
                        v[i, current] = q * v[i, current] - p * z;
                    }
                }
                else
                {
                    d[current - 1] = x + p;
                    d[current] = x + p;
                    e[current - 1] = z;
                    e[current] = -z;
                }

                current -= 2;
                iter = 0;
            }
            else
            {
                x = hess[current, current];
                y = 0.0;
                w = 0.0;
                if (l < current)
                {
                    y = hess[current - 1, current - 1];
                    w = hess[current, current - 1] * hess[current - 1, current];
                }

                if (iter == 10 || iter == 20)
                {
                    exshift += x;
                    for (var i = low; i <= current; i++)
                    {
                        hess[i, i] -= x;
                    }
                    s = Math.Abs(hess[current, current - 1]) + Math.Abs(hess[current - 1, current - 2]);
                    x = y = 0.75 * s;
                    w = -0.4375 * s * s;
                }

                if (iter == 30)
                {
                    s = (y - x) / 2.0;
                    s = s * s + w;
                    if (s > 0)
                    {
                        s = Math.Sqrt(s);
                        if (y < x)
                        {
                            s = -s;
                        }
                        s = x - w / ((y - x) / 2.0 + s);
                        for (var i = low; i <= current; i++)
                        {
                            hess[i, i] -= s;
                        }
                        exshift += s;
                        x = y = w = 0.964;
                    }
                }

                iter++;
                if (iter > 250)
                {
                    throw new Exception("Eigenvalue iteration did not converge");
                }

                var m = current - 2;
                while (m >= l)
                {
                    z = hess[m, m];
                    r = x - z;
                    s = y - z;
                    p = (r * s - w) / hess[m + 1, m] + hess[m, m + 1];
                    q = hess[m + 1, m + 1] - z - r - s;
                    r = hess[m + 2, m + 1];
                    s = Math.Abs(p) + Math.Abs(q) + Math.Abs(r);
                    p /= s;
                    q /= s;
                    r /= s;
                    if (m == l)
                    {
                        break;
                    }
                    if (Math.Abs(hess[m, m - 1]) * (Math.Abs(q) + Math.Abs(r)) <
                        eps * (Math.Abs(p) * (Math.Abs(hess[m - 1, m - 1]) + Math.Abs(z) + Math.Abs(hess[m + 1, m + 1]))))
                    {
                        break;
                    }
                    m--;
                }

                for (var i = m + 2; i <= current; i++)
                {
                    hess[i, i - 2] = 0.0;
                    if (i > m + 2)
                    {
                        hess[i, i - 3] = 0.0;
                    }
                }

                for (var k = m; k <= current - 1; k++)
                {
                    var notLast = k != current - 1;
                    if (k != m)
                    {
                        p = hess[k, k - 1];
                        q = hess[k + 1, k - 1];
                        r = notLast ? hess[k + 2, k - 1] : 0.0;
                        x = Math.Abs(p) + Math.Abs(q) + Math.Abs(r);
                        if (x == 0.0)
                        {
                            continue;
                        }
                        p /= x;
                        q /= x;
                        r /= x;
                    }

                    s = Math.Sqrt(p * p + q * q + r * r);
                    if (p < 0)
                    {
                        s = -s;
                    }

                    if (s != 0)
                    {
                        if (k != m)
                        {
                            hess[k, k - 1] = -s * x;
                        }
                        else if (l != m)
                        {
                            hess[k, k - 1] = -hess[k, k - 1];
                        }

                        p += s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q /= p;
                        r /= p;

                        for (var j = k; j < nn; j++)
                        {
                            p = hess[k, j] + q * hess[k + 1, j];
                            if (notLast)
                            {
                                p += r * hess[k + 2, j];
                                hess[k + 2, j] -= p * z;
                            }
                            hess[k, j] -= p * x;
                            hess[k + 1, j] -= p * y;
                        }

                        for (var i = 0; i <= Math.Min(current, k + 3); i++)
                        {
                            p = x * hess[i, k] + y * hess[i, k + 1];
                            if (notLast)
                            {
                                p += z * hess[i, k + 2];
                                hess[i, k + 2] -= p * r;
                            }
                            hess[i, k] -= p;
                            hess[i, k + 1] -= p * q;
                        }

                        for (var i = low; i <= high; i++)
                        {
                            p = x * v[i, k] + y * v[i, k + 1];
                            if (notLast)
                            {
                                p += z * v[i, k + 2];
                                v[i, k + 2] -= p * r;
                            }
                            v[i, k] -= p;
                            v[i, k + 1] -= p * q;
                        }
                    }
                }
            }
        }

        if (norm == 0.0)
        {
            return;
        }

        for (current = nn - 1; current >= 0; current--)
        {
            p = d[current];
            q = e[current];

            if (q == 0)
            {
                var l = current;
                hess[current, current] = 1.0;
                for (var i = current - 1; i >= 0; i--)
                {
                    w = hess[i, i] - p;
                    r = 0.0;
                    for (var j = l; j <= current; j++)
                    {
                        r += hess[i, j] * hess[j, current];
                    }

                    if (e[i] < 0.0)
                    {
                        z = w;
                        s = r;
                    }
                    else
                    {
                        l = i;
                        if (e[i] == 0.0)
                        {
                            hess[i, current] = w != 0.0 ? -r / w : -r / (eps * norm);
                        }
                        else
                        {
                            x = hess[i, i + 1];
                            y = hess[i + 1, i];
                            q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                            t = (x * s - z * r) / q;
                            hess[i, current] = t;
                            hess[i + 1, current] = Math.Abs(x) > Math.Abs(z)
                                ? (-r - w * t) / x
                                : (-s - y * t) / z;
                        }

                        t = Math.Abs(hess[i, current]);
                        if (eps * t * t > 1)
                        {
                            for (var j = i; j <= current; j++)
                            {
                                hess[j, current] /= t;
                            }
                        }
                    }
                }
            }
            else if (q < 0)
            {
                var l = current - 1;

                if (Math.Abs(hess[current, current - 1]) > Math.Abs(hess[current - 1, current]))
                {
                    hess[current - 1, current - 1] = q / hess[current, current - 1];
                    hess[current - 1, current] = -(hess[current, current] - p) / hess[current, current - 1];
                }
                else
                {
                    ComplexDivide(0.0, -hess[current - 1, current], hess[current - 1, current - 1] - p, q);
                    hess[current - 1, current - 1] = cdivr;
                    hess[current - 1, current] = cdivi;
                }

                hess[current, current - 1] = 0.0;
                hess[current, current] = 1.0;
                for (var i = current - 2; i >= 0; i--)
                {
                    var ra = 0.0;
                    var sa = 0.0;
                    for (var j = l; j <= current; j++)
                    {
                        ra += hess[i, j] * hess[j, current - 1];
                        sa += hess[i, j] * hess[j, current];
                    }
                    w = hess[i, i] - p;

                    if (e[i] < 0.0)
                    {
                        z = w;
                        r = ra;
                        s = sa;
                    }
                    else
                    {
                        l = i;
                        if (e[i] == 0)
                        {
                            ComplexDivide(-ra, -sa, w, q);
                            hess[i, current - 1] = cdivr;
                            hess[i, current] = cdivi;
                        }
                        else
                        {
                            x = hess[i, i + 1];
                            y = hess[i + 1, i];
                            var vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                            var vi = (d[i] - p) * 2.0 * q;
                            if (vr == 0.0 && vi == 0.0)
                            {
                                vr = eps * norm * (Math.Abs(w) + Math.Abs(q) + Math.Abs(x) + Math.Abs(y) + Math.Abs(z));
                            }
                            ComplexDivide(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                            hess[i, current - 1] = cdivr;
                            hess[i, current] = cdivi;

                            if (Math.Abs(x) > Math.Abs(z) + Math.Abs(q))
                            {
                                hess[i + 1, current - 1] = (-ra - w * hess[i, current - 1] + q * hess[i, current]) / x;
                                hess[i + 1, current] = (-sa - w * hess[i, current] - q * hess[i, current - 1]) / x;
                            }
                            else
                            {
                                ComplexDivide(-r - y * hess[i, current - 1], -s - y * hess[i, current], z, q);
                                hess[i + 1, current - 1] = cdivr;
                                hess[i + 1, current] = cdivi;
                            }
                        }

                        t = Math.Max(Math.Abs(hess[i, current - 1]), Math.Abs(hess[i, current]));
                        if (eps * t * t > 1)
                        {
                            for (var j = i; j <= current; j++)
                            {
                                hess[j, current - 1] /= t;
                                hess[j, current] /= t;
                            }
                        }
                    }
                }
            }
        }

        for (var j = nn - 1; j >= low; j--)
        {
            for (var i = low; i <= high; i++)
            {
                z = 0.0;
                for (var k = low; k <= Math.Min(j, high); k++)
                {
                    z += v[i, k] * hess[k, j];
                }
                v[i, j] = z;
            }
        }
    }

    private static double Hypot(double a, double b)
    {
        if (Math.Abs(a) > Math.Abs(b))
        {
            var r = b / a;
            return Math.Abs(a) * Math.Sqrt(1 + r * r);
        }
        if (b != 0)
        {
            var r = a / b;
            return Math.Abs(b) * Math.Sqrt(1 + r * r);
        }
        return 0.0;
    }
}
