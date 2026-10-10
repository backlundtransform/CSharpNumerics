using System;

namespace CSharpNumerics.Numerics.RootFinding;

/// <summary>
/// Scalar root finders for f(x) = 0.
/// </summary>
/// <remarks>
/// <para>
/// Every method reports whether it converged rather than silently returning a bad value —
/// see <see cref="RootResult"/>.
/// </para>
/// <para>
/// Choosing a method: <see cref="Brent"/> is the default choice when a bracket is known.
/// It cannot fail to converge and is nearly as fast as Newton's method. <see cref="Bisection"/>
/// is slower but the most robust. <see cref="Secant"/> and <see cref="Newton"/> need only a
/// starting point rather than a bracket, but may diverge.
/// </para>
/// </remarks>
public static class RootFinder
{
    /// <summary>Default convergence tolerance.</summary>
    public const double DefaultTolerance = 1e-12;

    /// <summary>Default iteration budget.</summary>
    public const int DefaultMaxIterations = 100;

    /// <summary>
    /// Machine epsilon for <see cref="double"/> (2⁻⁵²).
    /// Note this is not <see cref="double.Epsilon"/>, which is the smallest subnormal value.
    /// </summary>
    private const double MachineEpsilon = 2.220446049250313e-16;

    /// <summary>
    /// Step used by the default finite-difference derivative in <see cref="Newton(Func{double, double}, double, double, int)"/>.
    /// Matches <see cref="DerivativeExtensions.Derivate(Func{double, double}, double, int)"/> at order 1.
    /// </summary>
    private const double DerivativeStep = 10.0 * DerivativeExtensions.h;

    /// <summary>
    /// Bisection: repeatedly halves a bracketing interval. Converges for any continuous
    /// function that changes sign across the bracket, but only linearly — roughly one
    /// extra correct bit per iteration.
    /// </summary>
    /// <param name="function">The function whose root is sought.</param>
    /// <param name="lower">Lower end of a bracketing interval.</param>
    /// <param name="upper">Upper end of a bracketing interval.</param>
    /// <param name="tolerance">Convergence tolerance on the half-width of the interval.</param>
    /// <param name="maxIterations">Maximum number of iterations.</param>
    /// <exception cref="ArgumentNullException"><paramref name="function"/> is null.</exception>
    /// <exception cref="ArgumentException">The interval is degenerate, or does not bracket a sign change.</exception>
    public static RootResult Bisection(
        Func<double, double> function,
        double lower,
        double upper,
        double tolerance = DefaultTolerance,
        int maxIterations = DefaultMaxIterations)
    {
        ValidateBracket(function, ref lower, ref upper, tolerance, maxIterations, out var fLower, out var fUpper);

        if (fLower == 0.0)
        {
            return new RootResult(lower, true, 0, 0.0, nameof(Bisection));
        }

        if (fUpper == 0.0)
        {
            return new RootResult(upper, true, 0, 0.0, nameof(Bisection));
        }

        var a = lower;
        var b = upper;
        var fa = fLower;
        var midpoint = a;
        var fMid = fa;

        for (var iteration = 1; iteration <= maxIterations; iteration++)
        {
            midpoint = 0.5 * (a + b);
            fMid = function(midpoint);

            if (fMid == 0.0 || 0.5 * (b - a) <= tolerance)
            {
                return new RootResult(midpoint, true, iteration, Math.Abs(fMid), nameof(Bisection));
            }

            if (Math.Sign(fMid) == Math.Sign(fa))
            {
                a = midpoint;
                fa = fMid;
            }
            else
            {
                b = midpoint;
            }
        }

        return new RootResult(midpoint, false, maxIterations, Math.Abs(fMid), nameof(Bisection));
    }

    /// <summary>
    /// Secant method: Newton's method with the derivative replaced by a finite difference
    /// over the two previous iterates. Needs no derivative and no bracket, but may diverge.
    /// </summary>
    /// <param name="function">The function whose root is sought.</param>
    /// <param name="first">First starting point.</param>
    /// <param name="second">Second starting point, distinct from the first.</param>
    /// <param name="tolerance">Convergence tolerance on the step size and residual.</param>
    /// <param name="maxIterations">Maximum number of iterations.</param>
    /// <exception cref="ArgumentNullException"><paramref name="function"/> is null.</exception>
    /// <exception cref="ArgumentException">The two starting points are equal.</exception>
    public static RootResult Secant(
        Func<double, double> function,
        double first,
        double second,
        double tolerance = DefaultTolerance,
        int maxIterations = DefaultMaxIterations)
    {
        ValidateCommon(function, tolerance, maxIterations);

        if (first == second)
        {
            throw new ArgumentException("The two starting points must be distinct.", nameof(second));
        }

        var previous = first;
        var current = second;
        var fPrevious = function(previous);
        var fCurrent = function(current);

        for (var iteration = 1; iteration <= maxIterations; iteration++)
        {
            if (Math.Abs(fCurrent) <= tolerance)
            {
                return new RootResult(current, true, iteration - 1, Math.Abs(fCurrent), nameof(Secant));
            }

            var denominator = fCurrent - fPrevious;

            // The secant through two equal function values is horizontal: no step is defined.
            if (denominator == 0.0)
            {
                return new RootResult(current, false, iteration - 1, Math.Abs(fCurrent), nameof(Secant));
            }

            var step = fCurrent * (current - previous) / denominator;
            var next = current - step;

            if (!IsFinite(next))
            {
                return new RootResult(current, false, iteration - 1, Math.Abs(fCurrent), nameof(Secant));
            }

            previous = current;
            fPrevious = fCurrent;
            current = next;
            fCurrent = function(current);

            if (Math.Abs(step) <= tolerance * (1.0 + Math.Abs(current)))
            {
                return new RootResult(current, true, iteration, Math.Abs(fCurrent), nameof(Secant));
            }
        }

        return new RootResult(current, false, maxIterations, Math.Abs(fCurrent), nameof(Secant));
    }

    /// <summary>
    /// Brent's method: combines bisection, the secant method and inverse quadratic
    /// interpolation. Keeps the root bracketed at all times, so it always converges,
    /// while achieving superlinear speed on well-behaved functions. The default choice
    /// when a bracket is available.
    /// </summary>
    /// <param name="function">The function whose root is sought.</param>
    /// <param name="lower">Lower end of a bracketing interval.</param>
    /// <param name="upper">Upper end of a bracketing interval.</param>
    /// <param name="tolerance">Convergence tolerance on the root position.</param>
    /// <param name="maxIterations">Maximum number of iterations.</param>
    /// <exception cref="ArgumentNullException"><paramref name="function"/> is null.</exception>
    /// <exception cref="ArgumentException">The interval is degenerate, or does not bracket a sign change.</exception>
    public static RootResult Brent(
        Func<double, double> function,
        double lower,
        double upper,
        double tolerance = DefaultTolerance,
        int maxIterations = DefaultMaxIterations)
    {
        ValidateBracket(function, ref lower, ref upper, tolerance, maxIterations, out var fa, out var fb);

        if (fa == 0.0)
        {
            return new RootResult(lower, true, 0, 0.0, nameof(Brent));
        }

        if (fb == 0.0)
        {
            return new RootResult(upper, true, 0, 0.0, nameof(Brent));
        }

        var a = lower;
        var b = upper;
        var c = a;
        var fc = fa;
        var d = b - a;
        var e = d;

        for (var iteration = 1; iteration <= maxIterations; iteration++)
        {
            // Keep c on the opposite side of the root from b.
            if (Math.Sign(fb) == Math.Sign(fc))
            {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }

            // Keep b as the best estimate so far.
            if (Math.Abs(fc) < Math.Abs(fb))
            {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            var tolerance1 = 2.0 * MachineEpsilon * Math.Abs(b) + 0.5 * tolerance;
            var bisectionStep = 0.5 * (c - b);

            if (Math.Abs(bisectionStep) <= tolerance1 || fb == 0.0)
            {
                return new RootResult(b, true, iteration, Math.Abs(fb), nameof(Brent));
            }

            if (Math.Abs(e) >= tolerance1 && Math.Abs(fa) > Math.Abs(fb))
            {
                // Attempt interpolation: linear when only two points are distinct,
                // inverse quadratic when three are.
                double numerator;
                double denominator;
                var s = fb / fa;

                if (a == c)
                {
                    numerator = 2.0 * bisectionStep * s;
                    denominator = 1.0 - s;
                }
                else
                {
                    var q = fa / fc;
                    var r = fb / fc;
                    numerator = s * (2.0 * bisectionStep * q * (q - r) - (b - a) * (r - 1.0));
                    denominator = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (numerator > 0.0)
                {
                    denominator = -denominator;
                }

                numerator = Math.Abs(numerator);

                // Accept the interpolated step only if it stays inside the bracket and
                // improves on the step before last; otherwise fall back to bisection.
                var bound = 3.0 * bisectionStep * denominator - Math.Abs(tolerance1 * denominator);
                var previousStep = Math.Abs(e * denominator);

                if (2.0 * numerator < Math.Min(bound, previousStep))
                {
                    e = d;
                    d = numerator / denominator;
                }
                else
                {
                    d = bisectionStep;
                    e = d;
                }
            }
            else
            {
                d = bisectionStep;
                e = d;
            }

            a = b;
            fa = fb;
            b += Math.Abs(d) > tolerance1
                ? d
                : bisectionStep >= 0.0 ? tolerance1 : -tolerance1;
            fb = function(b);
        }

        return new RootResult(b, false, maxIterations, Math.Abs(fb), nameof(Brent));
    }

    /// <summary>
    /// Newton's method using a finite-difference derivative. Converges quadratically near
    /// a simple root, but can diverge from a poor starting point and stalls where the
    /// derivative vanishes — both are reported rather than hidden.
    /// </summary>
    /// <param name="function">The function whose root is sought.</param>
    /// <param name="initialGuess">Starting point.</param>
    /// <param name="tolerance">Convergence tolerance on the step size and residual.</param>
    /// <param name="maxIterations">Maximum number of iterations.</param>
    /// <remarks>
    /// The derivative is approximated by the same backward difference that
    /// <see cref="DerivativeExtensions.Derivate(Func{double, double}, double, int)"/> uses at
    /// order 1, so results match. Prefer the overload taking an analytic derivative when one
    /// is available: it is both faster and more accurate.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="function"/> is null.</exception>
    public static RootResult Newton(
        Func<double, double> function,
        double initialGuess = 1.0,
        double tolerance = DefaultTolerance,
        int maxIterations = DefaultMaxIterations)
    {
        if (function == null)
        {
            throw new ArgumentNullException(nameof(function));
        }

        return Newton(
            function,
            x => (function(x) - function(x - DerivativeStep)) / DerivativeStep,
            initialGuess,
            tolerance,
            maxIterations);
    }

    /// <summary>
    /// Newton's method with an analytic derivative.
    /// </summary>
    /// <param name="function">The function whose root is sought.</param>
    /// <param name="derivative">The derivative of <paramref name="function"/>.</param>
    /// <param name="initialGuess">Starting point.</param>
    /// <param name="tolerance">Convergence tolerance on the step size and residual.</param>
    /// <param name="maxIterations">Maximum number of iterations.</param>
    /// <exception cref="ArgumentNullException"><paramref name="function"/> or <paramref name="derivative"/> is null.</exception>
    public static RootResult Newton(
        Func<double, double> function,
        Func<double, double> derivative,
        double initialGuess,
        double tolerance = DefaultTolerance,
        int maxIterations = DefaultMaxIterations)
    {
        ValidateCommon(function, tolerance, maxIterations);

        if (derivative == null)
        {
            throw new ArgumentNullException(nameof(derivative));
        }

        var x = initialGuess;
        var fx = function(x);

        for (var iteration = 1; iteration <= maxIterations; iteration++)
        {
            if (Math.Abs(fx) <= tolerance)
            {
                return new RootResult(x, true, iteration - 1, Math.Abs(fx), nameof(Newton));
            }

            var slope = derivative(x);

            // A vanishing or non-finite slope gives no usable step — stop instead of
            // dividing into an infinity and iterating on NaN.
            if (slope == 0.0 || !IsFinite(slope))
            {
                return new RootResult(x, false, iteration - 1, Math.Abs(fx), nameof(Newton));
            }

            var step = fx / slope;
            var next = x - step;

            if (!IsFinite(next))
            {
                return new RootResult(x, false, iteration - 1, Math.Abs(fx), nameof(Newton));
            }

            x = next;
            fx = function(x);

            if (Math.Abs(step) <= tolerance * (1.0 + Math.Abs(x)))
            {
                return new RootResult(x, true, iteration, Math.Abs(fx), nameof(Newton));
            }
        }

        return new RootResult(x, false, maxIterations, Math.Abs(fx), nameof(Newton));
    }

    private static void ValidateCommon(Func<double, double> function, double tolerance, int maxIterations)
    {
        if (function == null)
        {
            throw new ArgumentNullException(nameof(function));
        }

        if (tolerance <= 0.0 || !IsFinite(tolerance))
        {
            throw new ArgumentOutOfRangeException(nameof(tolerance), "Tolerance must be positive and finite.");
        }

        if (maxIterations < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(maxIterations), "At least one iteration is required.");
        }
    }

    private static void ValidateBracket(
        Func<double, double> function,
        ref double lower,
        ref double upper,
        double tolerance,
        int maxIterations,
        out double fLower,
        out double fUpper)
    {
        ValidateCommon(function, tolerance, maxIterations);

        if (!IsFinite(lower) || !IsFinite(upper))
        {
            throw new ArgumentException("The bracket endpoints must be finite.");
        }

        if (lower == upper)
        {
            throw new ArgumentException("The bracket endpoints must be distinct.");
        }

        if (lower > upper)
        {
            var swap = lower;
            lower = upper;
            upper = swap;
        }

        fLower = function(lower);
        fUpper = function(upper);

        if (!IsFinite(fLower) || !IsFinite(fUpper))
        {
            throw new ArgumentException("The function must be finite at both bracket endpoints.");
        }

        // A sign change guarantees a root for a continuous function; without one these
        // methods have nothing to bisect towards.
        if (fLower != 0.0 && fUpper != 0.0 && Math.Sign(fLower) == Math.Sign(fUpper))
        {
            throw new ArgumentException(
                $"The interval [{lower:G6}, {upper:G6}] does not bracket a sign change: " +
                $"f(lower) = {fLower:G6} and f(upper) = {fUpper:G6} have the same sign.");
        }
    }

    private static bool IsFinite(double value) => !double.IsNaN(value) && !double.IsInfinity(value);
}
