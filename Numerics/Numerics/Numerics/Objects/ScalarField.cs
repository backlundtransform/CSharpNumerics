using CSharpNumerics.Numerics.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.Objects
{
    /// <summary>
    /// Represents a scalar field f : ℝ³ → ℝ, the scalar counterpart of <see cref="VectorField"/>.
    /// <para>
    /// Wraps a <see cref="Func{Vector, Double}"/> and provides instance methods for the
    /// standard vector-calculus operations (∇, ∇², ∂/∂xᵢ) that are otherwise scattered
    /// across extension classes.
    /// </para>
    /// <para>
    /// An implicit conversion to <c>Func&lt;Vector, double&gt;</c> is provided so that a
    /// <see cref="ScalarField"/> can be passed to any API that accepts the raw delegate.
    /// </para>
    /// </summary>
    public struct ScalarField
    {
        public Func<Vector, double> f;


        public ScalarField(Func<Vector, double> func)
        {
            f = func;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Evaluation
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Evaluates the field at a point: f(r).
        /// </summary>
        public double Evaluate(Vector point) => f(point);

        /// <summary>
        /// Evaluates the field at a point given as a tuple.
        /// </summary>
        public double Evaluate((double, double, double) point) => f(new Vector(point));

        /// <summary>
        /// Evaluates the field over a diagonal line of points.
        /// </summary>
        public IDictionary<Vector, double> EvaluateRange(
            double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var result = new Dictionary<Vector, double>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                var point = new Vector(xmin + i, ymin + i, zmin + i);
                result.Add(point, f(point));
            }
            return result;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Differential operators
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes ∇f at a point: (∂f/∂x, ∂f/∂y, ∂f/∂z).
        /// </summary>
        public Vector Gradient((double, double, double) point) =>
            f.Gradient(point);

        /// <summary>
        /// Creates a gradient vector field over a diagonal line of points.
        /// </summary>
        public IDictionary<Vector, Vector> Gradient(
            double xmin, double ymin, double zmin, double stepSize, double maxSteps) =>
            f.Gradient(xmin, ymin, zmin, stepSize, maxSteps);

        /// <summary>
        /// Computes ∇²f at a point: ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z².
        /// </summary>
        public double Laplacian((double, double, double) point) =>
            f.Laplacian(point);

        /// <summary>
        /// Computes a partial derivative ∂ⁿf/∂xᵢⁿ at a point.
        /// </summary>
        public double Derivate(Vector variables, Cartesian cartesian, int order = 1) =>
            f.Derivate(variables, cartesian, order);

        /// <summary>
        /// Returns a <see cref="VectorField"/> representing ∇f — the gradient field.
        /// Each component of the returned field is ∂f/∂xᵢ evaluated lazily.
        /// </summary>
        public VectorField GradientField()
        {
            var func = f;
            return new VectorField(
                r => func.Derivate(r, Cartesian.x),
                r => func.Derivate(r, Cartesian.y),
                r => func.Derivate(r, Cartesian.z));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Integration
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Monte Carlo integration of the field over a 3D box.
        /// </summary>
        public double Integrate(Vector lowerLimit, Vector upperLimit) =>
            f.Integrate(lowerLimit, upperLimit);

        // ═══════════════════════════════════════════════════════════════
        //  Arithmetic
        // ═══════════════════════════════════════════════════════════════

        public static ScalarField operator +(ScalarField a, ScalarField b) =>
            new ScalarField(r => a.f(r) + b.f(r));

        public static ScalarField operator -(ScalarField a, ScalarField b) =>
            new ScalarField(r => a.f(r) - b.f(r));

        public static ScalarField operator -(ScalarField a) =>
            new ScalarField(r => -a.f(r));

        public static ScalarField operator *(ScalarField a, ScalarField b) =>
            new ScalarField(r => a.f(r) * b.f(r));

        public static ScalarField operator *(double a, ScalarField b) =>
            new ScalarField(r => a * b.f(r));

        public static ScalarField operator *(ScalarField a, double b) =>
            new ScalarField(r => a.f(r) * b);

        public static ScalarField operator /(ScalarField a, ScalarField b) =>
            new ScalarField(r => a.f(r) / b.f(r));

        public static ScalarField operator /(ScalarField a, double b) =>
            new ScalarField(r => a.f(r) / b);

        // ═══════════════════════════════════════════════════════════════
        //  Conversion
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Implicit conversion to the underlying delegate so that a <see cref="ScalarField"/>
        /// can be passed to any method that accepts <c>Func&lt;Vector, double&gt;</c>.
        /// </summary>
        public static implicit operator Func<Vector, double>(ScalarField field) => field.f;
    }
}
