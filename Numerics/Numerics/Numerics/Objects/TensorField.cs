using CSharpNumerics.Numerics.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.Objects
{
    /// <summary>
    /// Represents a rank-2 tensor field T : ℝ³ → ℝ³ˣ³, mapping each spatial point
    /// to a 3×3 <see cref="Matrix"/>.
    /// <para>
    /// Completes the field hierarchy:
    /// <see cref="ScalarField"/> (rank 0) → <see cref="VectorField"/> (rank 1) → <see cref="TensorField"/> (rank 2).
    /// </para>
    /// <para>
    /// Typical uses: stress/strain tensors, inertia tensor fields, the Maxwell stress tensor,
    /// and the Jacobian of a vector field.
    /// </para>
    /// </summary>
    public struct TensorField
    {
        public Func<Vector, Matrix> f;


        public TensorField(Func<Vector, Matrix> func)
        {
            f = func;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Evaluation
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Evaluates the tensor field at a point: T(r).
        /// </summary>
        public Matrix Evaluate(Vector point) => f(point);

        /// <summary>
        /// Evaluates the tensor field at a point given as a tuple.
        /// </summary>
        public Matrix Evaluate((double, double, double) point) => f(new Vector(point));

        /// <summary>
        /// Evaluates the field over a diagonal line of points.
        /// </summary>
        public IDictionary<Vector, Matrix> EvaluateRange(
            double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var result = new Dictionary<Vector, Matrix>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                var point = new Vector(xmin + i, ymin + i, zmin + i);
                result.Add(point, f(point));
            }
            return result;
        }

        // ═══════════════════════════════════════════════════════════════
        //  Component access
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Extracts component Tᵢⱼ as a <see cref="ScalarField"/>.
        /// </summary>
        public ScalarField Component(int i, int j)
        {
            var func = f;
            return new ScalarField(r => func(r).values[i, j]);
        }

        /// <summary>
        /// Extracts row i as a <see cref="VectorField"/>: (Tᵢ₁, Tᵢ₂, Tᵢ₃).
        /// </summary>
        public VectorField Row(int i)
        {
            var func = f;
            return new VectorField(
                r => func(r).values[i, 0],
                r => func(r).values[i, 1],
                r => func(r).values[i, 2]);
        }

        // ═══════════════════════════════════════════════════════════════
        //  Differential operators
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes the divergence of the tensor field: (∇·T)ᵢ = Σⱼ ∂Tᵢⱼ/∂xⱼ.
        /// Returns a <see cref="VectorField"/>.
        /// </summary>
        public VectorField Divergence()
        {
            var func = f;
            var t00 = Cmp(func, 0, 0); var t01 = Cmp(func, 0, 1); var t02 = Cmp(func, 0, 2);
            var t10 = Cmp(func, 1, 0); var t11 = Cmp(func, 1, 1); var t12 = Cmp(func, 1, 2);
            var t20 = Cmp(func, 2, 0); var t21 = Cmp(func, 2, 1); var t22 = Cmp(func, 2, 2);

            return new VectorField(
                r => t00.Derivate(r, Cartesian.x) + t01.Derivate(r, Cartesian.y) + t02.Derivate(r, Cartesian.z),
                r => t10.Derivate(r, Cartesian.x) + t11.Derivate(r, Cartesian.y) + t12.Derivate(r, Cartesian.z),
                r => t20.Derivate(r, Cartesian.x) + t21.Derivate(r, Cartesian.y) + t22.Derivate(r, Cartesian.z));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Contractions
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Computes the trace: Tr(T) = Σᵢ Tᵢᵢ. Returns a <see cref="ScalarField"/>.
        /// </summary>
        public ScalarField Trace()
        {
            var func = f;
            return new ScalarField(r =>
            {
                var m = func(r);
                double sum = 0;
                int n = Math.Min(m.rowLength, m.columnLength);
                for (int i = 0; i < n; i++)
                    sum += m.values[i, i];
                return sum;
            });
        }

        /// <summary>
        /// Contracts the tensor field with a <see cref="VectorField"/>: (T·v)ᵢ = Σⱼ Tᵢⱼ vⱼ.
        /// Returns a <see cref="VectorField"/>.
        /// </summary>
        public VectorField Contract(VectorField v)
        {
            var func = f;
            return new VectorField(
                r => { var m = func(r); return m.values[0, 0] * v.fx(r) + m.values[0, 1] * v.fy(r) + m.values[0, 2] * v.fz(r); },
                r => { var m = func(r); return m.values[1, 0] * v.fx(r) + m.values[1, 1] * v.fy(r) + m.values[1, 2] * v.fz(r); },
                r => { var m = func(r); return m.values[2, 0] * v.fx(r) + m.values[2, 1] * v.fy(r) + m.values[2, 2] * v.fz(r); });
        }

        /// <summary>
        /// Double contraction T : S = Σᵢⱼ Tᵢⱼ Sᵢⱼ. Returns a <see cref="ScalarField"/>.
        /// </summary>
        public ScalarField DoubleContract(TensorField other)
        {
            var funcA = f;
            var funcB = other.f;
            return new ScalarField(r =>
            {
                var a = funcA(r);
                var b = funcB(r);
                double sum = 0;
                for (int i = 0; i < a.rowLength; i++)
                    for (int j = 0; j < a.columnLength; j++)
                        sum += a.values[i, j] * b.values[i, j];
                return sum;
            });
        }

        // ═══════════════════════════════════════════════════════════════
        //  Transpose
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Returns the pointwise transpose: Tᵀᵢⱼ(r) = Tⱼᵢ(r).
        /// </summary>
        public TensorField Transpose()
        {
            var func = f;
            return new TensorField(r => func(r).Transpose());
        }

        // ═══════════════════════════════════════════════════════════════
        //  Factories
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Creates a <see cref="TensorField"/> representing the Jacobian of a vector field:
        /// (∇F)ᵢⱼ = ∂Fᵢ/∂xⱼ.
        /// </summary>
        public static TensorField FromJacobian(VectorField F)
        {
            return new TensorField(r => new Matrix(new double[,]
            {
                { F.fx.Derivate(r, Cartesian.x), F.fx.Derivate(r, Cartesian.y), F.fx.Derivate(r, Cartesian.z) },
                { F.fy.Derivate(r, Cartesian.x), F.fy.Derivate(r, Cartesian.y), F.fy.Derivate(r, Cartesian.z) },
                { F.fz.Derivate(r, Cartesian.x), F.fz.Derivate(r, Cartesian.y), F.fz.Derivate(r, Cartesian.z) }
            }));
        }

        /// <summary>
        /// Creates a constant identity tensor field: Tᵢⱼ(r) = δᵢⱼ.
        /// </summary>
        public static TensorField Identity()
        {
            return new TensorField(_ => new Matrix(new double[,]
            {
                { 1, 0, 0 },
                { 0, 1, 0 },
                { 0, 0, 1 }
            }));
        }

        // ═══════════════════════════════════════════════════════════════
        //  Arithmetic
        // ═══════════════════════════════════════════════════════════════

        public static TensorField operator +(TensorField a, TensorField b)
        {
            var fa = a.f; var fb = b.f;
            return new TensorField(r =>
            {
                var ma = fa(r); var mb = fb(r);
                var result = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        result[i, j] = ma.values[i, j] + mb.values[i, j];
                return new Matrix(result);
            });
        }

        public static TensorField operator -(TensorField a, TensorField b)
        {
            var fa = a.f; var fb = b.f;
            return new TensorField(r =>
            {
                var ma = fa(r); var mb = fb(r);
                var result = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        result[i, j] = ma.values[i, j] - mb.values[i, j];
                return new Matrix(result);
            });
        }

        public static TensorField operator -(TensorField a)
        {
            var fa = a.f;
            return new TensorField(r =>
            {
                var m = fa(r);
                var result = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        result[i, j] = -m.values[i, j];
                return new Matrix(result);
            });
        }

        public static TensorField operator *(double a, TensorField b)
        {
            var fb = b.f;
            return new TensorField(r =>
            {
                var m = fb(r);
                var result = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        result[i, j] = a * m.values[i, j];
                return new Matrix(result);
            });
        }

        public static TensorField operator *(TensorField a, double b) => b * a;

        public static TensorField operator *(ScalarField s, TensorField t)
        {
            var fs = s.f; var ft = t.f;
            return new TensorField(r =>
            {
                var sv = fs(r); var m = ft(r);
                var result = new double[3, 3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        result[i, j] = sv * m.values[i, j];
                return new Matrix(result);
            });
        }

        // ═══════════════════════════════════════════════════════════════
        //  Conversion
        // ═══════════════════════════════════════════════════════════════

        /// <summary>
        /// Implicit conversion to the underlying delegate.
        /// </summary>
        public static implicit operator Func<Vector, Matrix>(TensorField field) => field.f;


        private static Func<Vector, double> Cmp(Func<Vector, Matrix> func, int i, int j) =>
            r => func(r).values[i, j];
    }
}
