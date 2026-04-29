using CSharpNumerics.Numerics.Enums;
using System;
using System.Collections.Generic;

namespace CSharpNumerics.Numerics.Objects
{
    public struct VectorField
    {
        public Func<Vector, double> fx;

        public Func<Vector, double> fy;

        public Func<Vector, double> fz;


        public VectorField(Func<Vector, double> u, Func<Vector, double> v, Func<Vector, double> w)
        {
            fx = u;
            fy = v;
            fz = w;

        }
        public IDictionary<Vector, Vector> Curl(double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var vectorField = new Dictionary<Vector, Vector>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                vectorField.Add(new Vector(xmin + i, ymin + i, zmin + i), Curl((xmin + i, ymin + i, zmin + i)));
            }
            return vectorField;

        }

        public IDictionary<Vector, Vector> EvaluateRange(double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var vectorField = new Dictionary<Vector, Vector>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                var parameters = new Vector(xmin + i, ymin + i, zmin + i);
                vectorField.Add(parameters, new Vector(fx(parameters), fy(parameters), fz(parameters)));
            }
            return vectorField;
        }

        /// <summary>
        /// Evaluates the vector field on a structured 2D grid (z = 0).
        /// Returns vectors at each (x, y) point for arrow-plot rendering.
        /// </summary>
        /// <param name="xmin">Minimum x coordinate.</param>
        /// <param name="xmax">Maximum x coordinate.</param>
        /// <param name="ymin">Minimum y coordinate.</param>
        /// <param name="ymax">Maximum y coordinate.</param>
        /// <param name="nx">Number of evaluation points in x.</param>
        /// <param name="ny">Number of evaluation points in y.</param>
        /// <returns>Dictionary mapping position vectors to field vectors on the 2D grid.</returns>
        public IDictionary<Vector, Vector> EvaluateGrid2D(double xmin, double xmax, double ymin, double ymax, int nx, int ny)
        {
            var vectorField = new Dictionary<Vector, Vector>();
            double dx = (nx > 1) ? (xmax - xmin) / (nx - 1) : 0.0;
            double dy = (ny > 1) ? (ymax - ymin) / (ny - 1) : 0.0;

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    var point = new Vector(xmin + i * dx, ymin + j * dy, 0.0);
                    vectorField.Add(point, new Vector(fx(point), fy(point), fz(point)));
                }
            }
            return vectorField;
        }


        public Vector Curl((double, double, double) points)
        {

            var dx = fz.Derivate(new Vector(points), Cartesian.y) - fy.Derivate(new Vector(points), Cartesian.z);
            var dy = fx.Derivate(new Vector(points), Cartesian.z) - fz.Derivate(new Vector(points), Cartesian.x);
            var dz = fy.Derivate(new Vector(points), Cartesian.x) - fx.Derivate(new Vector(points), Cartesian.y);
            return new Vector(dx, dy, dz);

        }

        public double Divergence((double, double, double) points) => fx.Derivate(new Vector(points), Cartesian.x) + fy.Derivate(new Vector(points), Cartesian.y) + fz.Derivate(new Vector(points), Cartesian.z);

    }
}