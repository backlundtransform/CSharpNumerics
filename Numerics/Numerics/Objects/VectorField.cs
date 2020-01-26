﻿using System;
using System.Collections.Generic;

namespace Numerics.Objects
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
        public IEnumerable<Vector> Curl(double xMin, double yMin, double zMin, double stepSize, double maxSteps)
        {
            var vectorField = new List<Vector>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                vectorField.Add(Curl((xMin + i, yMin + i, zMin + i)));
            }
            return vectorField;

        }

        public IEnumerable<Vector> EvaluateRange(double xMin, double yMin, double zMin, double stepSize, double maxSteps)
        {
            var vectorField = new List<Vector>();
            for (var i = 0.0; i < maxSteps; i+=stepSize)
            {
                    var parameters = new Vector(xMin + i, yMin + i, zMin + i);
                    vectorField.Add(new Vector(fx(parameters), fy(parameters), fz(parameters)));            
            }
            return vectorField;
        }


        public Vector Curl((double, double, double) points) {

            var dx = fz.Derivate(new Vector(points), Cartesian.y)- fy.Derivate(new Vector(points), Cartesian.z);
            var dy = fx.Derivate(new Vector(points), Cartesian.z) - fz.Derivate(new Vector(points), Cartesian.x);
            var dz = fy.Derivate(new Vector(points), Cartesian.x)- fx.Derivate(new Vector(points), Cartesian.y);
            return new Vector(dx, dy, dz);

        }

        public double Divergence((double, double, double) points)=> fx.Derivate(new Vector(points), Cartesian.x)+ fy.Derivate(new Vector(points), Cartesian.y)+fz.Derivate(new Vector(points), Cartesian.z);

    }
}
