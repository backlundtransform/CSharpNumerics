using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Linq;

namespace Numerics.Objects
{
    public struct VectorField
    {
        public Func<Vector, double, double> fx;

        public Func<Vector, double, double> fy;

        public Func<Vector, double, double> fz;

        public Func<double, double> xt;

        public Func<double, double> yt;

        public Func<double, double> zt;




        public VectorField(Func<Vector, double> u, Func<Vector,  double> v, Func<Vector, double> w)
        {
            fx = (Vector r, double t)=>u(r);
            fy = (Vector r, double t) => v(r);
            fz = (Vector r, double t) => w(r);
            xt = null;
            yt =null;
            zt= null;


        }

        public VectorField(Func<Vector, double, double> u, Func<Vector, double, double> v, Func<Vector, double, double> w, Func<double, double> x, Func<double, double> y, Func<double, double> z)
        {
            fx = u;
            fy = v;
            fz = w;
            xt = x;
            yt = y;
            zt = z;

        }

        public IDictionary<double,Vector> Velocity(double t)
        {

            var x0 = fx(new Vector(xt(t),yt(t),zt(t)), t);
            var y0 = fy(new Vector(xt(t), yt(t), zt(t)), t);
            var z0 = fz(new Vector(xt(t), yt(t), zt(t)), t);

            return new Dictionary<double, Vector>() { { t, new Vector(x0, y0, z0) } };

        }
        public IDictionary<double, Vector> Velocity(double tmin, double stepSize, double maxSteps)
        {

            var vectorField = new Dictionary<double, Vector>();
            for (var i = tmin; i < tmin + maxSteps; i += stepSize)
            {
                vectorField.Add(i,Velocity(i)[i]);
            }
            return vectorField;

        }

        public IDictionary<Vector, Vector> EvaluateRange(double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var vectorField = new Dictionary<Vector, Vector>();
            for (var i = 0.0; i < maxSteps; i+=stepSize)
            {
                    var parameters = new Vector(xmin + i, ymin + i, zmin + i);
                    vectorField.Add(parameters,new Vector(fx(parameters,0), fy(parameters,0), fz(parameters,0)));            
            }
            return vectorField;
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

        public IDictionary<double, Vector> Curl(double tmin, double stepSize, double maxSteps)
        {
            var vectorField = new Dictionary<double, Vector>();
            for (var i = tmin; i < tmin + maxSteps; i += stepSize)
            {
                var parameters = Velocity(i)[i];
                vectorField.Add(i, Curl((fx(parameters,0), fy(parameters,0), fz(parameters,0))));
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
