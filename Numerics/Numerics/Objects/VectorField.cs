using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Linq;

namespace Numerics.Objects
{
    public struct VectorField
    {
        public Func<Vector, double, double> fxt;

        public Func<Vector, double, double> fyt;

        public Func<Vector, double, double> fzt;

        public Func<double, double> xt;

        public Func<double, double> yt;

        public Func<double, double> zt;

        public Func<Vector, double> fx;

        public Func<Vector, double> fy;

        public Func<Vector, double> fz;




        public VectorField(Func<Vector, double> u, Func<Vector,  double> v, Func<Vector, double> w)
        {
            fx = u;
            fy =  v;
            fz =  w;
            xt = (double t) => 1;
            yt = (double t) => 1;
            zt = (double t) => 1;
            fxt = (Vector p, double t) => u(p);
            fyt = (Vector p, double t) => v(p);
            fzt = (Vector p, double t) => w(p);


        }

        public VectorField(Func<Vector, double, double> u, Func<Vector, double, double> v, Func<Vector, double, double> w, Func<double, double> x, Func<double, double> y, Func<double, double> z)
        {
            fxt = u;
            fyt = v;
            fzt = w;
            xt = x;
            yt = y;
            zt = z;
            fx = (Vector p)=>u(p,0);
            fy = (Vector p) => v(p, 0); 
            fz = (Vector p) => v(p, 0); 

        }

        public IDictionary<double,Vector> Velocity(double t)
        {

            var x0 = fxt(new Vector(xt(t),yt(t),zt(t)), t);
            var y0 = fyt(new Vector(xt(t), yt(t), zt(t)), t);
            var z0 = fzt(new Vector(xt(t), yt(t), zt(t)), t);

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
                    vectorField.Add(parameters,new Vector(fx(parameters), fy(parameters), fz(parameters)));            
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
                vectorField.Add(i, Curl((fx(parameters), fy(parameters), fz(parameters))));
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
