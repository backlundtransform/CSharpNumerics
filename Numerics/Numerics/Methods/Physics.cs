using Numerics.Objects;
using System;
using System.Collections.Generic;
using System.Text;

namespace CSharpNumerics.Methods
{
    public class Physics
    {

        public static Vector ConvectiveDerivative(VectorField velocity, double time)
        {
            var u = velocity.xt.Derivate(time);
            var v = velocity.yt.Derivate(time);
            var w = velocity.zt.Derivate(time);
            var velocityVector = new Vector(u, v, w);
            var convectiveAccX = velocityVector.Dot(velocity.fx.Gradient((velocity.xt(time), velocity.yt(time), velocity.zt(time))));
            var convectiveAccY = velocityVector.Dot(velocity.fy.Gradient((velocity.xt(time), velocity.yt(time), velocity.zt(time))));
            var convectiveAccZ = velocityVector.Dot(velocity.fz.Gradient((velocity.xt(time), velocity.yt(time), velocity.zt(time))));

            Func<double, double> fx = (double t) => velocity.fxt(velocityVector, t);
            Func<double, double> fy = (double t) => velocity.fxt(velocityVector, t);
            Func<double, double> fz = (double t) => velocity.fxt(velocityVector, t);

            return new Vector(fx.Derivate(time) + convectiveAccX, fy.Derivate(time) + convectiveAccY, fz.Derivate(time) + convectiveAccZ);
        }


    }
}
