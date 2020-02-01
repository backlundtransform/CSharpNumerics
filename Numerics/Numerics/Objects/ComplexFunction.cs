using System;
using System.Collections.Generic;
using System.Text;

namespace Numerics.Objects
{


    public class ComplexFunction
    {
        public Func<(double x, double y), double> u;

        public Func<(double x, double y), double> v;

      

        public ComplexFunction(Func<(double x, double y), double> re, Func<(double x, double y), double> im)
        {
            u = re;
            v = im;
        }

        public Matrix Jacobian((double x, double y) point) {


            var udx = u.Derivate((point.x, point.y), Cartesian.x);
            var vdy = v.Derivate((point.x, point.y), Cartesian.y);
            var udy = u.Derivate((point.x, point.y), Cartesian.y);
            var vdx = v.Derivate((point.x, point.y), Cartesian.x);

            return new Matrix(new double[,] { { udx, udy }, { vdx, vdy } });
        }

        public bool IsAnalytical((double x, double y) point)
        {

            var jacobian = Jacobian(point);

            var udx = u.Derivate((point.x, point.y), Cartesian.x);
            var vdy = v.Derivate((point.x, point.y), Cartesian.y);
            var udy = u.Derivate((point.x, point.y), Cartesian.y);
            var vdx = v.Derivate((point.x, point.y), Cartesian.x);       

            return Math.Round(jacobian.values[0,0], 2)== Math.Round(jacobian.values[1,1], 2)
                 && Math.Round(jacobian.values[0, 1], 2) == -Math.Round(jacobian.values[1, 0], 2);

        }
      

    }

}
