using Numerics;
using Numerics.Enums;
using Numerics.Models;
using Numerics.Objects;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace System
{
    public static class CalculusExtensions
    {
        public const double h = 0.000001;

        public static double LeftLimit(this Func<double, double> function, double x)
        {
            return function(x - double.Epsilon);
        }

        public static double RightLimit(this Func<double, double> function, double x)
        {
            return function(x + double.Epsilon);
        }

        public static double Limit(this Func<double, double> function, double x)
        {
            var right =  function.LeftLimit(x);

            var left = function.RightLimit(x);

            return (right == left) ? right : double.NaN;
        }
        public static double Derivate(this Func<double, double> func, double x, int order = 1)
        {

            var h0 = Math.Pow(10, order) * h;

            var d = 0.0;

            var pascalMatrix = new Matrix(new double[order + 1, order + 1]).Pascal();

            for (var i = 1; i <= order + 1; i++)
            {

                d += Sign(i) * pascalMatrix.values[order, i - 1] * func(x + (order - i) * h0);

            }

            return d / Math.Pow(h0, order);
        }

        public static double Derivate(this Func<double[], double> func, double[] variables, int index, int order = 1)
        {
            Func<double, double> funcIndex = (double v) => func(variables.Select((c, i) => { if (i == index) { return v; } return c; }).ToArray());
            return funcIndex.Derivate(variables[index], order);

        }

        public static double Derivate(this Func<double, double> funcF,  Func<double, double> funcG,  double x, DerivateOperator derivateOperator = DerivateOperator.Chain, int order = 1)
        {

            switch (derivateOperator)
            {
               
                case DerivateOperator.Product:
                   return funcF.Derivate(x, order) * funcG(x) + funcF(x)*funcG.Derivate(x, order);
                case DerivateOperator.Quotient:
                    return (funcF.Derivate(x, order) * funcG(x) -funcF(x) * funcG.Derivate(x, order)) / Math.Pow(funcG(x),2);
                case DerivateOperator.Chain:
                default:
                    return funcF.Derivate(funcG(x), order) * funcG.Derivate(x, order);

            }

          
        }

        public static double Derivate(this Func<Vector, double> func, Vector variables, Cartesian cartesian, int order = 1)
        {

            switch (cartesian)
            {
                case Cartesian.x:
                    Func<double, double> funcX = (double v) => func(new Vector(v, variables.y, variables.z));
                    return funcX.Derivate(variables.x, order);

                case Cartesian.y:
                    Func<double, double> funcY = (double v) => func(new Vector(variables.x, v, variables.z));
                    return funcY.Derivate(variables.y, order);
                case Cartesian.z:
                    Func<double, double> funcZ = (double v) => func(new Vector(variables.x, variables.y, v));
                    return funcZ.Derivate(variables.z, order);

            }

            return func(variables);
        }


        public static double Derivate(this Func<(double x, double y), double> func, (double x, double y) variables, Cartesian cartesian, int order = 1)
        {

            switch (cartesian)
            {
                case Cartesian.x:
                    Func<double, double> funcX = (double v) => func((v, variables.y));
                    return funcX.Derivate(variables.x, order);

                case Cartesian.y:
                    Func<double, double> funcY = (double v) => func((variables.x, v));
                    return funcY.Derivate(variables.y, order);

            }

            return func(variables);
        }


        public static ComplexNumber Derivate(this ComplexFunction func, ComplexNumber variables, int order = 1)
        {

            var du = func.u.Derivate((variables.realPart, variables.imaginaryPart), Cartesian.x, order);

            var dv = func.v.Derivate((variables.realPart, variables.imaginaryPart), Cartesian.x, order);

            return (new ComplexNumber(du, dv));

        }

        public static IEnumerable<TimeSerie> Derivate<T>(this List<T> data, Func<T, (DateTime timeStamp, double value)> func)
        {
            return data.Select(func).Select(p => new TimeSerie() { TimeStamp = p.timeStamp, Value = p.value }).Derivate();
        }

        public static IEnumerable<Serie> Derivate<T>(this List<T> data, Func<T, (double index, double value)> func)
        {
            return data.Select(func).Select(p => new Serie() { Index = p.index, Value = p.value }).Derivate();
        }


        public static IEnumerable<Serie> Derivate(this IEnumerable<Serie> series)
        {
            var derivativesSeries = new List<Serie>();

            foreach (var serie in series) 
            {
                var deltaY = series.LinearInterpolation(p => (p.Index, p.Value), serie.Index + h) - serie.Value;
    
                derivativesSeries.Add(new Serie() { Value = deltaY / h,  Index =serie.Index});
            }
            return derivativesSeries;


        }


        public static IEnumerable<TimeSerie> Derivate(this IEnumerable<TimeSerie> series)
        {
            var derivativesSeries = new List<TimeSerie>();
            foreach (var serie in series)
            {
                var deltaY = series.LinearInterpolation(p => (p.TimeStamp.Ticks, p.Value), serie.TimeStamp.Ticks + h) - serie.Value;

                derivativesSeries.Add(new TimeSerie() { Value = deltaY / h, TimeStamp = serie.TimeStamp });
            }
            return derivativesSeries;


        }

        public static IDictionary<Vector, Vector> Gradient(this Func<Vector, double> func, double xmin, double ymin, double zmin, double stepSize, double maxSteps)
        {
            var vectorField = new Dictionary<Vector, Vector>();
            for (var i = 0.0; i < maxSteps; i += stepSize)
            {
                vectorField.Add(new Vector((xmin + i, ymin + i, zmin + i)), func.Gradient((xmin + i, ymin + i, zmin + i)));
            }
            return vectorField;

        }

        public static Vector Gradient(this Func<Vector, double> func, (double, double, double) points)
        {
            var dx = func.Derivate(new Vector(points), Cartesian.x);
            var dy = func.Derivate(new Vector(points), Cartesian.y);
            var dz = func.Derivate(new Vector(points), Cartesian.z);
            return new Vector(dx, dy, dz);

        }

    

        public static double Laplacian(this Func<Vector, double> func, (double, double, double) points)
        {
            var dx2 = func.Derivate(new Vector(points), Cartesian.x, 2);
            var dy2 = func.Derivate(new Vector(points), Cartesian.y, 2);
            var dz2 = func.Derivate(new Vector(points), Cartesian.z, 2);
            return dx2 + dy2 + dz2;

        }

        public static double Integrate(this Func<double, double> func, double lowerLimit, double upperLimit)
        {
            var result = 0.0;
            for (var i = lowerLimit; i <= upperLimit; i += h)
            {
                result += h * (func(i) + func(i + h)) / 2;
            }
            return result;
        }


        public static double Integrate(this Func<(double x, double y), double> func, (double lowerLimit, double upperLimit) xlimit, (double lowerLimit, double upperLimit) ylimit)
        {
            Func<Vector, double> funcv = (Vector v) => func((v.x, v.y));

            return funcv.Integrate(new Vector(xlimit.lowerLimit, ylimit.lowerLimit, 1), new Vector(xlimit.upperLimit, ylimit.upperLimit, 2));
        }


        public static double Integrate(this Func<Vector, double> func, Vector lowerLimit, Vector upperLimit)
        {

            var rnd = new Random();

            var result = 0.0;
            var throws = 999999;


            for (var i = 0; i < throws; i++)
            {

                var x = rnd.NextDouble() * (upperLimit.x - lowerLimit.x) + lowerLimit.x;
                var y = rnd.NextDouble() * (upperLimit.y - lowerLimit.y) + lowerLimit.y;
                var z = rnd.NextDouble() * (upperLimit.z - lowerLimit.z) + lowerLimit.z;

                result += func(new Vector(x, y, z));

            }
            return (upperLimit.x - lowerLimit.x) * (upperLimit.y - lowerLimit.y) * (upperLimit.z - lowerLimit.z) * result / 999999;

        }


        public static ComplexNumber Integrate(this ComplexFunction func, (double lowerLimit, double upperLimit) xlimit, (double lowerLimit, double upperLimit) ylimit)
        {
            Func<Vector, double> funcu = (Vector u) => func.u((u.x, u.y));
            Func<Vector, double> funcv = (Vector v) => func.v((v.x, v.y));

            var re = funcu.Integrate(new Vector(xlimit.lowerLimit, ylimit.lowerLimit, 1), new Vector(xlimit.upperLimit, ylimit.upperLimit, 2));

            var img = funcv.Integrate(new Vector(xlimit.lowerLimit, ylimit.lowerLimit, 1), new Vector(xlimit.upperLimit, ylimit.upperLimit, 2));

            return new ComplexNumber(re, img);
        }

        public static double Integrate<T>(this List<T> data, Func<T, (DateTime timeStamp, double value)> func)
        {
            return data.Select(func).Select(p => new TimeSerie() { TimeStamp = p.timeStamp, Value = p.value }).Integrate();
        }

        public static double Integrate<T>(this List<T> data, Func<T, (double index, double value)> func)
        {
            return data.Select(func).Select(p => new Serie() { Index = p.index, Value = p.value }).Integrate();
        }


        public static double Integrate(this List<TimeSerie> data, DateTime start, DateTime end)
        {
            return data.Select(p => new Serie() { Index = (p.TimeStamp-start).TotalSeconds, Value = p.Value }).ToList().Integrate(0.0, (end - start).TotalSeconds);
        }

        public static double Integrate(this List<Serie> data, double start, double end)
        {
            var sum = 0.0;

            if (data.Count == 1)
            {
                return (end - start)*data[0].Value;
            }

            sum += ((data[0].Index - start) + (data[1].Index - data[0].Index)/ 2) * data[0].Value;

            for (var i = 1; i < data.Count - 1; i++)
            {
                sum += ((data[i].Index - data[i - 1].Index + (data[i + 1].Index - data[i].Index)) / 2) * data[i].Value;
            }

            sum += (data[data.Count - 1].Index - data[data.Count - 2].Index) / 2 + (end - data[data.Count - 1].Index) * data[data.Count - 1].Value;

            return sum;
        }



        public static double Integrate(this IEnumerable<TimeSerie> data)
        {
            return data.Select(p => new Serie() { Index = (p.TimeStamp - data.First().TimeStamp).TotalSeconds, Value = p.Value }).ToList().Integrate();

        }

        public static double Integrate(this IEnumerable<Serie> data)
        {
            var sum = 0.0;
            var dataList = data.ToList();

            for (var i = 0; i < dataList.Count - 1; i++)
            {
                sum += (dataList[i + 1].Index - dataList[i].Index) * (dataList[i].Value + dataList[i + 1].Value) / 2;
            }

            return sum;
        }

        private static int Sign(int index) => index % 2 == 0 ? -1 : 1;
    }
}

