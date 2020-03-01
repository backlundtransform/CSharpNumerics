namespace System
{
    public static class DifferentialEquationExtensions
    {

        public static double RungeKutta(this Func<(double t, double y), double> func, double min, double max, double stepSize, double yInitial)
        {

            var y = yInitial;

            for (var t = min+stepSize; t <= max; t += stepSize)
            {
                var k1 = func((t, y));

                var k2 = func((t + 2.0 / 3.0 * stepSize, y + 2.0 / 3.0 * stepSize * k1));
                y += stepSize * (1.0 / 4.0 * k1 + 3.0 / 4.0 * k2);

            }
            return y;

        }
    }
}
