using CSharpNumerics.Statistics.Data;
using Xunit.Sdk;


namespace NumericsTests
{
    
    [TestClass]
    public class CalculusTests
    {
        
        private const double timeZero = 5;
       
        
        private static readonly Func<double, double> acceleration =(double time) => 9.8;
        private static readonly Func<double, double> velocity= (double time) => acceleration(time) * time;
        private static readonly Func<double, double> displacement = (double time) => acceleration(time) * Math.Pow(time, 2) / 2;

        //Derivative
        [TestMethod]
        public void TestVelocityToAcceleration()
        {
            var result = velocity.Derivate(timeZero);
            Assert.IsTrue(Math.Round(result, 1) == acceleration(timeZero));
        }

        [TestMethod]
        public void TestVelocityToAccelerationSeries()
        {
          
            var series = velocity.GetSeries(0, 10, 1000).Derivate();

            var (slope, intercept, correlation) = series.LinearRegression(p => (p.Index, p.Value));

            Assert.IsTrue(Math.Round(slope, 1) == 0 && Math.Round(intercept, 1) == acceleration(timeZero));
        }


        [TestMethod]
        public void TestDisplacementToVelocity()
        {
         
            var result = displacement.Derivate(timeZero);

            Assert.IsTrue(Math.Round(result, 1) == Math.Round(acceleration(timeZero) * timeZero, 1));
        }

        [TestMethod]
        public void TestDisplacementVelocitySeries()
        {
          
            var series = displacement.GetSeries(0, 10, 1000).Derivate();

            var (slope, intercept, correlation) = series.LinearRegression(p => (p.Index, p.Value));

            Assert.IsTrue(Math.Round(slope, 1) == acceleration(timeZero) && Math.Round(intercept) == 0);
        }

        //Integral

        [TestMethod]
        public void TestAccelerationToVelocity()
        {
        
            var result = acceleration.Integrate(0,timeZero);
            Assert.IsTrue(Math.Round(result, 1) == Math.Round(velocity(timeZero),1));
        }

        [TestMethod]
        public void TestAccelerationToVelocitySeries()
        {
            var result = acceleration.GetSeries(0, timeZero, 1000).Integrate();
            Assert.IsTrue(Math.Round(result, 1) == Math.Round(velocity(timeZero), 1));
        }



        [TestMethod]
        public void TestVelocityToDisplacement()
        {
          
            var result = velocity.Integrate(0, timeZero);
            Assert.IsTrue(Math.Round(result, 1) == Math.Round(displacement(timeZero), 1));
        }

        [TestMethod]
        public void TestVelocityToDisplacementSeries()
        {
            var result = velocity.GetSeries(0, timeZero, 1000).Integrate();
            Assert.IsTrue(Math.Round(result, 1) == Math.Round(displacement(timeZero), 1)); 
        }


        //Differential equation

        [TestMethod]
        public void TestVelocityWithOutDrag()
        {

            Func<(double time, double velocity), double> func = ((double time, double velocity) v) => acceleration(v.time);
            var result = func.EulerMetod(0, timeZero, 0.00025, 0);
  
            Assert.IsTrue(Math.Round(result, 1) == Math.Round(velocity(timeZero), 1));

        }

        [TestMethod]
        public void TestDisplacementWithOutDrag()
        {

            Func<(double time, double displacement), double> func = ((double time, double displacement) v) => velocity(v.time);
            var result = func.RungeKutta(0, timeZero, 0.00025, 0);

            Assert.IsTrue(Math.Round(result, 1) == Math.Round(displacement(timeZero), 1));

        }
        [TestMethod]
        public void TestVelocityWithDrag()
        {

            var dragCoefficient = 0.5;
            var mass= 5;
            Func<(double time, double velocity), double> func = ((double time, double velocity) v) =>acceleration(v.time) - dragCoefficient * v.velocity / mass;
            var result = func.RungeKutta(0, timeZero, 0.00025, 0);

            var velocityDragFunction = mass * acceleration(timeZero) / dragCoefficient * (1 - Math.Exp(- dragCoefficient * timeZero/mass));

            Assert.IsTrue(Math.Round(result, 1) == Math.Round(velocityDragFunction, 1));    

        }

 
    }
}
