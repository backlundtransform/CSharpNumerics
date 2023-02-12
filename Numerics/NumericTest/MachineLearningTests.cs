using Numerics.Objects;
using Xunit.Sdk;

namespace NumericTest
{
    [TestClass]
    public class MachineLearningTests
    {
        [TestMethod]
        public void MachineLearning()
        {
            Func<double, double> activateFunction = (x) => Math.Tanh(x);

            var layers = new int[] { 1, 0, 1 };


            var features = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

            var labels = new double[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

            var net = new NeuralNetwork();

            var model = net.Train(features, labels, 100);

            var result = net.Predict(model, 11);

            Assert.AreEqual(result, 13, 1);
        }
    }
}