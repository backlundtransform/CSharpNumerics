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
            var features = new Tensor(new double[,] { { 1, 2 }, { 3, 4 } , { 5, 6 } });
            var labels = new Tensor(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } });
            int[] layers = new int[] { 4, 4, 4 };
            Func<double, double> activate = (x) => Math.Tanh(x);
          


            var neuralNetwork = new NeuralNetwork(layers, 0.1, activate);
            neuralNetwork.Train(features, labels, 1);
         

            var input = new Tensor(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } });
            var result = neuralNetwork.Predict(input);


            Assert.AreEqual(7, result[0],0.2);

        }

      
    }
}