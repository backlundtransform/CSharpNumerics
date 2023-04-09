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
            var features = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            var labels = new double[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
            int[] layers = new int[] { 10, 10, 10 };
            Func<double, double> activate = (x) => 1 / (1 + Math.Exp(-x));
          

            var neuralNetwork = new NeuralNetwork(layers, 0.01, activate);
            neuralNetwork.Train(features, labels, 1);
         

            var input =new double[] { 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
            var result = neuralNetwork.Predict(input);


            Assert.AreEqual(13, result[0],0.2);

        }

      
    }
}