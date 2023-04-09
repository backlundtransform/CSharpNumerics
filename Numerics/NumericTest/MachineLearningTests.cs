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
            var features = new double[] { 0, 0 };
            var labels = new double[] { 1 };
            int[] layers = new int[] { 2, 1, 1 };
            Func<double, double> activate = (x) => x;
          

            var neuralNetwork = new NeuralNetwork(layers, 0.1, activate);
            neuralNetwork.Train(features, labels, 1);
         

            var input =new double[] { 0, 1};
            var result = neuralNetwork.Predict(input);


            Assert.AreEqual(1, result[0],0.2);

        }

      
    }
}