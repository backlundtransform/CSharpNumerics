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
            NeuralNetwork neuralNetwork = new NeuralNetwork(new int[]{ 10, 1, 10 }, 0.1, x => 1 / (1 + Math.Exp(-x)));

            for (int i = 0; i < 10000; i++)
            {
                neuralNetwork.Train(new double[] {  1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, new double[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },0);
              
            }

            var result = neuralNetwork.Predict(new double[] { 11, 12, 13,14,15,16,17,18,19,20,21});

            Assert.AreEqual(result[0], 13, 1);
        }
    }
}