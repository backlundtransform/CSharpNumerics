using CSharpNumerics.ML.Scalers;


namespace NumericTest
{

    [TestClass]
    public class StandardScalerTests
    {
        [TestMethod]
        public void Fit_ComputeMeansAndStdDevs_Correctly()
        {
            // Arrange
            double[,] data = {
            { 1.0, 2.0, 3.0 },
            { 4.0, 5.0, 6.0 },
            { 7.0, 8.0, 9.0 }
        };
            StandardScaler scaler = new StandardScaler();

            // Act
            scaler.Fit(data);

            // Assert
            Assert.IsTrue(4.0 == scaler.Means[0]);
            Assert.IsTrue(5.0 == scaler.Means[1]);
            Assert.IsTrue(6.0 == scaler.Means[2]);
            Assert.IsTrue(2.4495 == Math.Round(scaler.StdDevs[0], 4)); 
            Assert.IsTrue(2.4495 == Math.Round(scaler.StdDevs[1], 4)); 
            Assert.IsTrue(2.4495 == Math.Round(scaler.StdDevs[2], 4)); 
        }
        [TestMethod]
        public void Transform_StandardizeData_Correctly()
        {
            // Arrange
            double[,] data = {
            { 1.0, 2.0, 3.0 },
            { 4.0, 5.0, 6.0 },
            { 7.0, 8.0, 9.0 }
        };
            StandardScaler scaler = new StandardScaler();
            scaler.Fit(data);

            // Act
            double[,] transformedData = scaler.Transform(data);

            // Assert
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 0], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 1], 4)); 
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 2], 4)); 
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 0], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 1], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 2], 1));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 0], 4)); 
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 1], 4));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 2], 4)); 
        }

        [TestMethod]
        public void FitTransform_StandardizeData_Correctly()
        {
            // Arrange
            double[,] data = {
            { 1.0, 2.0, 3.0 },
            { 4.0, 5.0, 6.0 },
            { 7.0, 8.0, 9.0 }
        };
            StandardScaler scaler = new StandardScaler();

            // Act
            double[,] transformedData = scaler.FitTransform(data);

            // Assert
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 0], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 1], 4));
            Assert.IsTrue(-1.2247 == Math.Round(transformedData[0, 2], 4)); 
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 0], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 1], 1));
            Assert.IsTrue(0.0 == Math.Round(transformedData[1, 2], 1));
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 0], 4)); 
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 1], 4)); 
            Assert.IsTrue(1.2247 == Math.Round(transformedData[2, 2], 4)); 
        }

    }
}
