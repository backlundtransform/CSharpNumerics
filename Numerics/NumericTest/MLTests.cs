using CSharpNumerics.ML;
using CSharpNumerics.ML.Models.Regression;
using CSharpNumerics.ML.Scalers;
using CSharpNumerics.ML.Selector;
using CSharpNumerics.Objects;
using Numerics.Objects;


namespace NumericTest
{

    [TestClass]
    public class MLTests
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
            Matrix transformedDataMatrix = scaler.FitTransform(new Matrix(data));

            var transformedData = transformedDataMatrix.values;

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
            Matrix transformedDataMatrix = scaler.FitTransform(new Matrix(data));

            var transformedData = transformedDataMatrix.values;

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
        public void TestCrossValidatorPipeline()
        {
           
            int nSamples = 100;
            int nFeatures = 2;

            Matrix X = new Matrix(nSamples, nFeatures);
            VectorN y = new VectorN(nSamples);

            Random rnd = new Random(123);
            for (int i = 0; i < nSamples; i++)
            {
                double x1 = rnd.NextDouble() * 10;
                double x2 = rnd.NextDouble() * 10;
                X.values[i, 0] = x1;
                X.values[i, 1] = x2;
                y[i] = 3 * x1 + 2 * x2 + rnd.NextDouble() * 0.1;
            }

            var pipeline = new Pipeline(new Linear(), new() { ["LearningRate"] = 0.01 }, selector: new SelectKBest(), selectorParams: new() { ["K"] = 1 });

            var cv = new RollingCrossValidator([pipeline], 5);
            var results = cv.Run(X, y);


            // 5. Skriv ut resultat
            Console.WriteLine("Cross validation results:");
           

        }

    }
}
