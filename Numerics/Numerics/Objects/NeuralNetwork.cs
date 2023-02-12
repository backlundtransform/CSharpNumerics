using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Objects
{

    public struct NeuralNetwork
    {
      
        public NeuralNetwork()
        {         

        }
      

        public double[] Train(double[] features, double[] labels, int epoch)
        {
            var weights = new double[features.Length];
            var bias = 0.0;
            var learningRate = 0.1;

            for (int i = 0; i < epoch; i++)
            {
                for (int j = 0; j < features.Length; j++)
                {
                   var prediction = Predict(weights, features[j]) + bias;
                   var error = labels[j] - prediction;

                    for (int k = 0; k < weights.Length; k++)
                    {
                        weights[k] = weights[k] + error * features[j] * learningRate;
                    }

                    bias += error * learningRate;
                }
            }

            return weights;
        }

        public double Predict(double[] model, double input)
        {
           var prediction = 0.0;

            for (var i = 0; i < model.Length; i++)
            {
                prediction += input * model[i];
            }

            return prediction;
        }


    }
}
