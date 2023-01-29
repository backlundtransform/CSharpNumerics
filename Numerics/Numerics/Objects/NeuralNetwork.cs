using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Objects
{

    public struct NeuralNetwork
    {
        public int[] layers;
        public Tensor neurons;
        public Tensor biases;
        public Tensor weights;
        public double fitness;
        public Func<double, double> activateFunction;
        public NeuralNetwork(int[] layers, double fitness, Func<double,double> activateFunction)
        {

            this.fitness = fitness;
            this.activateFunction =activateFunction;
            this.layers = new int[layers.Length];
       
            this.activateFunction = activateFunction;
            this.layers = new int[layers.Length];
            neurons = new Tensor(Array.Empty<double>());
            biases = new Tensor(Array.Empty<double>());

            weights = new Tensor(Array.Empty<double>());

        }
      

        public double[] Train(double[] features, double[] label, int epoch)
        {
            return new double[] { 1.0, 2.0 };
        }

        public double Predict(double[] model, double input)
        {
            return 13;
        }


    }
}
