using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace Numerics.Objects
{

    public struct NeuralNetwork
    {
        public int[] layers;
        public List<double[]> neurons;
        public List<double[]> biases;
        public List<double[][]> weights;
        public double fitness;
        public Func<double, double> activateFunction;
        public  NeuralNetwork(int[] layers, double fitness, Func<double,double> activateFunction)
        {

            this.fitness = fitness;
            this.activateFunction =activateFunction;
            this.layers = new int[layers.Length];
            for (var i = 0; i < layers.Length; i++)
            {
                this.layers[i] = layers[i];
            }
            this.neurons = new List<double[]>();
            for (var i = 0; i < layers.Length; i++)
            {
                neurons.Add(new double[layers[i]]);
            }
            this.biases = new List<double[]>();
            for (var i = 0; i < layers.Length; i++)
            {
                var bias = new double[layers[i]];
                for (var j = 0; j < layers[i]; j++)
                {
                    var rnd = new Random();
                    bias[j] = rnd.Next(-1, 0) + 0.5f;
                }
                biases.Add(bias);
            }


            this.weights = new List<double[][]>();
            for (var i = 1; i < layers.Length; i++)
            {
                var layerWeightsList = new List<double[]>();
                var  neuronsInPreviousLayer = layers[i - 1];
                for (var j = 0; j < neurons[i].Length; j++)
                {
                    var neuronWeights = new double[neuronsInPreviousLayer];
                    for (var k = 0; k < neuronsInPreviousLayer; k++)
                    {
                        
                        neuronWeights[k] = RandomDouble(-0.5f, 0.5f);
                    }
                    layerWeightsList.Add(neuronWeights);
                }
                weights.Add(layerWeightsList.ToArray());
            }

        }
        public int CompareTo(NeuralNetwork other)
        {

            if (fitness > other.fitness)
            {
                return 1;
            }
                  
            if (fitness < other.fitness)
            {
                return -1;
            }
                 
                return 0;
        }
        public double[] FeedForward(double[] inputs)
        {
            for (var i = 0; i < inputs.Length; i++)
            {
                neurons[0][i] = inputs[i];
            }
            for (var i = 1; i < layers.Length; i++)
            {

                for (var j = 0; j < neurons[i].Length; j++)
                {
                    var value = 0.0;
                    for (var k = 0; k < neurons[i - 1].Length; k++)
                    {
                        value += weights[i - 1][j][k] * neurons[i - 1][k];
                    }
                    neurons[i][j] = activateFunction(value + biases[i][j]);
                }
            }
            return neurons[neurons.Count - 1];
        }

 
        public void Mutate(int chance, double val)
        {
            for (int i = 0; i < biases.Count; i++)
            {
                for (int j = 0; j < biases[i].Length; j++)
                {
  
                    biases[i][j] = (RandomDouble(0, chance) <= 5) ? biases[i][j] += RandomDouble(-val, val) : biases[i][j];
                }
            }

            for (int i = 0; i < weights.Count; i++)
            {
                for (int j = 0; j < weights[i].Length; j++)
                {
                    for (int k = 0; k < weights[i][j].Length; k++)
                    {
                        weights[i][j][k] = (RandomDouble(0f, chance) <= 5) ? weights[i][j][k] += RandomDouble(-val, val) : weights[i][j][k];

                    }
                }
            }
        }

        private double RandomDouble(double min, double max) 
        {
           var rnd = new Random();
           return rnd.NextDouble() * (max - min) + min;
        }

    }
}
