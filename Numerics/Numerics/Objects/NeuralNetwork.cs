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
            for (var i = 0; i < layers.Length; i++)
            {
                this.layers[i] = layers[i];
            }

            var neuronList = new List<double[]>();

            for (var i = 0; i < layers.Length; i++)
            {
                neuronList.Add(new double[layers[i]]);
            }
            neurons = new Tensor(neuronList.ToArray());
            var biasesList = new List<double[]>();
            for (var i = 0; i < layers.Length; i++)
            {
                var bias = new double[layers[i]];
                for (var j = 0; j < layers[i]; j++)
                {
                   
                    bias[j] = new Random().RandomDouble(-0.5f, 0.5f); 
                }
                biasesList.Add(bias);
            }

           biases = new Tensor(biasesList.ToArray());
           var weightsList = new List<double[][]>();

            for (var i = 1; i < layers.Length; i++)
            {
                var layerWeightsList = new List<double[]>();
                var  neuronsInPreviousLayer = layers[i - 1];

                for (var j = 0; j < neurons.values.GetLength(i); j++)
                {
                    var neuronWeights = new double[neuronsInPreviousLayer];
                    for (var k = 0; k < neuronsInPreviousLayer; k++)
                    {
                        
                        neuronWeights[k] = new Random().RandomDouble(-0.5f, 0.5f);
                    }
                    layerWeightsList.Add(neuronWeights);
                }
                weightsList.Add(layerWeightsList.ToArray());
            }
            weights = new Tensor(weightsList.ToArray());

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
                neurons.values.SetValue(inputs[i],0, i);
            }
            for (var i = 1; i < layers.Length; i++)
            {

                for (var j = 0; j < neurons.values.GetLength(i); j++)
                {
                    var value = 0.0;
                    for (var k = 0; k < neurons.values.GetLength(i - 1); k++)
                    {
                        value += Convert.ToDouble(weights.values.GetValue(i - 1, j, k)) * Convert.ToDouble(neurons.values.GetValue(i - 1, k));
                    }
                    neurons.values.SetValue(activateFunction(value + Convert.ToDouble(biases.values.GetValue(i,j))),i,j) ;
                }
            }
            return (double[])neurons.values.GetValue(neurons.values.Length - 1);
        }


        public double[] Train(double[] features, double label, int epoch)
        {
            throw new Exception();
        }


    }
}
