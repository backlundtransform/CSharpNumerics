using System;
using System.Linq;

public  class NeuralNetwork
{
    int[] layers;
    double learningRate;
    Func<double, double> activateFunction;

    double weight;
    double bias;

    public NeuralNetwork(int[] layers, double learningRate, Func<double, double> activateFunction)
    {
        this.layers = layers;
        this.learningRate = learningRate;
        this.activateFunction = activateFunction;
    }

    public void Train(double[] features, double[] labels, int epochs)
    {
        int n = features.Length;

        // Initialize model parameters
        weight = 0.0;
        bias = 0.0;

        // Perform gradient descent
        for (int i = 0; i < epochs; i++)
        {
            // Compute predictions
            double[] predictions = new double[n];
            for (int j = 0; j < n; j++)
            {
                predictions[j] = Predict(features[j]);
            }

            // Compute gradients
            double dw = 0.0;
            double db = 0.0;
            for (int j = 0; j < n; j++)
            {
                dw += (predictions[j] - labels[j]) * features[j];
                db += (predictions[j] - labels[j]);
            }
            dw /= n;
            db /= n;

            // Update model parameters
            weight -= learningRate * dw;
            bias -= learningRate * db;
        }
    }

    public double Predict(double input)
    {
        double output = input * weight + bias;
        return output;
    }
}