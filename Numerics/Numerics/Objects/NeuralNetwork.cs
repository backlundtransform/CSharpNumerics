using System;

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
        var n = features.Length;

       
        weight = 0.0;
        bias = 0.0;

     
        for (int i = 0; i < epochs; i++)
        {
           
            var predictions = new double[n];
            for (int j = 0; j < n; j++)
            {
                predictions[j] = Predict(features[j]);
            }

            var dw = 0.0;
            var db = 0.0;
            for (var j = 0; j < n; j++)
            {
                dw += (predictions[j] - labels[j]) * features[j];
                db += (predictions[j] - labels[j]);
            }
            dw /= n;
            db /= n;

            weight -= learningRate * dw;
            bias -= learningRate * db;
        }
    }


    public double Predict(double input)
    {
        var output = input * weight + bias;
        return output;
    }
}