using System;
using System.Linq;

public class NeuralNetwork
{

    int[] layers;
    double learningRate;
    Func<double, double> activationFunction;
    private double[][] inputs;
    private double[][] hidden;
    private double[][] outputs;

    private double[][] hiddenWeights;
    private double[][] outputWeights;

    public NeuralNetwork(int[] layers, double learningRate, Func<double, double> activateFunction)
    {
        this.layers = layers;
        this.learningRate = learningRate;
        this.activationFunction = activateFunction;
        int inputCount = layers[0];
        int hiddenCount = layers[1];
        int outputCount = layers[2];
        inputs = new double[inputCount][];
        hidden = new double[hiddenCount][];
        outputs = new double[outputCount][];

        hiddenWeights = new double[inputCount][];
        outputWeights = new double[hiddenCount][];

        for (int i = 0; i < inputCount; i++)
        {
            inputs[i] = new double[1];
            hiddenWeights[i] = new double[hiddenCount];
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            hidden[i] = new double[1];
            outputWeights[i] = new double[outputCount];
        }

        for (int i = 0; i < outputCount; i++)
        {
            outputs[i] = new double[1];
        }

        InitializeWeights();
    }

    private void InitializeWeights()
    {
        var random = new Random();
        int inputCount = layers[0];
        int hiddenCount = layers[1];
        int outputCount = layers[2];

        for (int i = 0; i < inputCount; i++)
        {
            for (int j = 0; j < hiddenCount; j++)
            {
                hiddenWeights[i][j] = random.NextDouble();
            }
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            for (int j = 0; j < outputCount; j++)
            {
                outputWeights[i][j] = random.NextDouble();
            }
        }
    }



    private void ForwardPropagation(double[] featureValues)
    {

        int inputCount = layers[0];
        int hiddenCount = layers[1];
        int outputCount = layers[2];



        for (int i = 0; i < inputCount; i++)
        {
            inputs[i] = new double[1];
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            hidden[i] = new double[1];
        }

        for (int i = 0; i < outputCount; i++)
        {
            outputs[i] = new double[1];
        }

        for (int i = 0; i < inputCount; i++)
        {
            inputs[i][0] = featureValues[i];
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            hidden[i][0] = 0;

            for (int j = 0; j < inputCount; j++)
            {
                hidden[i][0] += inputs[j][0] * hiddenWeights[j][i];
            }

            hidden[i][0] = activationFunction(hidden[i][0]);
        }
        for (int i = 0; i < inputCount; i++)
        {
            inputs[i][0] = featureValues[i];
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            hidden[i][0] = 0;

            for (int j = 0; j < inputCount; j++)
            {
                hidden[i][0] += inputs[j][0] * hiddenWeights[j][i];
            }

            hidden[i][0] = activationFunction(hidden[i][0]);
        }

        for (int i = 0; i < outputCount; i++)
        {
            outputs[i][0] = 0;

            for (int j = 0; j < hiddenCount; j++)
            {
                outputs[i][0] += hidden[j][0] * outputWeights[j][i];
            }

            outputs[i][0] = activationFunction(outputs[i][0]);
        }
    }
    private void BackwardPropagation(double[] expectedValues)
    {
        int inputCount = layers[0];
        int hiddenCount = layers[1];
        int outputCount = layers[2];


        double error = 0;

        for (int i = 0; i < outputCount; i++)
        {
            var delta = expectedValues[i] - outputs[i][0];
            error += delta * delta;
        }
        error /= outputCount;
        double[] outputError = new double[outputCount];
        double[] hiddenError = new double[hiddenCount];

        for (int i = 0; i < outputCount; i++)
        {
            outputError[i] = outputs[i][0] * (1 - outputs[i][0]) * (expectedValues[i] - outputs[i][0]);
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            hiddenError[i] = hidden[i][0] * (1 - hidden[i][0]);

            for (int j = 0; j < outputCount; j++)
            {
                hiddenError[i] *= outputError[j] * outputWeights[i][j];
            }
        }

        for (int i = 0; i < hiddenCount; i++)
        {
            for (int j = 0; j < outputCount; j++)
            {
                outputWeights[i][j] += learningRate * outputError[j] * hidden[i][0];
            }
        }

        for (int i = 0; i < inputCount; i++)
        {
            for (int j = 0; j < hiddenCount; j++)
            {
                hiddenWeights[i][j] += learningRate * hiddenError[j] * inputs[i][0];
            }
        }
    }

    public void Train(double[] features, double[] labels, int epochs)
    {
     
            ForwardPropagation(features);
            BackwardPropagation(labels);

       

    }
    public double[] Predict(double[] featureValues)
    {
        ForwardPropagation(featureValues);
        return outputs.Select(o => o[0]).ToArray();
    }
}