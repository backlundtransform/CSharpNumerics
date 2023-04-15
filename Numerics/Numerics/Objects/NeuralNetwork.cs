using Numerics.Objects;
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


    private void ForwardPropagation(Tensor featureValues)
    {
        double[][] inputs = new double[featureValues.values.GetLength(0)][];
        for (int i = 0; i < featureValues.values.GetLength(0); i++)
        {
            inputs[i] = new double[1];
            inputs[i][0] = Convert.ToDouble(featureValues.values.GetValue(i));
        }

        for (int i = 0; i < hidden.Length; i++)
        {
            double sum = 0;
            for (int j = 0; j < inputs.Length; j++)
            {
                sum += inputs[j][0] * hiddenWeights[j][i];
            }
            hidden[i][0] = activationFunction(sum);
        }

  
        for (int i = 0; i < outputs.Length; i++)
        {
            double sum = 0;
            for (int j = 0; j < hidden.Length; j++)
            {
                sum += hidden[j][0] * outputWeights[j][i];
            }
            outputs[i][0] = activationFunction(sum);
        }
    }


    public double MeanAbsoluteError(double[][] predictedOutput, double[][] actualOutput)
    {
        int outputCount = predictedOutput.Length;
        double errorSum = 0;

        for (int i = 0; i < outputCount; i++)
        {
            double error = Math.Abs(predictedOutput[i][0] - actualOutput[i][0]);
            errorSum += error;
        }

        double meanError = errorSum / outputCount;
        return meanError;
    }


    private void BackwardPropagation(Tensor expectedValues)
    {

        var values = new double[expectedValues.values.GetLength(0)][];

        int outputCount = layers[2];
        // Calculate output error
        double[][] outputError = new double[outputCount][];
        for (int i = 0; i < outputCount; i++)
        {
            outputError[i] = new double[1];
            outputError[i][0] = values[i][0] - outputs[i][0];
        }

        // Calculate output gradient
        double[][] outputGradient = new double[outputCount][];
        for (int i = 0; i < outputCount; i++)
        {
            outputGradient[i] = new double[1];
            outputGradient[i][0] = -(2.0 / outputCount) * (values[i][0] - outputs[i][0]) * activationFunction(outputs[i][0]);
        }

        // Calculate hidden gradient
        int hiddenCount = layers[1];
        double[][] hiddenGradient = new double[hiddenCount][];
        for (int i = 0; i < hiddenCount; i++)
        {
            hiddenGradient[i] = new double[1];
            double sum = 0;
            for (int j = 0; j < outputCount; j++)
            {
                sum += outputGradient[j][0] * outputWeights[i][j];
            }
            hiddenGradient[i][0] = (1.0 / hiddenCount) * sum * activationFunction(hidden[i][0]);
        }

        // Update output weights
        for (int i = 0; i < hiddenCount; i++)
        {
            for (int j = 0; j < outputCount; j++)
            {
                outputWeights[i][j] -= learningRate * outputGradient[j][0] * hidden[i][0];
            }
        }

        // Update hidden weights
        int inputCount = layers[0];
        for (int i = 0; i < inputCount; i++)
        {
            for (int j = 0; j < hiddenCount; j++)
            {
                hiddenWeights[i][j] -= learningRate * hiddenGradient[j][0] * inputs[i][0];
            }
        }
    }

    public void Train(Tensor features, Tensor labels, int epochs)
    {
        for (var i = 0; i < epochs; i++)
        {
            ForwardPropagation(features);
            BackwardPropagation(labels);

        }
    }
    public double[] Predict(Tensor featureValues)
    {
        ForwardPropagation(featureValues);
        return outputs.Select(o => o[0]).ToArray();
    }
}