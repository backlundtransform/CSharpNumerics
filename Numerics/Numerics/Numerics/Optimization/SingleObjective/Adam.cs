using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.Numerics.Optimization.SingleObjective;

/// <summary>
/// Adam optimiser (Kingma &amp; Ba, 2015).
/// Maintains per-parameter first and second moment estimates with bias correction.
/// Also supports AdamW (decoupled weight decay).
/// </summary>
public class Adam : IOptimizer
{
    public double LearningRate { get; set; }
    public double Beta1 { get; set; }
    public double Beta2 { get; set; }
    public double Epsilon { get; set; }
    public double WeightDecay { get; set; }
    public bool DecoupledWeightDecay { get; set; }

    private VectorN _m; // first moment
    private VectorN _v; // second moment
    private int _t;     // timestep
    private bool _initialised;

    public Adam(double learningRate = 0.001, double beta1 = 0.9, double beta2 = 0.999,
        double epsilon = 1e-8, double weightDecay = 0.0, bool decoupledWeightDecay = false)
    {
        LearningRate = learningRate;
        Beta1 = beta1;
        Beta2 = beta2;
        Epsilon = epsilon;
        WeightDecay = weightDecay;
        DecoupledWeightDecay = decoupledWeightDecay;
    }

    public VectorN Step(VectorN parameters, VectorN gradient)
    {
        int n = parameters.Length;

        if (!_initialised || _m.Length != n)
        {
            _m = new VectorN(n);
            _v = new VectorN(n);
            _t = 0;
            _initialised = true;
        }

        _t++;

        // Optional L2 added to gradient (classic Adam) vs decoupled (AdamW)
        VectorN g = (!DecoupledWeightDecay && WeightDecay > 0)
            ? gradient + WeightDecay * parameters
            : gradient;

        // Update biased first/second moment estimates
        for (int i = 0; i < n; i++)
        {
            _m[i] = Beta1 * _m[i] + (1 - Beta1) * g[i];
            _v[i] = Beta2 * _v[i] + (1 - Beta2) * g[i] * g[i];
        }

        // Bias correction
        double bc1 = 1.0 - Math.Pow(Beta1, _t);
        double bc2 = 1.0 - Math.Pow(Beta2, _t);

        var result = new double[n];
        for (int i = 0; i < n; i++)
        {
            double mHat = _m[i] / bc1;
            double vHat = _v[i] / bc2;
            result[i] = parameters[i] - LearningRate * mHat / (Math.Sqrt(vHat) + Epsilon);
        }

        // AdamW decoupled weight decay
        if (DecoupledWeightDecay && WeightDecay > 0)
        {
            for (int i = 0; i < n; i++)
                result[i] -= LearningRate * WeightDecay * parameters[i];
        }

        return new VectorN(result);
    }

    public void Reset()
    {
        _initialised = false;
        _m = default;
        _v = default;
        _t = 0;
    }
}
