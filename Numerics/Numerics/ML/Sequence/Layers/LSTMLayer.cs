using CSharpNumerics.ML.NeuralNetwork;
using CSharpNumerics.ML.NeuralNetwork.Layers.Interfaces;
using CSharpNumerics.Numerics.Objects;
using CSharpNumerics.Numerics.Optimization.Interfaces;
using System;

namespace CSharpNumerics.ML.Sequence.Layers;

/// <summary>
/// LSTM layer implementing full forward pass and BPTT (Backpropagation Through Time).
/// Gate equations per timestep t:
///   concat_t = [h_{t-1} ; x_t]
///   f_t = σ(W_f · concat_t + b_f)          forget gate
///   i_t = σ(W_i · concat_t + b_i)          input gate
///   c̃_t = tanh(W_c · concat_t + b_c)       candidate cell
///   o_t = σ(W_o · concat_t + b_o)          output gate
///   C_t = f_t ⊙ C_{t-1} + i_t ⊙ c̃_t
///   h_t = o_t ⊙ tanh(C_t)
/// </summary>
public class LSTMLayer : ILayer
{
    private readonly int _inputSize;
    private readonly int _hiddenSize;
    private readonly bool _returnSequences;
    private readonly double _clipNorm;

    // Gate weights: each (concatSize, hiddenSize) where concatSize = hiddenSize + inputSize
    private Matrix _Wf, _Wi, _Wc, _Wo;
    // Gate biases: each (hiddenSize)
    private VectorN _bf, _bi, _bc, _bo;

    // Gradient accumulators
    private Matrix _dWf, _dWi, _dWc, _dWo;
    private VectorN _dbf, _dbi, _dbc, _dbo;

    // Per-layer optimizer state
    private IOptimizer _wfOpt, _wiOpt, _wcOpt, _woOpt;
    private IOptimizer _bfOpt, _biOpt, _bcOpt, _boOpt;

    // Cached forward-pass values for BPTT (indexed by timestep)
    private VectorN[] _xCache;          // input at each timestep
    private VectorN[] _hCache;          // hidden state (h[0] = initial h)
    private VectorN[] _cCache;          // cell state   (c[0] = initial C)
    private VectorN[] _fCache;          // forget gate activations
    private VectorN[] _iCache;          // input gate activations
    private VectorN[] _candidateCache;  // candidate cell (c̃) values
    private VectorN[] _oCache;          // output gate activations
    private VectorN[] _tanhCCache;      // tanh(C_t) cached for backward
    private int _timeSteps;

    public LSTMLayer(int inputSize, int hiddenSize, bool returnSequences = true, double clipNorm = 5.0, int seed = 123)
    {
        if (inputSize <= 0) throw new ArgumentOutOfRangeException(nameof(inputSize));
        if (hiddenSize <= 0) throw new ArgumentOutOfRangeException(nameof(hiddenSize));

        _inputSize = inputSize;
        _hiddenSize = hiddenSize;
        _returnSequences = returnSequences;
        _clipNorm = clipNorm;

        int concatSize = hiddenSize + inputSize;
        var rng = new Random(seed);

        _Wf = CreateXavierMatrix(concatSize, hiddenSize, rng);
        _Wi = CreateXavierMatrix(concatSize, hiddenSize, rng);
        _Wc = CreateXavierMatrix(concatSize, hiddenSize, rng);
        _Wo = CreateXavierMatrix(concatSize, hiddenSize, rng);

        // Forget gate bias initialized to 1.0 (standard practice to reduce vanishing gradients)
        _bf = CreateConstantVector(hiddenSize, 1.0);
        _bi = new VectorN(hiddenSize);
        _bc = new VectorN(hiddenSize);
        _bo = new VectorN(hiddenSize);

        ResetGradients();
    }

    public int InputSize => _inputSize;

    public int HiddenSize => _hiddenSize;

    public int ParameterCount
    {
        get
        {
            int concatSize = _hiddenSize + _inputSize;
            return 4 * ((concatSize * _hiddenSize) + _hiddenSize);
        }
    }

    public VectorN[] Forward(VectorN[] input, bool training = true)
    {
        if (input is null || input.Length == 0)
            throw new ArgumentException("Input sequence must not be empty.", nameof(input));

        _timeSteps = input.Length;
        int T = _timeSteps;

        // Allocate caches: index 0 = initial state, 1..T = timestep outputs
        _xCache = new VectorN[T];
        _hCache = new VectorN[T + 1];
        _cCache = new VectorN[T + 1];
        _fCache = new VectorN[T];
        _iCache = new VectorN[T];
        _candidateCache = new VectorN[T];
        _oCache = new VectorN[T];
        _tanhCCache = new VectorN[T];

        // Initialize h_0 and C_0 to zeros
        _hCache[0] = new VectorN(_hiddenSize);
        _cCache[0] = new VectorN(_hiddenSize);

        var output = _returnSequences ? new VectorN[T] : new VectorN[1];

        for (int t = 0; t < T; t++)
        {
            if (input[t].Length != _inputSize)
                throw new ArgumentException($"Input at timestep {t} has width {input[t].Length}, expected {_inputSize}.");

            _xCache[t] = new VectorN(input[t].Values);

            // concat_t = [h_{t-1} ; x_t]
            var concat = _hCache[t].Concat(input[t]);

            // Gate computations
            var fRaw = (_Wf.Transpose() * concat) + _bf;
            var iRaw = (_Wi.Transpose() * concat) + _bi;
            var cRaw = (_Wc.Transpose() * concat) + _bc;
            var oRaw = (_Wo.Transpose() * concat) + _bo;

            var ft = Sigmoid(fRaw);
            var it = Sigmoid(iRaw);
            var ct = Tanh(cRaw);
            var ot = Sigmoid(oRaw);

            // Cell and hidden state
            var Ct = ft.Hadamard(_cCache[t]) + it.Hadamard(ct);
            var tanhCt = Tanh(Ct);
            var ht = ot.Hadamard(tanhCt);

            // Cache for backward
            _fCache[t] = ft;
            _iCache[t] = it;
            _candidateCache[t] = ct;
            _oCache[t] = ot;
            _tanhCCache[t] = tanhCt;
            _hCache[t + 1] = ht;
            _cCache[t + 1] = Ct;

            if (_returnSequences)
                output[t] = ht;
        }

        if (!_returnSequences)
            output[0] = _hCache[T];

        return output;
    }

    public VectorN[] Backward(VectorN[] gradOutput)
    {
        if (_xCache is null)
            throw new InvalidOperationException("Call Forward before Backward.");

        int T = _timeSteps;

        // If returnSequences=false, gradOutput has length 1 — assign it to the last timestep
        var gradH = new VectorN[T];
        if (_returnSequences)
        {
            if (gradOutput.Length != T)
                throw new ArgumentException("Gradient output length must match sequence length when returnSequences=true.");
            for (int t = 0; t < T; t++)
                gradH[t] = gradOutput[t];
        }
        else
        {
            if (gradOutput.Length != 1)
                throw new ArgumentException("Gradient output must have length 1 when returnSequences=false.");
            for (int t = 0; t < T - 1; t++)
                gradH[t] = new VectorN(_hiddenSize);
            gradH[T - 1] = gradOutput[0];
        }

        var gradInput = new VectorN[T];
        var dh_next = new VectorN(_hiddenSize); // gradient flowing from future timestep into h
        var dC_next = new VectorN(_hiddenSize); // gradient flowing from future timestep into C

        int concatSize = _hiddenSize + _inputSize;

        for (int t = T - 1; t >= 0; t--)
        {
            var dh = gradH[t] + dh_next; // combine external gradient with recurrent gradient

            // d(loss)/d(o_t) = dh ⊙ tanh(C_t)
            var dOt = dh.Hadamard(_tanhCCache[t]);

            // d(loss)/d(C_t) = dh ⊙ o_t ⊙ (1 - tanh²(C_t)) + dC_next
            var dtanhC = dh.Hadamard(_oCache[t]);
            var tanhCDeriv = TanhDerivative(_tanhCCache[t]);
            var dCt = dtanhC.Hadamard(tanhCDeriv) + dC_next;

            // d(loss)/d(f_t) = dC_t ⊙ C_{t-1}
            var dFt = dCt.Hadamard(_cCache[t]);

            // d(loss)/d(i_t) = dC_t ⊙ c̃_t
            var dIt = dCt.Hadamard(_candidateCache[t]);

            // d(loss)/d(c̃_t) = dC_t ⊙ i_t
            var dCandidateT = dCt.Hadamard(_iCache[t]);

            // Derivatives through gate activations (σ and tanh)
            var dFtRaw = dFt.Hadamard(SigmoidDerivative(_fCache[t]));
            var dItRaw = dIt.Hadamard(SigmoidDerivative(_iCache[t]));
            var dCandidateRaw = dCandidateT.Hadamard(TanhDerivative(_candidateCache[t]));
            var dOtRaw = dOt.Hadamard(SigmoidDerivative(_oCache[t]));

            // Concat input for this timestep
            var concat = _hCache[t].Concat(_xCache[t]);

            // Accumulate weight gradients
            _dWf += concat.Outer(dFtRaw);
            _dWi += concat.Outer(dItRaw);
            _dWc += concat.Outer(dCandidateRaw);
            _dWo += concat.Outer(dOtRaw);

            _dbf += dFtRaw;
            _dbi += dItRaw;
            _dbc += dCandidateRaw;
            _dbo += dOtRaw;

            // Gradient w.r.t. concat = [dh_prev; dx_t]
            var dConcat = (_Wf * dFtRaw) + (_Wi * dItRaw) + (_Wc * dCandidateRaw) + (_Wo * dOtRaw);

            // Split dConcat into dh_prev and dx_t
            dh_next = SliceVector(dConcat, 0, _hiddenSize);
            gradInput[t] = SliceVector(dConcat, _hiddenSize, _inputSize);

            // dC flowing backward through forget gate
            dC_next = dCt.Hadamard(_fCache[t]);
        }

        // Apply gradient clipping to accumulated gradients
        if (_clipNorm > 0)
            ClipGradients();

        return gradInput;
    }

    public void ApplyGradients(IOptimizer weightOptimizer, IOptimizer biasOptimizer, int batchSize)
    {
        if (batchSize <= 0)
            throw new ArgumentOutOfRangeException(nameof(batchSize));

        _wfOpt ??= OptimizerFactory.Clone(weightOptimizer);
        _wiOpt ??= OptimizerFactory.Clone(weightOptimizer);
        _wcOpt ??= OptimizerFactory.Clone(weightOptimizer);
        _woOpt ??= OptimizerFactory.Clone(weightOptimizer);
        _bfOpt ??= OptimizerFactory.Clone(biasOptimizer);
        _biOpt ??= OptimizerFactory.Clone(biasOptimizer);
        _bcOpt ??= OptimizerFactory.Clone(biasOptimizer);
        _boOpt ??= OptimizerFactory.Clone(biasOptimizer);

        _Wf = ApplyMatrixGradient(_Wf, _dWf, _wfOpt, batchSize);
        _Wi = ApplyMatrixGradient(_Wi, _dWi, _wiOpt, batchSize);
        _Wc = ApplyMatrixGradient(_Wc, _dWc, _wcOpt, batchSize);
        _Wo = ApplyMatrixGradient(_Wo, _dWo, _woOpt, batchSize);

        _bf = _bfOpt.Step(_bf, _dbf / batchSize);
        _bi = _biOpt.Step(_bi, _dbi / batchSize);
        _bc = _bcOpt.Step(_bc, _dbc / batchSize);
        _bo = _boOpt.Step(_bo, _dbo / batchSize);

        ResetGradients();
    }

    private Matrix ApplyMatrixGradient(Matrix weights, Matrix gradients, IOptimizer optimizer, int batchSize)
    {
        int rows = weights.rowLength;
        int cols = weights.columnLength;
        var flatW = new double[rows * cols];
        var flatG = new double[rows * cols];

        for (int r = 0; r < rows; r++)
        {
            for (int c = 0; c < cols; c++)
            {
                int idx = (r * cols) + c;
                flatW[idx] = weights.values[r, c];
                flatG[idx] = gradients.values[r, c] / batchSize;
            }
        }

        var updated = optimizer.Step(new VectorN(flatW), new VectorN(flatG));
        var result = new double[rows, cols];
        for (int r = 0; r < rows; r++)
        {
            for (int c = 0; c < cols; c++)
            {
                result[r, c] = updated[(r * cols) + c];
            }
        }

        return new Matrix(result);
    }

    private void ClipGradients()
    {
        double totalNorm = 0.0;
        totalNorm += MatrixSquaredNorm(_dWf);
        totalNorm += MatrixSquaredNorm(_dWi);
        totalNorm += MatrixSquaredNorm(_dWc);
        totalNorm += MatrixSquaredNorm(_dWo);
        totalNorm += VectorSquaredNorm(_dbf);
        totalNorm += VectorSquaredNorm(_dbi);
        totalNorm += VectorSquaredNorm(_dbc);
        totalNorm += VectorSquaredNorm(_dbo);
        totalNorm = Math.Sqrt(totalNorm);

        if (totalNorm > _clipNorm)
        {
            double scale = _clipNorm / totalNorm;
            _dWf = ScaleMatrix(_dWf, scale);
            _dWi = ScaleMatrix(_dWi, scale);
            _dWc = ScaleMatrix(_dWc, scale);
            _dWo = ScaleMatrix(_dWo, scale);
            _dbf = ScaleVector(_dbf, scale);
            _dbi = ScaleVector(_dbi, scale);
            _dbc = ScaleVector(_dbc, scale);
            _dbo = ScaleVector(_dbo, scale);
        }
    }

    private void ResetGradients()
    {
        int concatSize = _hiddenSize + _inputSize;
        _dWf = new Matrix(concatSize, _hiddenSize);
        _dWi = new Matrix(concatSize, _hiddenSize);
        _dWc = new Matrix(concatSize, _hiddenSize);
        _dWo = new Matrix(concatSize, _hiddenSize);
        _dbf = new VectorN(_hiddenSize);
        _dbi = new VectorN(_hiddenSize);
        _dbc = new VectorN(_hiddenSize);
        _dbo = new VectorN(_hiddenSize);
    }

    private static VectorN Sigmoid(VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = 1.0 / (1.0 + Math.Exp(-v[i]));
        return new VectorN(result);
    }

    private static VectorN SigmoidDerivative(VectorN activated)
    {
        // σ'(x) = σ(x)·(1 - σ(x)), where activated = σ(x)
        var result = new double[activated.Length];
        for (int i = 0; i < activated.Length; i++)
            result[i] = activated[i] * (1.0 - activated[i]);
        return new VectorN(result);
    }

    private static VectorN Tanh(VectorN v)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = Math.Tanh(v[i]);
        return new VectorN(result);
    }

    private static VectorN TanhDerivative(VectorN activated)
    {
        // tanh'(x) = 1 - tanh²(x), where activated = tanh(x)
        var result = new double[activated.Length];
        for (int i = 0; i < activated.Length; i++)
            result[i] = 1.0 - (activated[i] * activated[i]);
        return new VectorN(result);
    }

    private static VectorN SliceVector(VectorN v, int start, int length)
    {
        var result = new double[length];
        for (int i = 0; i < length; i++)
            result[i] = v[start + i];
        return new VectorN(result);
    }

    private static VectorN CreateConstantVector(int size, double value)
    {
        var values = new double[size];
        for (int i = 0; i < size; i++)
            values[i] = value;
        return new VectorN(values);
    }

    private static Matrix CreateXavierMatrix(int rows, int cols, Random rng)
    {
        var values = new double[rows, cols];
        double limit = Math.Sqrt(6.0 / (rows + cols));
        for (int r = 0; r < rows; r++)
        {
            for (int c = 0; c < cols; c++)
            {
                values[r, c] = (rng.NextDouble() * 2.0 - 1.0) * limit;
            }
        }

        return new Matrix(values);
    }

    private static double MatrixSquaredNorm(Matrix m)
    {
        double sum = 0.0;
        for (int r = 0; r < m.rowLength; r++)
        {
            for (int c = 0; c < m.columnLength; c++)
            {
                double v = m.values[r, c];
                sum += v * v;
            }
        }

        return sum;
    }

    private static double VectorSquaredNorm(VectorN v)
    {
        double sum = 0.0;
        for (int i = 0; i < v.Length; i++)
            sum += v[i] * v[i];
        return sum;
    }

    private static Matrix ScaleMatrix(Matrix m, double scale)
    {
        var result = new double[m.rowLength, m.columnLength];
        for (int r = 0; r < m.rowLength; r++)
        {
            for (int c = 0; c < m.columnLength; c++)
            {
                result[r, c] = m.values[r, c] * scale;
            }
        }

        return new Matrix(result);
    }

    private static VectorN ScaleVector(VectorN v, double scale)
    {
        var result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = v[i] * scale;
        return new VectorN(result);
    }
}
